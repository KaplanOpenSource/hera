"""
Native pytest tests for the HighFreqToolKit, analysis calculators, and turbulence statistics.

Data is loaded into the test project via ``conftest.hera_test_project``.
The ``hf_toolkit`` fixture (session-scoped, from conftest) provides
a real ``HighFreqToolKit`` backed by MongoDB – no direct file-path access.

Replaces:
  - hera/measurements/meteorology/highfreqdata/test_unit_highfreq_toolkit.py
  - hera/tests/json_definitions/test_definitions_highfreq.json
"""

import os

import dask.dataframe as dd
import pandas as pd
import pytest

from hera.measurements.meteorology.highfreqdata.analysis.abstractcalculator import AbstractCalculator
from hera.measurements.meteorology.highfreqdata.analysis.meandatacalculator import MeanDataCalculator, AveragingCalculator
from hera.measurements.meteorology.highfreqdata.analysis.analysislayer import RawdataAnalysis
from hera.measurements.meteorology.highfreqdata.analysis.turbulencestatistics import singlePointTurbulenceStatistics


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def sonic_df(hf_toolkit):
    """Computed Pandas DataFrame from the sonic datasource in the project."""
    ds = hf_toolkit.getDataSourceData("slicedYamim_sonic")
    if ds is None:
        pytest.skip("slicedYamim_sonic datasource not loaded in project")
    # dask → pandas
    return ds.compute() if hasattr(ds, "compute") else ds


@pytest.fixture(scope="module")
def trh_df(hf_toolkit):
    """Computed Pandas DataFrame from the TRH datasource in the project."""
    ds = hf_toolkit.getDataSourceData("slicedYamim_TRH")
    if ds is None:
        pytest.skip("slicedYamim_TRH datasource not loaded in project")
    return ds.compute() if hasattr(ds, "compute") else ds


@pytest.fixture(scope="module")
def sonic_df_indexed(sonic_df):
    """sonic_df with a DatetimeIndex on the 'Time' column."""
    df = sonic_df.copy()
    # If the index is already a DatetimeIndex (parquet stores Time as index)
    if isinstance(df.index, pd.DatetimeIndex):
        if df.index.name is None:
            df.index.name = "Time"
        return df.sort_index()
    if "Time" in df.columns:
        df["Time"] = pd.to_datetime(df["Time"], errors="coerce")
        df = df.dropna(subset=["Time"]).set_index("Time").sort_index()
    elif "timestamp" in df.columns:
        df["timestamp"] = pd.to_datetime(df["timestamp"], errors="coerce")
        df = df.dropna(subset=["timestamp"]).set_index("timestamp").sort_index()
    else:
        pytest.fail("No Time or timestamp column in sonic data")
    return df


# ---------------------------------------------------------------------------
# Toolkit property tests
# ---------------------------------------------------------------------------

class TestHighFreqToolkitInit:
    def test_docType_property(self, hf_toolkit):
        assert hf_toolkit.docType == "highFreqMeteorology_HighFreqData"


# ---------------------------------------------------------------------------
# Data reading tests
# ---------------------------------------------------------------------------

class TestReadData:
    def test_read_sonic_data(self, sonic_df):
        assert isinstance(sonic_df, pd.DataFrame)
        assert not sonic_df.columns.empty

    def test_read_trh_data(self, trh_df):
        assert isinstance(trh_df, pd.DataFrame)
        assert not trh_df.columns.empty

    def test_read_nonexistent_datasource(self, hf_toolkit):
        result = hf_toolkit.getDataSourceData("totally_nonexistent_source")
        assert result is None


class TestTimeRange:
    def _get_time_series(self, df_pd):
        """Return the time series — either from the index or a column."""
        if isinstance(df_pd.index, pd.DatetimeIndex):
            return df_pd.index
        for col in ("timestamp", "Time"):
            if col in df_pd.columns:
                return df_pd[col]
        return None

    def test_sonic_time_range(self, sonic_df):
        ts = self._get_time_series(sonic_df)
        assert ts is not None, "No time column or datetime index found in sonic data"
        assert ts.min() < ts.max()

    def test_trh_time_range(self, trh_df):
        ts = self._get_time_series(trh_df)
        assert ts is not None, "No time column or datetime index found in TRH data"
        assert ts.min() < ts.max()


class TestSpecificPoints:
    def test_sonic_first_row(self, sonic_df):
        first = sonic_df.iloc[0]
        assert abs(first["u"] - (-2.65)) < 0.01
        assert abs(first["v"] - 3.05) < 0.01
        assert abs(first["w"] - 1.96) < 0.01
        assert abs(first["T"] - 24.93) < 0.01

    def test_trh_first_row(self, trh_df):
        first = trh_df.iloc[0]
        assert abs(first["TC_T"] - 24.71) < 0.01
        assert abs(first["RH"] - 70.2) < 0.01


class TestErrorPaths:
    def test_campbelToParquet_nonexistent(self, hf_toolkit):
        with pytest.raises(ValueError):
            hf_toolkit.campbelToParquet(binaryFile="/path/to/nonexistent_file.dat")

    def test_asciiToParquet_nonexistent(self, hf_toolkit):
        with pytest.raises(FileNotFoundError):
            hf_toolkit.asciiToParquet(path="/path/to/nonexistent_file.txt")


# ---------------------------------------------------------------------------
# AbstractCalculator tests
# ---------------------------------------------------------------------------

class TestAbstractCalculator:
    def test_init_basic(self, sonic_df):
        ac = AbstractCalculator(
            rawData=sonic_df, metadata={"samplingWindow": "10s", "isMissingData": False}
        )
        assert ac.RawData is not None
        assert ac.metaData["samplingWindow"] == "10s"

    def test_sampling_window(self, sonic_df):
        ac = AbstractCalculator(
            rawData=sonic_df, metadata={"samplingWindow": "10s", "isMissingData": False}
        )
        assert ac.SamplingWindow == "10s"

    def test_compute_methods_exist(self, sonic_df):
        ac = AbstractCalculator(
            rawData=sonic_df, metadata={"samplingWindow": "10s", "isMissingData": False}
        )
        assert hasattr(ac, "compute")
        assert hasattr(ac, "_compute")

    def test_set_save_properties(self, sonic_df):
        ac = AbstractCalculator(rawData=sonic_df, metadata={"samplingWindow": "10s"})
        props = ac.set_saveProperties(dataFormat="testFormat")
        assert props is None or isinstance(props, dict)


# ---------------------------------------------------------------------------
# MeanDataCalculator tests
# ---------------------------------------------------------------------------

class TestMeanDataCalculator:
    @pytest.fixture()
    def calc(self, sonic_df_indexed):
        df = sonic_df_indexed.copy()
        df["TC_T_bar"] = df["T"]
        metadata = {
            "samplingWindow": "10s",
            "isMissingData": False,
            "start": df.index.min(),
            "end": df.index.max(),
        }
        turb = singlePointTurbulenceStatistics(df, metadata=metadata)
        return MeanDataCalculator(TurbCalcOrData=turb)

    def test_calculate_mean(self, calc):
        assert calc.MeanData is not None
        assert not calc.MeanData.empty

    def test_hour_and_timeWithinDay(self, calc):
        calc.hour()
        calc.timeWithinDay()
        assert "hour" in calc.MeanData.columns
        assert "timeWithinDay" in calc.MeanData.columns

    def test_horizontalSpeed(self, calc):
        calc.horizontalSpeed()
        assert "horizontal_speed_bar" in calc.MeanData.columns

    def test_sigma_sigmaH(self, calc):
        calc.sigma()
        calc.sigmaH()
        for col in ["sigmaU", "sigmaV", "sigmaW", "sigmaH"]:
            assert col in calc.MeanData.columns

    def test_Ustar_and_uStarOverWindSpeed(self, calc):
        calc.Ustar()
        calc.uStarOverWindSpeed()
        assert "Ustar" in calc.MeanData.columns
        assert "uStarOverWindSpeed" in calc.MeanData.columns

    def test_compute_returns_dataframe(self, calc):
        result = calc.compute()
        assert isinstance(result, pd.DataFrame)


class TestMeanDataCalculatorAdvanced:
    def test_TKE(self, sonic_df_indexed):
        df = sonic_df_indexed.head(100).copy()
        metadata = {
            "samplingWindow": "10s",
            "isMissingData": False,
            "start": df.index.min(),
            "end": df.index.max(),
        }
        turb = singlePointTurbulenceStatistics(df, metadata)
        turb.secondMoments()
        calc = MeanDataCalculator(TurbCalcOrData=turb, metadata=metadata)
        calc.TKE()
        assert "TKE" in calc.MeanData.columns

    def test_MOLength(self, sonic_df_indexed):
        df = sonic_df_indexed.head(100).copy()
        metadata = {
            "samplingWindow": "10s",
            "isMissingData": False,
            "start": df.index.min(),
            "end": df.index.max(),
        }
        turb = singlePointTurbulenceStatistics(df, metadata)
        turb.secondMoments()
        calc = MeanDataCalculator(TurbCalcOrData=turb, metadata=metadata)
        calc.MeanData["TC_T_bar"] = 20
        calc.MeanData["Ustar"] = 0.3
        calc.MOLength()
        assert "L" in calc.MeanData.columns


# ---------------------------------------------------------------------------
# RawdataAnalysis tests
# ---------------------------------------------------------------------------

class TestRawdataAnalysis:
    def test_singlePointTurbulenceStatistics_returns_instance(self, hf_toolkit, sonic_df_indexed):
        df = sonic_df_indexed.copy()
        analysis = RawdataAnalysis(datalayer=hf_toolkit)
        result = analysis.singlePointTurbulenceStatistics(
            sonicData=df,
            samplingWindow="10s",
            start=df.index.min(),
            end=df.index.max(),
            height=10,
            buildingHeight=5,
            averagedHeight=7,
            inmemory=True,
            isMissingData=False,
        )
        assert result is not None

    def test_singlePointTurbulenceStatistics_raises_on_invalid(self, hf_toolkit):
        analysis = RawdataAnalysis(datalayer=hf_toolkit)
        with pytest.raises(TypeError):
            analysis.singlePointTurbulenceStatistics(sonicData="not_a_dataframe")

    def test_AveragingCalculator(self, hf_toolkit, sonic_df_indexed):
        df = sonic_df_indexed.copy()
        analysis = RawdataAnalysis(datalayer=hf_toolkit)
        calc = analysis.AveragingCalculator(
            df, "5min", df.index[0], df.index[-1], 2, 10, 8,
        )
        result = calc.compute()
        assert isinstance(result, pd.DataFrame)
        assert not result.empty

    def test_AveragingCalculator_raises_on_invalid(self, hf_toolkit):
        analysis = RawdataAnalysis(datalayer=hf_toolkit)
        with pytest.raises(TypeError):
            analysis.AveragingCalculator(sonicData=None)


# ---------------------------------------------------------------------------
# singlePointTurbulenceStatistics unit tests (synthetic data)
# ---------------------------------------------------------------------------

class TestSinglePointTurbulenceStatistics:
    @pytest.fixture()
    def turb_calc(self):
        df = pd.DataFrame({
            "u": [1.0, 2.0, 3.0],
            "v": [0.5, 1.5, 2.5],
            "w": [0.1, 0.2, 0.3],
            "T": [20.0, 21.0, 22.0],
        })
        metadata = {
            "samplingWindow": "10s",
            "start": pd.Timestamp("2020-01-01 00:00:00"),
            "end": pd.Timestamp("2020-01-01 00:00:02"),
            "height": 10,
            "buildingHeight": 2,
            "averagedHeight": 6,
            "isMissingData": False,
        }
        df["Time"] = pd.date_range(start=metadata["start"], periods=len(df), freq="1s")
        df = df.set_index("Time")
        return singlePointTurbulenceStatistics(df, metadata)

    def test_instantiation(self, turb_calc):
        assert turb_calc is not None
        assert hasattr(turb_calc, "RawData")
        assert hasattr(turb_calc, "metaData")

    def test_invalid_input_type(self):
        with pytest.raises(ValueError):
            singlePointTurbulenceStatistics("not_a_dataframe", {"samplingWindow": "10s"})

    def test_fluctuations(self, turb_calc):
        turb_calc.fluctuations()
        for col in ["u_bar", "v_bar", "w_bar", "T_bar", "up", "vp", "wp", "Tp"]:
            assert col in turb_calc.RawData.columns

    def test_secondMoments(self, turb_calc):
        turb_calc.secondMoments()
        for col in ["uu", "vv", "ww", "uw", "uv", "uT", "vT", "wT", "TT", "vw"]:
            assert col in turb_calc.TemporaryData.columns

    def test_sigma(self, turb_calc):
        turb_calc.sigma()
        for col in ["sigmaU", "sigmaV", "sigmaW"]:
            assert col in turb_calc.TemporaryData.columns

    def test_horizontalSpeed(self, turb_calc):
        turb_calc.horizontalSpeed()
        assert "horizontal_speed_bar" in turb_calc.TemporaryData.columns

    def test_Ustar(self, turb_calc):
        turb_calc.Ustar()
        assert "Ustar" in turb_calc.TemporaryData.columns

    def test_TKE(self, turb_calc):
        turb_calc.uu().vv().ww().TKE()
        assert "TKE" in turb_calc.TemporaryData.columns

    def test_MOLength_Sonic(self, turb_calc):
        turb_calc.wT().Ustar().MOLength_Sonic()
        assert "L_Sonic" in turb_calc.TemporaryData.columns
