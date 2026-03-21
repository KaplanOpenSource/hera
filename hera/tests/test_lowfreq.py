"""
Native pytest tests for the lowFreqToolKit and its analysis/presentation layers.

Data is loaded into the test project via ``conftest.hera_test_project``.
The ``lf_toolkit`` fixture (session-scoped, from conftest) provides
a real ``lowFreqToolKit`` backed by MongoDB – no direct file-path access.

Replaces:
  - hera/measurements/meteorology/lowfreqdata/test_unit_lowfreq_toolkit.py
  - hera/tests/json_definitions/lowfreq_toolkit.json
"""

import os
import tempfile

import numpy as np
import pandas as pd
import pytest

from hera.measurements.meteorology.lowfreqdata.presentationLayer import presenation
from hera.utils.statistics import calcDist2d


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def lowfreq_df(lf_toolkit):
    """Load the YAVNEEL dataframe through the project's datasource.

    The ``parquet`` data handler returns a dask DataFrame; we ``.compute()``
    to get a pandas DataFrame and fix the datetime column.
    """
    ds = lf_toolkit.getDataSourceData("YAVNEEL")
    if ds is None:
        pytest.skip("YAVNEEL datasource not loaded in project")

    # dask → pandas
    df = ds.compute() if hasattr(ds, "compute") else ds

    if "datetime" not in df.columns:
        pytest.fail("'datetime' column missing in dataframe")
    df["datetime"] = pd.to_datetime(df["datetime"], utc=True)
    return df


# ---------------------------------------------------------------------------
# Toolkit structure tests
# ---------------------------------------------------------------------------

class TestLowFreqToolkitInit:
    def test_has_analysis(self, lf_toolkit):
        assert lf_toolkit.analysis is not None

    def test_has_presentation(self, lf_toolkit):
        assert lf_toolkit.presentation is not None

    def test_has_docType(self, lf_toolkit):
        assert hasattr(lf_toolkit, "docType")

    def test_docType_value(self, lf_toolkit):
        doc_type = lf_toolkit.docType
        assert isinstance(doc_type, str)
        assert "LOWFREQ" in doc_type.upper()


# ---------------------------------------------------------------------------
# Analysis tests
# ---------------------------------------------------------------------------

class TestAnalysisAddDatesColumns:
    def test_basic(self, lf_toolkit, lowfreq_df):
        df = lowfreq_df.copy()
        enriched = lf_toolkit.analysis.addDatesColumns(df, datecolumn="datetime")
        assert enriched is not None
        assert len(enriched) > 0


class TestAnalysisCalcHourlyDist:
    def test_max_normalized(self, lf_toolkit, lowfreq_df):
        df = lowfreq_df.copy()
        x, y, M = lf_toolkit.analysis.calcHourlyDist(
            df, Field="RH", normalization="max_normalized"
        )
        assert x is not None
        assert y is not None
        assert M is not None
        assert len(y) == M.shape[0]
        assert len(x) == M.shape[1]

    def test_density(self, lf_toolkit, lowfreq_df):
        df = lowfreq_df.copy()
        x, y, M = lf_toolkit.analysis.calcHourlyDist(
            df, Field="RH", normalization="density"
        )
        assert (M >= 0).all(), "Density matrix should not contain negative values"


class TestCalcDist2dYNormalized:
    def test_y_normalized_behaviour(self):
        x = np.array([0.5, 1.5, 2.5, 1.5])
        y = np.array([1, 1, 1, 2])
        x_range = (0, 3)
        y_range = (0, 3)

        x_mid, y_mid, M = calcDist2d(
            x, y, bins=3, normalization="y_normalized",
            x_range=x_range, y_range=y_range,
        )

        # Row 0: no data -> all zeros
        assert np.allclose(M[0], 0)
        # Row 1: multiple values -> sum to 1
        assert np.isclose(M[1].sum(), 1.0)
        assert np.count_nonzero(M[1]) > 1
        # Row 2: single value -> sum to 1 with one non-zero
        assert np.isclose(M[2].sum(), 1.0)
        assert np.count_nonzero(M[2]) == 1


class TestResampleSecondMoments:
    def test_basic(self, lf_toolkit, lowfreq_df):
        df = lowfreq_df.copy().set_index("datetime")
        df["X"] = df["WS"]
        df["Y"] = df["RH"]
        df["XX"] = df["X"] ** 2
        df["YY"] = df["Y"] ** 2
        df["XY"] = df["X"] * df["Y"]
        df["X_bar"] = df["X"]
        df["Y_bar"] = df["Y"]

        result = lf_toolkit.analysis.resampleSecondMoments(
            df,
            SamplingWindow="D",
            fieldsFirstMoments=["X_bar", "Y_bar"],
            fieldsSecondMoments=["X", "Y"],
        )
        assert isinstance(result, pd.DataFrame)
        for col in ["XX", "XY", "YY"]:
            assert col in result.columns


# ---------------------------------------------------------------------------
# Presentation tests
# ---------------------------------------------------------------------------

class TestPresentationDailyPlots:
    def test_plotScatter(self, lf_toolkit, lowfreq_df):
        df = lowfreq_df.copy()
        ax = lf_toolkit.presentation.dailyPlots.plotScatter(df, plotField="RH")
        assert ax is not None

    def test_dateLinePlot(self, lf_toolkit, lowfreq_df):
        df = lowfreq_df.copy()
        example_date = df["datetime"].dt.date.astype(str).iloc[0]
        ax, line = lf_toolkit.presentation.dailyPlots.dateLinePlot(
            df, plotField="RH", date=example_date
        )
        assert ax is not None
        assert len(line) > 0

    def test_plotProbContourf(self, lf_toolkit, lowfreq_df):
        df = lowfreq_df.copy()
        CS, CFS, ax = lf_toolkit.presentation.dailyPlots.plotProbContourf(
            df, plotField="RH"
        )
        assert CS is not None
        assert CFS is not None
        assert ax is not None


class TestPresentationDataMatchPlots:
    def test_dateLinePlot_matches_data(self, lf_toolkit, lowfreq_df):
        df = lowfreq_df.copy()
        df["datetime"] = pd.to_datetime(df["datetime"], utc=True)
        df = df.set_index("datetime")

        example_date = df.index.date[0].strftime("%Y-%m-%d")
        ax, line = lf_toolkit.presentation.dailyPlots.dateLinePlot(
            df, plotField="RH", date=example_date
        )

        line_obj = line[0]
        y_vals = line_obj.get_ydata()

        daily = df.copy()
        daily["RH"] = daily["RH"].where(daily["RH"] > -5000)
        daily = daily.assign(houronly=daily.index.hour + daily.index.minute / 60.0)
        filtered = daily[daily.index.date.astype(str) == example_date].dropna(subset=["RH"])

        assert len(y_vals) == len(filtered)
        np.testing.assert_array_almost_equal(y_vals, filtered["RH"].values, decimal=2)

    def test_plotScatter_matches_data(self, lf_toolkit, lowfreq_df):
        df = lowfreq_df.copy()
        ax = lf_toolkit.presentation.dailyPlots.plotScatter(df, plotField="RH")

        y_vals = []
        for coll in ax.collections:
            offsets = coll.get_offsets()
            if len(offsets) > 0:
                y_vals.extend(offsets[:, 1])

        filtered = df["RH"].where(df["RH"] > -5000).dropna()
        matches = sum(np.isin(np.round(y_vals, 2), np.round(filtered.values, 2)))
        assert matches > 0


class TestPresentationEdgeCases:
    def test_scatter_empty_dataframe(self, lf_toolkit):
        df = pd.DataFrame(columns=["datetime", "RH"])
        df["datetime"] = pd.to_datetime(df["datetime"])
        ax = lf_toolkit.presentation.dailyPlots.plotScatter(df, plotField="RH")
        total = sum(len(coll.get_offsets()) for coll in ax.collections)
        assert total == 0

    def test_scatter_nan_and_outliers(self, lf_toolkit, lowfreq_df):
        df = lowfreq_df.copy()
        df.iloc[0, df.columns.get_loc("RH")] = np.nan
        df.iloc[1, df.columns.get_loc("RH")] = -9999

        ax = lf_toolkit.presentation.dailyPlots.plotScatter(df, plotField="RH")
        y_vals = []
        for coll in ax.collections:
            offsets = coll.get_offsets()
            if len(offsets) > 0:
                y_vals.extend(offsets[:, 1])

        assert len(y_vals) < len(df)

    def test_scatter_WS_field(self, lf_toolkit, lowfreq_df):
        df = lowfreq_df.copy()
        ax = lf_toolkit.presentation.dailyPlots.plotScatter(df, plotField="WS")
        y_vals = []
        for coll in ax.collections:
            offsets = coll.get_offsets()
            if len(offsets) > 0:
                y_vals.extend(offsets[:, 1])

        filtered = df["WS"].where(df["WS"] > -5000).dropna()
        matches = sum(np.isin(np.round(y_vals, 2), np.round(filtered.values, 2)))
        assert matches > 0


class TestPresentationSeasonalPlots:
    def test_plotProbContourf_bySeason(self, lf_toolkit, lowfreq_df):
        df = lowfreq_df.copy()
        ax = lf_toolkit.presentation.seasonalPlots.plotProbContourf_bySeason(
            df, plotField="RH"
        )
        assert ax is not None
        assert ax.shape == (2, 2)

    def test_contourf_distribution_ranges(self, lf_toolkit, lowfreq_df):
        df = lowfreq_df.copy()
        df["datetime"] = pd.to_datetime(df["datetime"], utc=True)
        df = df.set_index("datetime")

        CS, CFS, ax = lf_toolkit.presentation.dailyPlots.plotProbContourf(
            df, plotField="RH"
        )
        y_min, y_max = ax.get_ylim()
        filtered = df["RH"].where(df["RH"] > -5000).dropna()
        assert filtered.min() <= y_max
        assert filtered.max() >= y_min


class TestPresentationSavePlot:
    def test_scatter_creates_non_empty_image(self, lf_toolkit, lowfreq_df):
        df = lowfreq_df.copy()
        ax = lf_toolkit.presentation.dailyPlots.plotScatter(df, plotField="RH")
        with tempfile.NamedTemporaryFile(suffix=".png") as tmp:
            ax.figure.savefig(tmp.name)
            assert os.path.getsize(tmp.name) > 1000
