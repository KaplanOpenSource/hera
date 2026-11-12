"""measurements/meteorology/highfreqdata/analysis/analysislayer.py:
``RawdataAnalysis``, the factory layer the high-frequency toolkit exposes as
``toolkit.analysis``.

It has no numerics of its own -- it resolves a data-source name to a frame
and hands that frame, plus an identifier dict assembled from the measurement
geometry, to one of the three calculators.  So the assertions here are about
the contract: what ``_resolveData`` accepts and rejects, exactly which keys
the identifier carries, and that the calculators come back wired to the
resolved data.

Most tests drive the layer through a minimal stand-in data layer, because
``_resolveData``'s whole job is the two-way branch on its argument.  One
class goes through the real ``HighFreqToolKit`` from the ``unit_toolkit_factory``
fixture, so the by-name path is exercised against a genuinely registered
data source rather than a stub.

Two defects surfaced, both of the same shape -- a documented parameter that
the body never reads:

* B301: ``singlePointTurbulenceStatistics`` and ``AveragingCalculator``
  both accept ``inmemory`` and document it as "Whether to use pandas
  (in-memory) mode", but neither body mentions it again.  It reaches
  neither the identifier nor the calculator, so the two modes are
  indistinguishable and a caller asking for in-memory silently gets
  whatever the frame already was.
* B302: ``AveragingCalculator`` accepts and documents ``isMissingData``
  but -- unlike its sibling ``singlePointTurbulenceStatistics``, whose
  identifier does carry it -- drops it on the floor.  The resulting
  metadata cannot be told apart from a complete-data run, which also makes
  the two runs collide in the Cache lookup ``AbstractCalculator`` builds
  from that metadata.
"""
import numpy
import pandas
import pytest
import dask.dataframe

from hera.measurements.meteorology.highfreqdata.analysis import analysislayer
from hera.measurements.meteorology.highfreqdata.analysis.analysislayer import (
    RawdataAnalysis,
)
from hera.measurements.meteorology.highfreqdata.analysis.meandatacalculator import (
    AveragingCalculator,
    MeanDataCalculator,
)
from hera.measurements.meteorology.highfreqdata.analysis.turbulencestatistics import (
    singlePointTurbulenceStatistics,
)

PROJECT = "UNIT_TEST_PROJECT"
GEOMETRY = dict(samplingWindow="30s", start=pandas.Timestamp("2020-01-01"),
                end=pandas.Timestamp("2020-01-01 00:02:00"),
                height=10.0, buildingHeight=5.0, averagedHeight=3.0)

IDENTIFIER_KEYS = {"projectName", "samplingWindow", "height", "buildingHeight",
                   "averagedHeight", "start", "end", "filters",
                   "dataSource1", "dataSource2"}


class StubDataLayer:
    """The smallest thing ``RawdataAnalysis`` needs: a name and a loader."""

    def __init__(self, sources=None, projectName=PROJECT):
        self.projectName = projectName
        self._sources = dict(sources or {})
        self.requested = []

    def getDataSourceData(self, name):
        self.requested.append(name)
        return self._sources.get(name)


def _sonic(n=120, seed=0):
    rng = numpy.random.default_rng(seed)
    index = pandas.date_range("2020-01-01 00:00:00", periods=n, freq="1s", name="Time")
    return pandas.DataFrame({
        "u": 3.0 + rng.normal(0, 0.4, n),
        "v": 1.0 + rng.normal(0, 0.4, n),
        "w": rng.normal(0, 0.2, n),
        "T": 20.0 + rng.normal(0, 0.1, n),
    }, index=index)


@pytest.fixture()
def sonic():
    return _sonic()


@pytest.fixture()
def analyser(sonic):
    return RawdataAnalysis(StubDataLayer({"SONIC_1": sonic}))


@pytest.mark.unit
class TestConstruction:
    def test_the_datalayer_argument_is_exposed_as_the_datalayer(self):
        layer = StubDataLayer()
        assert RawdataAnalysis(layer).datalayer is layer

    def test_the_datalayer_starts_out_unset_on_the_class(self):
        assert RawdataAnalysis._datalayer is None

    def test_two_analysers_keep_their_own_datalayers(self):
        first, second = StubDataLayer(projectName="A"), StubDataLayer(projectName="B")
        assert RawdataAnalysis(first).datalayer.projectName == "A"
        assert RawdataAnalysis(second).datalayer.projectName == "B"

    def test_the_real_toolkit_wires_itself_in_as_the_datalayer(self, unit_toolkit_factory):
        from hera import toolkitHome

        toolkit = unit_toolkit_factory(toolkitHome.METEOROLOGY_HIGHFREQ)
        assert isinstance(toolkit.analysis, RawdataAnalysis)
        assert toolkit.analysis.datalayer is toolkit


@pytest.mark.unit
class TestResolveData:
    def test_a_pandas_frame_is_passed_straight_through(self, analyser, sonic):
        assert analyser._resolveData(sonic) is sonic

    def test_a_dask_frame_is_passed_straight_through(self, analyser, sonic):
        lazy = dask.dataframe.from_pandas(sonic, npartitions=1)
        assert analyser._resolveData(lazy) is lazy

    def test_a_name_is_loaded_from_the_data_layer(self, analyser, sonic):
        assert analyser._resolveData("SONIC_1") is sonic
        assert analyser.datalayer.requested == ["SONIC_1"]

    def test_an_unknown_name_names_itself_and_the_project_in_the_error(self, analyser):
        with pytest.raises(ValueError, match="MISSING_SONIC"):
            analyser._resolveData("MISSING_SONIC")
        with pytest.raises(ValueError, match=PROJECT):
            analyser._resolveData("MISSING_SONIC")

    @pytest.mark.parametrize("bad", [None, 7, 3.5, ["u", "v"], {"u": 1}])
    def test_anything_that_is_neither_a_name_nor_a_frame_is_rejected(self, analyser, bad):
        with pytest.raises(ValueError, match="data source name"):
            analyser._resolveData(bad)

    def test_the_rejection_message_reports_the_offending_type(self, analyser):
        with pytest.raises(ValueError, match="int"):
            analyser._resolveData(7)

    def test_a_series_is_not_accepted_as_a_frame(self, analyser, sonic):
        with pytest.raises(ValueError, match="data source name"):
            analyser._resolveData(sonic["u"])


@pytest.mark.unit
class TestResolveDataAgainstTheRealToolkit:
    @pytest.fixture()
    def wired(self, unit_toolkit_factory, tmp_path, sonic):
        from hera import toolkitHome

        toolkit = unit_toolkit_factory(toolkitHome.METEOROLOGY_HIGHFREQ)
        target = tmp_path / "sonic.parquet"
        sonic.to_parquet(target)
        toolkit.addDataSource("SONIC_1", str(target), "parquet", version=[0, 0, 1])
        toolkit.setDataSourceDefaultVersion("SONIC_1", [0, 0, 1])
        return toolkit

    def test_a_registered_data_source_is_resolved_by_name(self, wired, sonic):
        resolved = wired.analysis._resolveData("SONIC_1")
        assert list(resolved.columns) == list(sonic.columns)
        assert len(resolved) == len(sonic)

    def test_an_unregistered_name_is_reported_against_the_real_project(self, wired):
        with pytest.raises(ValueError, match="NOT_REGISTERED"):
            wired.analysis._resolveData("NOT_REGISTERED")


@pytest.mark.unit
class TestSinglePointTurbulenceStatisticsFactory:
    def _build(self, analyser, data, **kwargs):
        merged = dict(GEOMETRY)
        merged.update(kwargs)
        return analyser.singlePointTurbulenceStatistics(sonicData=data, **merged)

    def test_it_returns_a_turbulence_calculator(self, analyser, sonic):
        assert isinstance(self._build(analyser, sonic), singlePointTurbulenceStatistics)

    def test_the_calculator_holds_the_frame_it_was_given(self, analyser, sonic):
        assert self._build(analyser, sonic).RawData is sonic

    def test_a_data_source_name_is_resolved_before_construction(self, analyser, sonic):
        assert self._build(analyser, "SONIC_1").RawData is sonic

    def test_the_project_name_comes_from_the_data_layer(self, analyser, sonic):
        assert self._build(analyser, sonic).metaData["projectName"] == PROJECT

    def test_the_measurement_geometry_is_recorded_verbatim(self, analyser, sonic):
        meta = self._build(analyser, sonic).metaData
        for key, value in GEOMETRY.items():
            assert meta[key] == value

    def test_the_identifier_declares_the_missing_data_flag(self, analyser, sonic):
        assert self._build(analyser, sonic).metaData["isMissingData"] is False
        assert self._build(analyser, sonic, isMissingData=True).metaData["isMissingData"] is True

    def test_the_missing_data_flag_reaches_the_calculator_property(self, analyser, sonic):
        assert self._build(analyser, sonic, isMissingData=True).isMissingData is True

    def test_the_filter_and_data_source_slots_start_empty(self, analyser, sonic):
        meta = self._build(analyser, sonic).metaData
        assert meta["filters"] is None
        assert meta["dataSource1"] is None
        assert meta["dataSource2"] is None

    def test_the_identifier_carries_exactly_the_declared_keys(self, analyser, sonic):
        meta = self._build(analyser, sonic).metaData
        assert set(meta) == IDENTIFIER_KEYS | {"isMissingData"}

    def test_extra_keywords_are_added_to_the_identifier(self, analyser, sonic):
        meta = self._build(analyser, sonic, station="YAVNEEL", instrument="CSAT3").metaData
        assert meta["station"] == "YAVNEEL"
        assert meta["instrument"] == "CSAT3"

    def test_an_extra_keyword_can_override_a_declared_slot(self, analyser, sonic):
        meta = self._build(analyser, sonic, filters=["u > 0"]).metaData
        assert meta["filters"] == ["u > 0"]

    def test_two_calculators_do_not_share_an_identifier(self, analyser, sonic):
        first = self._build(analyser, sonic, height=10.0)
        second = self._build(analyser, sonic, height=20.0)
        assert first.metaData["height"] == 10.0
        assert second.metaData["height"] == 20.0
        assert first.metaData is not second.metaData

    def test_the_calculator_is_ready_to_compute_fluctuations(self, analyser, sonic):
        calculator = self._build(analyser, sonic)
        calculator.fluctuations()
        assert "up" in calculator.getData().columns

    @pytest.mark.xfail(
        strict=True,
        reason="B301: `inmemory` is a documented parameter of "
               "RawdataAnalysis.singlePointTurbulenceStatistics but the body "
               "never reads it -- it reaches neither the identifier nor the "
               "calculator, so asking for in-memory mode has no effect at "
               "all. See the consolidated findings issue.",
    )
    def test_asking_for_in_memory_mode_should_change_something(self, analyser, sonic):
        default = self._build(analyser, sonic, inmemory=False)
        requested = self._build(analyser, sonic, inmemory=True)
        assert requested.metaData != default.metaData

    def test_the_in_memory_flag_currently_leaves_no_trace(self, analyser, sonic):
        """Characterisation of B301."""
        default = self._build(analyser, sonic, inmemory=False)
        requested = self._build(analyser, sonic, inmemory=True)
        assert requested.metaData == default.metaData
        assert "inmemory" not in requested.metaData


@pytest.mark.unit
class TestAveragingCalculatorFactory:
    def _build(self, analyser, data, **kwargs):
        merged = dict(GEOMETRY)
        merged.update(kwargs)
        return analyser.AveragingCalculator(deviceNameOrData=data, **merged)

    def test_it_returns_an_averaging_calculator(self, analyser, sonic):
        assert isinstance(self._build(analyser, sonic), AveragingCalculator)

    def test_a_data_source_name_is_resolved_before_construction(self, analyser, sonic):
        assert self._build(analyser, "SONIC_1").RawData is sonic

    def test_the_project_name_and_geometry_are_recorded(self, analyser, sonic):
        meta = self._build(analyser, sonic).metaData
        assert meta["projectName"] == PROJECT
        for key, value in GEOMETRY.items():
            assert meta[key] == value

    def test_the_identifier_carries_exactly_the_declared_keys(self, analyser, sonic):
        assert set(self._build(analyser, sonic).metaData) == IDENTIFIER_KEYS

    def test_extra_keywords_are_added_to_the_identifier(self, analyser, sonic):
        meta = self._build(analyser, sonic, station="YAVNEEL").metaData
        assert meta["station"] == "YAVNEEL"

    def test_the_windowed_means_are_computed_on_construction(self, analyser, sonic):
        averaged = self._build(analyser, sonic)
        expected = sonic.resample("30s").mean()
        assert list(averaged.getData().columns) == [c + "_bar" for c in sonic.columns]
        assert numpy.allclose(averaged.getData()["u_bar"], expected["u"])

    def test_one_row_per_sampling_window(self, analyser, sonic):
        averaged = self._build(analyser, sonic, samplingWindow="60s")
        assert len(averaged.getData()) == 2

    @pytest.mark.xfail(
        strict=True,
        reason="B302: RawdataAnalysis.AveragingCalculator accepts and "
               "documents `isMissingData` but, unlike its sibling "
               "singlePointTurbulenceStatistics, never puts it in the "
               "identifier -- so a gappy run is indistinguishable from a "
               "complete one and the two collide in the Cache query "
               "AbstractCalculator builds from that metadata. See the "
               "consolidated findings issue.",
    )
    def test_the_missing_data_flag_should_reach_the_identifier(self, analyser, sonic):
        assert self._build(analyser, sonic, isMissingData=True).metaData["isMissingData"] is True

    def test_the_missing_data_flag_is_currently_dropped(self, analyser, sonic):
        """Characterisation of B302."""
        gappy = self._build(analyser, sonic, isMissingData=True)
        complete = self._build(analyser, sonic, isMissingData=False)
        assert "isMissingData" not in gappy.metaData
        assert gappy.metaData == complete.metaData

    def test_the_in_memory_flag_currently_leaves_no_trace(self, analyser, sonic):
        """Characterisation of B301 on the averaging factory."""
        assert self._build(analyser, sonic, inmemory=True).metaData == \
               self._build(analyser, sonic, inmemory=False).metaData


@pytest.mark.unit
class TestMeanDataCalculatorFactory:
    def _turbulence(self, analyser, sonic):
        merged = dict(GEOMETRY)
        return analyser.singlePointTurbulenceStatistics(sonicData=sonic, **merged)

    def test_it_returns_a_mean_data_calculator_built_from_the_second_moments(self,
                                                                            analyser,
                                                                            sonic):
        calculator = analyser.MeanDataCalculator(self._turbulence(analyser, sonic))
        assert isinstance(calculator, MeanDataCalculator)
        for moment in ("uu", "vv", "ww", "uw", "wT"):
            assert moment in calculator.MeanData.columns

    def test_the_mean_data_is_clipped_to_the_metadata_interval(self, analyser, sonic):
        calculator = analyser.MeanDataCalculator(self._turbulence(analyser, sonic))
        assert (calculator.MeanData.index >= GEOMETRY["start"]).all()
        assert (calculator.MeanData.index < GEOMETRY["end"]).all()

    def test_the_turbulence_metadata_is_inherited(self, analyser, sonic):
        calculator = analyser.MeanDataCalculator(self._turbulence(analyser, sonic))
        assert calculator.metaData["height"] == GEOMETRY["height"]
        assert calculator.metaData["projectName"] == PROJECT

    def test_extra_metadata_overrides_what_the_calculator_carried(self, analyser, sonic):
        calculator = analyser.MeanDataCalculator(self._turbulence(analyser, sonic),
                                                 height=99.0)
        assert calculator.metaData["height"] == 99.0

    def test_every_argument_is_forwarded_in_the_documented_order(self, analyser,
                                                                 monkeypatch):
        seen = {}

        def recorder(turb, mode_turb, average, mode_average, **metadata):
            seen.update(turb=turb, mode_turb=mode_turb, average=average,
                        mode_average=mode_average, metadata=metadata)
            return "built"

        monkeypatch.setattr(analysislayer, "MeanDataCalculator", recorder)
        result = analyser.MeanDataCalculator("TURB", "from_db_and_save",
                                             "AVG", "not_from_db_and_save",
                                             station="YAVNEEL")
        assert result == "built"
        assert seen == dict(turb="TURB", mode_turb="from_db_and_save",
                            average="AVG", mode_average="not_from_db_and_save",
                            metadata=dict(station="YAVNEEL"))

    def test_the_default_modes_are_local_compute_and_an_unset_average_mode(self,
                                                                          analyser,
                                                                          monkeypatch):
        seen = {}

        def recorder(turb, mode_turb, average, mode_average, **metadata):
            seen.update(turb=turb, mode_turb=mode_turb, average=average,
                        mode_average=mode_average)
            return None

        monkeypatch.setattr(analysislayer, "MeanDataCalculator", recorder)
        analyser.MeanDataCalculator()
        assert seen == dict(turb=None, mode_turb="not_from_db_and_not_save",
                            average=None, mode_average=None)

    def test_an_unusable_turbulence_argument_is_reported(self, analyser):
        with pytest.raises(ValueError, match="singlePointTurbulenceStatistics"):
            analyser.MeanDataCalculator(TurbCalcOrData="not a calculator")
