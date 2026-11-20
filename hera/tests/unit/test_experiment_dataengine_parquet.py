"""measurements/experiment/dataEngine.py: ``parquetDataEngineHera``.

``test_experiment_analysis_dataengine.py`` covers the factory and the two
MongoDB engines (B76/B77); this file covers the third, hera-native one --
``__init__``, ``getData`` and ``getDataFromTrial``.

Construction strategy
---------------------
The real constructor is used throughout: ``parquetDataEngineHera`` is a
``datalayer.Project`` subclass, so it works unchanged against the unit
layer's mongomock database, and the device data it reads is registered as
ordinary ``type='Experiment_rawData'`` measurement documents holding
pickled frames (pandas for the common path, a one-partition Dask frame
for the ``autoCompute`` path). Only ``getDataFromTrial`` is driven with a
stand-in experiment object, because it never gets far enough to touch a
real one (B165) and the attributes it *would* read are all it needs.

``unit_project`` is requested by every fixture before the engine is
built: ``parquetDataEngineHera.__init__`` calls
``Project.__init__(projectName=...)`` without forwarding a
``filesDirectory``, so on a project with no persisted config it would
default to -- and create -- ``~/.hera/<projectName>``.

Deliberately not covered
------------------------
``pandasDataEngineDB`` / ``daskDataEngineDB``: they open raw ``pymongo``
connections to a real server, and their defects are already pinned as
B76/B77 in ``test_experiment_analysis_dataengine.py``.

Bugs pinned here
----------------
* B165: ``parquetDataEngineHera.getDataFromTrial`` is dead code. It
  reads ``self.experimentObj.experimentSetup.trialSet[trialSet][trialName]``,
  but the object it is given -- ``experimentSetupWithData`` -- has no
  ``experimentSetup`` attribute at all (it *is* the experiment setup, and
  exposes the trial sets directly as ``.trialSet``), so every call raises
  ``AttributeError`` on the very first statement that uses the argument.
  Two further defects sit behind it on the same code path: the
  ``trialSet=None`` default assigns ``self.experimentObj.trialSet``, i.e.
  the whole *dict* of trial sets, where a trial-set *name* is required
  (an unhashable key), and the ``withMetadata`` branch calls
  ``entitiesTable()`` although ``entitiesTable`` is a property, which
  would raise ``TypeError``. The same three defects are copy-pasted into
  ``pandasDataEngineDB.getDataFromTrial`` and
  ``daskDataEngineDB.getDataFromTrial``.
"""
import pandas
import pytest

from hera.measurements.experiment.dataEngine import parquetDataEngineHera
from hera.tests.unit.conftest import UNIT_PROJECT_NAME

TZ = "Asia/Jerusalem"
CONFIGURATION = {"experimentName": "UNITEXP"}


def _frame():
    """Twenty minutes of alternating sonic_1/sonic_2 rows, tz-aware."""
    index = pandas.date_range("2020-01-01 00:00:00", periods=20, freq="1min", tz=TZ)
    return pandas.DataFrame(
        {
            "deviceName": ["sonic_1" if i % 2 == 0 else "sonic_2" for i in range(20)],
            "TC": [20.0 + i for i in range(20)],
        },
        index=index,
    ).rename_axis("timestamp")


def _register(project, tmp_path, frame, name, deviceType="sonic", **desc):
    """Register a frame (pandas or dask) as an Experiment_rawData document."""
    import cloudpickle

    path = tmp_path / f"{name}.pkl"
    with open(str(path), "wb") as handle:
        cloudpickle.dump(frame, handle)
    project.addMeasurementsDocument(
        resource=str(path),
        dataFormat="pickle",
        type="Experiment_rawData",
        desc=dict(experimentName="UNITEXP", deviceType=deviceType, **desc),
    )


@pytest.fixture()
def engine(unit_project):
    """A real engine over the in-memory database, with no experiment object."""
    return parquetDataEngineHera(UNIT_PROJECT_NAME, CONFIGURATION, experimentObj=None)


@pytest.fixture()
def sharedFileData(unit_project, tmp_path):
    """One document holding every device of the type (``perDevice=False``)."""
    _register(unit_project, tmp_path, _frame(), "shared")


@pytest.mark.unit
class TestConstruction:
    def test_the_experiment_name_comes_from_the_data_source_configuration(self, engine):
        assert engine.experimentName == "UNITEXP"

    def test_the_project_name_is_forwarded_to_the_datalayer(self, engine):
        assert engine.projectName == UNIT_PROJECT_NAME

    def test_the_experiment_object_is_stored_as_given(self, unit_project):
        marker = object()
        engine = parquetDataEngineHera(UNIT_PROJECT_NAME, CONFIGURATION, experimentObj=marker)
        assert engine.experimentObj is marker

    def test_a_configuration_without_an_experiment_name_raises(self, unit_project):
        with pytest.raises(KeyError, match="experimentName"):
            parquetDataEngineHera(UNIT_PROJECT_NAME, {}, experimentObj=None)


@pytest.mark.unit
class TestGetData:
    def test_an_unknown_device_type_returns_an_empty_frame(self, engine, sharedFileData):
        data = engine.getData(deviceType="noSuchDevice")
        assert isinstance(data, pandas.DataFrame)
        assert data.empty

    def test_a_project_without_documents_returns_an_empty_frame(self, engine):
        assert engine.getData(deviceType="sonic").empty

    def test_it_returns_the_whole_document_when_nothing_is_filtered(self, engine, sharedFileData):
        assert len(engine.getData(deviceType="sonic")) == 20

    def test_a_device_name_filters_the_rows_of_the_shared_document(self, engine, sharedFileData):
        data = engine.getData(deviceType="sonic", deviceName="sonic_1")
        assert set(data.deviceName) == {"sonic_1"}
        assert len(data) == 10

    def test_a_start_time_alone_trims_the_beginning(self, engine, sharedFileData):
        data = engine.getData(
            deviceType="sonic", startTime=pandas.Timestamp("2020-01-01 00:15:00", tz=TZ))
        assert len(data) == 5

    def test_an_end_time_alone_trims_the_end(self, engine, sharedFileData):
        data = engine.getData(
            deviceType="sonic", endTime=pandas.Timestamp("2020-01-01 00:04:00", tz=TZ))
        assert len(data) == 5

    def test_both_ends_select_an_inclusive_window(self, engine, sharedFileData):
        data = engine.getData(
            deviceType="sonic",
            startTime=pandas.Timestamp("2020-01-01 00:05:00", tz=TZ),
            endTime=pandas.Timestamp("2020-01-01 00:07:00", tz=TZ),
        )
        assert len(data) == 3

    def test_the_device_filter_and_the_time_window_compose(self, engine, sharedFileData):
        data = engine.getData(
            deviceType="sonic",
            deviceName="sonic_2",
            startTime=pandas.Timestamp("2020-01-01 00:00:00", tz=TZ),
            endTime=pandas.Timestamp("2020-01-01 00:05:00", tz=TZ),
        )
        assert set(data.deviceName) == {"sonic_2"}
        assert len(data) == 3

    def test_extra_query_keywords_select_among_documents(self, unit_project, engine, tmp_path):
        _register(unit_project, tmp_path, _frame(), "calibrated", calibrated=True)
        _register(unit_project, tmp_path, _frame().head(4), "raw", calibrated=False)

        assert len(engine.getData(deviceType="sonic", calibrated=False)) == 4
        assert len(engine.getData(deviceType="sonic", calibrated=True)) == 20

    def test_it_reads_the_first_matching_document(self, unit_project, engine, tmp_path):
        _register(unit_project, tmp_path, _frame(), "first")
        _register(unit_project, tmp_path, _frame().head(3), "second")
        assert len(engine.getData(deviceType="sonic")) == 20


@pytest.mark.unit
class TestGetDataPerDevice:
    def test_the_per_device_branch_requires_a_device_name(self, engine, sharedFileData):
        with pytest.raises(AssertionError, match="deviceName should be defined"):
            engine.getData(deviceType="sonic", perDevice=True)

    def test_the_per_device_branch_selects_the_documents_of_that_device(
            self, unit_project, engine, tmp_path):
        _register(unit_project, tmp_path, _frame().head(6), "perDevice_1", deviceName="sonic_1")
        _register(unit_project, tmp_path, _frame().head(2), "perDevice_2", deviceName="sonic_2")

        data = engine.getData(deviceType="sonic", deviceName="sonic_1", perDevice=True)
        assert len(data) == 6

    def test_the_per_device_branch_does_not_filter_the_rows_again(
            self, unit_project, engine, tmp_path):
        """The whole per-device document is returned, mixed device names included."""
        _register(unit_project, tmp_path, _frame(), "perDevice_all", deviceName="sonic_1")
        data = engine.getData(deviceType="sonic", deviceName="sonic_1", perDevice=True)
        assert set(data.deviceName) == {"sonic_1", "sonic_2"}


@pytest.mark.unit
class TestGetDataAutoCompute:
    """autoCompute calls ``.compute()``, so it needs a Dask-backed document."""

    @pytest.fixture()
    def daskData(self, unit_project, tmp_path):
        import dask.dataframe

        _register(unit_project, tmp_path,
                  dask.dataframe.from_pandas(_frame(), npartitions=1), "dask")

    def test_without_autocompute_the_dask_frame_is_returned_lazily(self, engine, daskData):
        import dask.dataframe

        data = engine.getData(deviceType="sonic")
        assert isinstance(data, dask.dataframe.DataFrame)

    def test_autocompute_returns_a_pandas_frame(self, engine, daskData):
        data = engine.getData(deviceType="sonic", autoCompute=True)
        assert isinstance(data, pandas.DataFrame)
        assert len(data) == 20

    def test_autocompute_applies_the_time_window_first(self, engine, daskData):
        data = engine.getData(
            deviceType="sonic",
            startTime=pandas.Timestamp("2020-01-01 00:05:00", tz=TZ),
            endTime=pandas.Timestamp("2020-01-01 00:07:00", tz=TZ),
            autoCompute=True,
        )
        assert len(data) == 3


class _FakeTrial:
    def __init__(self):
        self.properties = {
            "TrialStart": pandas.Timestamp("2020-01-01 00:00:00", tz=TZ),
            "TrialEnd": pandas.Timestamp("2020-01-01 00:09:00", tz=TZ),
        }


class _FakeExperiment:
    """Exposes exactly what ``experimentSetupWithData`` exposes: the trial
    sets hang directly off the object, and there is no ``experimentSetup``."""

    def __init__(self):
        self.trialSet = {"Measurements": {"T1": _FakeTrial()}}


@pytest.mark.unit
class TestGetDataFromTrial:
    """B165: the method reaches for an attribute that does not exist."""

    @pytest.fixture()
    def trialEngine(self, unit_project):
        return parquetDataEngineHera(
            UNIT_PROJECT_NAME, CONFIGURATION, experimentObj=_FakeExperiment())

    def test_the_experiment_class_has_no_experiment_setup_attribute(self):
        """Characterisation of B165: the attribute the method reads."""
        from hera.measurements.experiment.experiment import experimentSetupWithData

        assert not hasattr(experimentSetupWithData, "experimentSetup")
        assert hasattr(experimentSetupWithData, "trialSet")

    @pytest.mark.xfail(
        strict=True,
        reason="B165: getDataFromTrial reads "
               "self.experimentObj.experimentSetup.trialSet[...], but the "
               "experiment object it is handed (experimentSetupWithData) has "
               "no experimentSetup attribute -- it exposes .trialSet "
               "directly -- so every call raises AttributeError. "
               "See the consolidated findings issue.",
    )
    def test_it_should_read_the_trials_data(self, trialEngine, sharedFileData):
        data = trialEngine.getDataFromTrial(
            deviceType="sonic", trialName="T1", trialSet="Measurements",
            withMetadata=False)
        assert len(data) == 10

    def test_it_currently_raises_an_attribute_error(self, trialEngine, sharedFileData):
        """Characterisation of B165."""
        with pytest.raises(AttributeError, match="experimentSetup"):
            trialEngine.getDataFromTrial(
                deviceType="sonic", trialName="T1", trialSet="Measurements",
                withMetadata=False)

    def test_the_default_trial_set_is_the_dict_of_trial_sets_not_a_name(self, trialEngine, sharedFileData):
        """Characterisation of B165: omitting trialSet is dead too.

        ``trialSet = self.experimentObj.trialSet if trialSet is None`` puts
        the whole mapping where a key is expected; the AttributeError on the
        next line is simply reached first.
        """
        with pytest.raises(AttributeError, match="experimentSetup"):
            trialEngine.getDataFromTrial(deviceType="sonic", trialName="T1")
