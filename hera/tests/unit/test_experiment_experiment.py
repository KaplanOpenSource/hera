"""measurements/experiment/experiment.py: the five data-aware argos
subclasses -- ``experimentSetupWithData``, ``TrialSetWithData``,
``TrialWithdata``, ``EntityTypeWithData`` and ``EntityWithData``.

Construction strategy
---------------------
No ``__new__`` bypass and no argos stubbing was needed. The real argos
package is importable, and ``argosDataObjects.ExperimentZipFile`` only
wants a zip archive containing a single ``data.json`` member, so every
test below drives the *real* constructor chain end to end:

* ``_experiment_setup()`` builds a minimal but genuine **v3.0.0** argos
  setup dict (``version``/``name``/``trialTypes``/``deviceTypes``/
  ``imageStandalone``), shaped after the real exports under
  ``argos/experimentSetup/example_exp`` -- one trial set with parsed
  ``TrialStart``/``TrialEnd`` ``datetime-local`` properties, one device
  type with a ``Constant``-scope ``StoreDataPerDevice`` attribute and two
  devices, and two devices-on-trial carrying a ``location`` plus a
  ``Number`` ``height`` attribute. ``_fix_json_version_3_0_0`` then does
  the real migration.
* ``_write_experiment()`` lays that out on disk exactly as
  ``experimentSetupWithData.__init__`` expects:
  ``<root>/runtimeExperimentData/Datasources_Configurations.json`` plus
  ``<root>/runtimeExperimentData/<experimentName>.zip``, both under
  ``tmp_path``.
* the data engine is the real ``parquetDataEngineHera`` against
  mongomock; the device data it serves is registered as ordinary
  measurement documents (``type='Experiment_rawData'``) holding pickled
  frames, so ``getData`` returns real pandas.

One ordering subtlety is load-bearing and is why the fixtures depend on
``unit_project``: ``parquetDataEngineHera.__init__`` calls
``datalayer.Project(projectName=...)`` *without* forwarding a
``filesDirectory``, so on a project that has no persisted
``filesDirectory`` config it falls back to ``~/.hera/<projectName>`` and
creates it -- which the unit layer's ``_no_home_writes`` guard rightly
fails. Building ``unit_project`` first persists the tmp files directory
in the (mongomock) project config, which the engine's ``Project`` then
picks up.

Covered
-------
``experimentSetupWithData``: ``__init__`` (both validation failures, the
configuration override, the cache directory, the toolkit name),
``configuration``, ``name``, ``analysis``, ``presentation``,
``getExperimentData``, ``_initTrialSets``, ``_initEntitiesTypes``,
``defaultTrialSet``, ``trialsOfDefaultTrialSet``,
``_initAnalysisAndPresentation``, ``getDataFromDateRange``,
``_process_row``, ``get_devices_image_coordinates``.
``TrialSetWithData``: ``__init__``, ``_initTrials``.
``TrialWithdata``: ``__init__``, ``getData`` (default and explicit time
range, device filter, the empty-window ``ValueError``, ``withMetadata``).
``EntityTypeWithData``: ``__init__``, ``_initEntities``, ``getData``,
``getDataTrial`` (both the shared-file and the per-device branch).
``EntityWithData``: ``__init__``, ``getData`` (both branches).

Deliberately not covered
------------------------
The ``PANDASDB``/``DASKDB`` engine variants of the same experiment: both
open raw ``pymongo`` connections to a real server in (or immediately
after) construction and their defects are already pinned as B76/B77 in
``test_experiment_analysis_dataengine.py``.
``experimentSetupWithData.getImage`` and the image-map plotting live in
argos / ``presentation.py``, not here.

Bugs pinned here
----------------
* B160: ``_process_row`` does ``pp = convertCRS(...)`` and then reads
  ``pp.x[0]`` / ``pp.y[0]``, but ``hera.measurements.GIS.utils.convertCRS``
  returns ``list(gdf.to_crs(...).geometry)`` -- a plain Python list of
  shapely Points, which has no ``.x``. Every call raises
  ``AttributeError: 'list' object has no attribute 'x'``, so the entire
  ``outputCRS=ITM`` branch of ``get_devices_image_coordinates`` is dead.
  The working access is ``pp[0].x`` / ``pp[0].y``.
* B161: ``get_devices_image_coordinates`` (and ``_process_row``) read
  the columns ``Latitude``/``Longitude``, but argos spreads a
  device-on-trial ``location`` into **lowercase** ``latitude``/
  ``longitude`` (``fillContained.spread_attributes``), so on any
  experiment whose device positions come from the map -- the normal case,
  and what every example export in argos does -- the lookup raises
  ``KeyError: 'Latitude'``. The method only ever works if the devices
  happen to carry free-text *attributes* literally named ``Latitude`` and
  ``Longitude``.
* B162: ``TrialWithdata.getData(withMetadata=True)`` merges the data
  against ``self.entitiesTable`` with ``right_on="entityName"``, but
  ``argos.Trial.entitiesTable`` names that column ``deviceItemName``
  (``entityName`` exists only on ``Experiment.entitiesTable`` and
  ``EntityType.entitiesTable``). Every metadata-joined trial read raises
  ``KeyError: 'entityName'``.
* B163: the two overridden initialisers index a metadata key that
  argos's own versions read defensively:
  ``TrialSetWithData._initTrials`` uses ``self._metadata["trials"]``
  (argos: ``.get('trials', [])``) and
  ``EntityTypeWithData._initEntities`` uses ``self._metadata["entities"]``
  (argos: ``.get('entities', [])``). A v3.0.0 export omits ``trials``
  entirely for a trial type that has none -- argos's own
  ``example_exp/exp_simple`` export does exactly that -- so loading such
  an experiment through hera raises ``KeyError: 'trials'`` from inside
  the constructor.
* B164: ``EntityTypeWithData.getDataTrial`` forwards
  ``perDevice=StoreDataPerDevice`` but never forwards a ``deviceName``,
  while ``parquetDataEngineHera.getData`` starts its per-device branch
  with ``assert deviceName, "If perDeivce=True then deviceName should be
  defined!"``. For every entity type whose data *is* stored per device --
  the only case in which the flag is set -- the call therefore raises
  ``AssertionError``. ``EntityWithData.getData`` gets this right, which
  is what makes the omission visible.
"""
import json
import os
import zipfile

import pandas
import pytest

from hera.measurements.experiment.experiment import (
    EntityTypeWithData,
    EntityWithData,
    TrialSetWithData,
    TrialWithdata,
    experimentSetupWithData,
)
from hera.measurements.GIS.utils import ITM, WSG84
from hera.tests.unit.conftest import UNIT_PROJECT_NAME

TRIAL_START = "2020-01-01 00:00:00"
TRIAL_END = "2020-01-01 00:09:00"
TZ = "Asia/Jerusalem"


# ---------------------------------------------------------------------------
# Synthetic (but genuine) argos v3.0.0 experiment fixtures
# ---------------------------------------------------------------------------

def _device_on_trial(name, longitude, latitude, coordinateAttributes=False):
    """One entry of a v3.0.0 ``devicesOnTrial`` list."""
    attributes = [{"name": "height", "value": "3"}]
    if coordinateAttributes:
        attributes += [
            {"name": "Latitude", "value": str(latitude)},
            {"name": "Longitude", "value": str(longitude)},
        ]
    return {
        "deviceTypeName": "sonic",
        "deviceItemName": name,
        "location": {"name": "map1", "coordinates": [longitude, latitude]},
        "attributes": attributes,
    }


def _experiment_setup(coordinateAttributes=False, withTrials=True, storePerDevice=False):
    """A minimal argos v3.0.0 setup dict.

    coordinateAttributes
        also give every device free-text ``Latitude``/``Longitude``
        attributes, which is the only shape ``get_devices_image_coordinates``
        can actually read (see B161).
    withTrials
        when False the trial type carries no ``trials`` key at all, which
        is what a real export looks like for an empty trial type.
    storePerDevice
        the value of the ``Constant``-scope ``StoreDataPerDevice``
        attribute.
    """
    attributeTypes = [
        {"name": "height", "type": "Number", "label": "height",
         "description": "", "scope": "Device"},
        {"name": "StoreDataPerDevice", "type": "String",
         "label": "StoreDataPerDevice", "description": "",
         "scope": "Constant", "defaultValue": storePerDevice},
    ]
    if coordinateAttributes:
        attributeTypes += [
            {"name": "Latitude", "type": "Number", "label": "Latitude",
             "description": "", "scope": "Device"},
            {"name": "Longitude", "type": "Number", "label": "Longitude",
             "description": "", "scope": "Device"},
        ]

    trialType = {
        "name": "Measurements",
        "description": "the measurement trials",
        "attributeTypes": [
            {"key": "TrialStart", "name": "TrialStart", "type": "datetime-local",
             "label": "TrialStart", "description": ""},
            {"key": "TrialEnd", "name": "TrialEnd", "type": "datetime-local",
             "label": "TrialEnd", "description": ""},
        ],
    }
    if withTrials:
        trialType["trials"] = [{
            "name": "T1",
            "createdDate": "2020-01-01T00:00:00.000Z",
            "properties": [
                {"key": "TrialStart", "val": TRIAL_START},
                {"key": "TrialEnd", "val": TRIAL_END},
            ],
            "devicesOnTrial": [
                _device_on_trial("sonic_1", 34.0, 32.0, coordinateAttributes),
                _device_on_trial("sonic_2", 34.5, 32.5, coordinateAttributes),
            ],
        }]

    return {
        "version": "3.0.0",
        "name": "UNITEXP",
        "description": "a synthetic experiment",
        "startDate": "2020-01-01T00:00:00.000Z",
        "endDate": "2020-01-02T00:00:00.000Z",
        "trialTypes": [trialType],
        "deviceTypes": [{
            "name": "sonic",
            "attributeTypes": attributeTypes,
            "devices": [
                {"name": "sonic_1", "attributes": [{"name": "height", "value": "3"}]},
                {"name": "sonic_2", "attributes": [{"name": "height", "value": "6"}]},
            ],
        }],
        "imageStandalone": [{
            "name": "map1", "filename": "map1.png",
            "left": 34.0, "right": 35.0, "lower": 32.0, "upper": 33.0,
            "width": 100, "height": 100,
        }],
    }


def _write_experiment(root, setup, experimentName="UNITEXP", configuration=None,
                      writeSetupFile=True, writeConfigurationFile=True):
    """Lay a setup dict out on disk the way ``experimentSetupWithData`` expects."""
    runtime = os.path.join(root, "runtimeExperimentData")
    os.makedirs(runtime, exist_ok=True)

    if writeConfigurationFile:
        config = {"experimentName": experimentName}
        config.update(configuration or {})
        with open(os.path.join(runtime, "Datasources_Configurations.json"), "w") as handle:
            json.dump(config, handle)

    if writeSetupFile:
        with zipfile.ZipFile(os.path.join(runtime, f"{experimentName}.zip"), "w") as archive:
            archive.writestr("data.json", json.dumps(setup))

    return root


@pytest.fixture()
def experimentFactory(tmp_path, unit_files_directory, unit_project):
    """Build a real ``experimentSetupWithData`` from a setup dict.

    ``unit_project`` is requested first on purpose: it persists the tmp
    files directory in the project config, so the ``Project`` that
    ``parquetDataEngineHera`` builds for itself (which forwards no
    ``filesDirectory``) does not fall back to ``~/.hera``.
    """
    counter = {"n": 0}

    def _build(setup=None, name="UNITEXP", defaultTrialSetName="Measurements", **kwargs):
        counter["n"] += 1
        root = str(tmp_path / f"exp{counter['n']}")
        setup = _experiment_setup() if setup is None else setup
        _write_experiment(root, setup, experimentName=name, **kwargs)
        return experimentSetupWithData(
            projectName=UNIT_PROJECT_NAME,
            pathToExperiment=root,
            filesDirectory=unit_files_directory,
            defaultTrialSetName=defaultTrialSetName,
        )

    return _build


@pytest.fixture()
def experiment(experimentFactory):
    return experimentFactory()


@pytest.fixture()
def coordinateExperiment(experimentFactory):
    """An experiment whose devices carry capitalised coordinate attributes."""
    return experimentFactory(_experiment_setup(coordinateAttributes=True))


def _deviceFrame():
    """Twenty minutes of alternating sonic_1/sonic_2 rows, tz-aware."""
    index = pandas.date_range("2020-01-01 00:00:00", periods=20, freq="1min", tz=TZ)
    return pandas.DataFrame(
        {
            "deviceName": ["sonic_1" if i % 2 == 0 else "sonic_2" for i in range(20)],
            "TC": [20.0 + i for i in range(20)],
        },
        index=index,
    ).rename_axis("timestamp")


def _registerDeviceData(project, tmp_path, frame, deviceType="sonic", **desc):
    """Register a pickled frame as an Experiment_rawData measurement document."""
    path = tmp_path / f"rawData_{deviceType}_{len(desc)}_{id(frame)}.pkl"
    frame.to_pickle(str(path))
    project.addMeasurementsDocument(
        resource=str(path),
        dataFormat="pickle",
        type="Experiment_rawData",
        desc=dict(experimentName="UNITEXP", deviceType=deviceType, **desc),
    )
    return frame


@pytest.fixture()
def deviceData(unit_project, tmp_path):
    """The shared-file (``perDevice=False``) raw-data document."""
    return _registerDeviceData(unit_project, tmp_path, _deviceFrame())


# ---------------------------------------------------------------------------
# experimentSetupWithData: construction
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestExperimentConstruction:
    def test_the_configuration_is_loaded_from_the_datasources_configuration_file(self, experiment):
        assert experiment.configuration["experimentName"] == "UNITEXP"

    def test_the_name_comes_from_the_configuration_and_not_from_the_setup(self, experimentFactory):
        setup = _experiment_setup()
        setup["name"] = "a different name inside the zip"
        exp = experimentFactory(setup)
        assert exp.name == "UNITEXP"

    def test_the_setup_is_the_migrated_v3_dictionary(self, experiment):
        assert sorted(experiment.setup) == ["entityTypes", "experiment", "maps", "shapes", "trialSets"]

    def test_a_missing_configuration_file_raises(self, experimentFactory):
        with pytest.raises(ValueError, match="configuration file does not exist"):
            experimentFactory(writeConfigurationFile=False)

    def test_a_missing_setup_zip_raises(self, experimentFactory):
        with pytest.raises(ValueError, match="setup file does not exist"):
            experimentFactory(writeSetupFile=False)

    def test_the_data_source_configuration_argument_is_merged_into_the_configuration(
            self, tmp_path, unit_files_directory, unit_project):
        root = _write_experiment(str(tmp_path / "override"), _experiment_setup())
        exp = experimentSetupWithData(
            projectName=UNIT_PROJECT_NAME,
            pathToExperiment=root,
            dataSourceConfiguration=dict(extraKey="extraValue"),
            filesDirectory=unit_files_directory,
        )
        assert exp.configuration["extraKey"] == "extraValue"
        assert exp.configuration["experimentName"] == "UNITEXP"

    def test_a_none_data_source_configuration_is_treated_as_an_empty_dict(
            self, tmp_path, unit_files_directory, unit_project):
        root = _write_experiment(str(tmp_path / "noneconf"), _experiment_setup())
        exp = experimentSetupWithData(
            projectName=UNIT_PROJECT_NAME,
            pathToExperiment=root,
            dataSourceConfiguration=None,
            filesDirectory=unit_files_directory,
        )
        assert exp.configuration["experimentName"] == "UNITEXP"

    def test_the_toolkit_name_is_derived_from_the_experiment_name(self, experiment):
        assert experiment.toolkitName == "UNITEXPToolkit"

    def test_the_cache_directory_is_created_under_the_given_files_directory(
            self, experiment, unit_files_directory):
        assert experiment.filesDirectory == os.path.join(unit_files_directory, "experimentCache")
        assert os.path.isdir(experiment.filesDirectory)

    def test_the_image_map_is_built_from_the_standalone_images(self, experiment):
        assert list(experiment.imageMap) == ["map1"]
        assert experiment.imageMap["map1"]["filename"] == "map1.png"


@pytest.mark.unit
class TestAnalysisAndPresentation:
    def test_both_layers_are_wired_up_at_construction(self, experiment):
        from hera.measurements.experiment.analysis import experimentAnalysis
        from hera.measurements.experiment.presentation import experimentPresentation

        assert isinstance(experiment.analysis, experimentAnalysis)
        assert isinstance(experiment.presentation, experimentPresentation)

    def test_the_presentation_layer_is_given_the_analysis_layer(self, experiment):
        assert experiment.presentation.analysis is experiment.analysis

    def test_init_analysis_and_presentation_replaces_both_layers(self, experiment):
        class _Analysis:
            def __init__(self, datalayer):
                self.datalayer = datalayer

        class _Presentation:
            def __init__(self, datalayer, analysis):
                self.datalayer = datalayer
                self.analysis = analysis

        experiment._initAnalysisAndPresentation(_Analysis, _Presentation)

        assert isinstance(experiment.analysis, _Analysis)
        assert isinstance(experiment.presentation, _Presentation)
        assert experiment.analysis.datalayer is experiment
        assert experiment.presentation.analysis is experiment.analysis


@pytest.mark.unit
class TestGetExperimentData:
    def test_the_default_engine_is_the_hera_parquet_engine(self, experiment):
        from hera.measurements.experiment.dataEngine import parquetDataEngineHera

        assert isinstance(experiment.getExperimentData(), parquetDataEngineHera)

    def test_the_engine_knows_the_experiment_name_and_object(self, experiment):
        engine = experiment.getExperimentData()
        assert engine.experimentName == "UNITEXP"
        assert engine.experimentObj is experiment


# ---------------------------------------------------------------------------
# _initTrialSets / TrialSetWithData / TrialWithdata construction
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestTrialSetInitialisation:
    def test_every_trial_set_in_the_setup_becomes_a_trial_set_with_data(self, experiment):
        assert list(experiment.trialSet) == ["Measurements"]
        assert isinstance(experiment.trialSet["Measurements"], TrialSetWithData)

    def test_the_trial_set_keeps_a_back_reference_to_the_experiment(self, experiment):
        assert experiment.trialSet["Measurements"].experiment is experiment

    def test_the_trial_set_shares_the_experiments_data_engine(self, experiment):
        trialSet = experiment.trialSet["Measurements"]
        assert trialSet._experimentData is experiment.getExperimentData()

    def test_the_trials_are_trial_with_data_objects(self, experiment):
        trial = experiment.trialSet["Measurements"]["T1"]
        assert isinstance(trial, TrialWithdata)
        assert trial.name == "T1"

    def test_every_trial_shares_the_same_data_engine(self, experiment):
        trial = experiment.trialSet["Measurements"]["T1"]
        assert trial._experimentData is experiment.getExperimentData()

    def test_the_trial_time_properties_are_parsed_from_the_setup(self, experiment):
        properties = experiment.trialSet["Measurements"]["T1"].properties
        assert properties["TrialStart"] == pandas.Timestamp(TRIAL_START, tz=TZ)
        assert properties["TrialEnd"] == pandas.Timestamp(TRIAL_END, tz=TZ)


@pytest.mark.unit
class TestTrialSetInitialisationIsNotDefensive:
    """B163: ``[...]`` where argos itself uses ``.get(..., [])``."""

    @pytest.mark.xfail(
        strict=True,
        reason="B163: TrialSetWithData._initTrials reads "
               "self._metadata['trials'] while argos's TrialSet._initTrials "
               "reads self._metadata.get('trials', []). A real v3.0.0 export "
               "omits the key for a trial type with no trials (see argos's "
               "own example_exp/exp_simple), so hera cannot load it. "
               "See the consolidated findings issue.",
    )
    def test_a_trial_set_without_trials_should_load_as_an_empty_trial_set(self, experimentFactory):
        exp = experimentFactory(_experiment_setup(withTrials=False))
        assert dict(exp.trialSet["Measurements"]) == {}

    def test_a_trial_set_without_trials_currently_fails_to_load(self, experimentFactory):
        """Characterisation of B163."""
        with pytest.raises(KeyError, match="trials"):
            experimentFactory(_experiment_setup(withTrials=False))

    @pytest.mark.xfail(
        strict=True,
        reason="B163: the same non-defensive indexing in "
               "EntityTypeWithData._initEntities (self._metadata['entities'] "
               "against argos's .get('entities', [])), reachable through a "
               "version-1.0.0 setup, which is passed through unmigrated. "
               "See the consolidated findings issue.",
    )
    def test_an_entity_type_without_entities_should_load_as_an_empty_entity_type(self, experimentFactory):
        exp = experimentFactory(_v1_setup_without_entities())
        assert dict(exp.entityType["sonic"]) == {}

    def test_an_entity_type_without_entities_currently_fails_to_load(self, experimentFactory):
        """Characterisation of B163, on the entity side."""
        with pytest.raises(KeyError, match="entities"):
            experimentFactory(_v1_setup_without_entities())


def _v1_setup_without_entities():
    """A version-1.0.0 setup (passed through unmigrated) whose entity type
    carries no ``entities`` key -- which argos's own ``_initEntities``
    tolerates."""
    return {
        "version": "1.0.0.",
        "name": "UNITEXP",
        "trialSets": [{"name": "Measurements", "description": "",
                       "attributeTypes": [], "trials": []}],
        "entityTypes": [{"name": "sonic", "attributeTypes": []}],
    }


@pytest.mark.unit
class TestDefaultTrialSet:
    def test_the_default_trial_set_name_is_the_one_given_at_construction(self, experiment):
        assert experiment.defaultTrialSet == "Measurements"

    def test_the_trials_of_the_default_trial_set_are_that_trial_set(self, experiment):
        assert experiment.trialsOfDefaultTrialSet is experiment.trialSet["Measurements"]

    def test_no_default_trial_set_name_leaves_the_property_none(self, experimentFactory):
        exp = experimentFactory(defaultTrialSetName=None)
        assert exp.defaultTrialSet is None

    def test_reading_the_trials_of_an_unset_default_trial_set_raises(self, experimentFactory):
        exp = experimentFactory(defaultTrialSetName=None)
        with pytest.raises(KeyError):
            exp.trialsOfDefaultTrialSet


# ---------------------------------------------------------------------------
# _initEntitiesTypes / EntityTypeWithData / EntityWithData construction
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestEntityTypeInitialisation:
    def test_every_entity_type_in_the_setup_becomes_an_entity_type_with_data(self, experiment):
        assert list(experiment.entityType) == ["sonic"]
        assert isinstance(experiment.entityType["sonic"], EntityTypeWithData)

    def test_the_entity_type_keeps_a_back_reference_to_the_experiment(self, experiment):
        assert experiment.entityType["sonic"].experiment is experiment

    def test_the_entity_type_shares_the_experiments_data_engine(self, experiment):
        assert experiment.entityType["sonic"]._experimentData is experiment.getExperimentData()

    def test_the_entities_are_entity_with_data_objects(self, experiment):
        entityType = experiment.entityType["sonic"]
        assert sorted(entityType) == ["sonic_1", "sonic_2"]
        assert isinstance(entityType["sonic_1"], EntityWithData)

    def test_every_entity_shares_the_same_data_engine(self, experiment):
        entity = experiment.entityType["sonic"]["sonic_1"]
        assert entity._experimentData is experiment.getExperimentData()

    def test_the_entity_properties_merge_constant_defaults_with_device_attributes(self, experiment):
        properties = experiment.entityType["sonic"]["sonic_2"].properties
        assert properties["StoreDataPerDevice"] is False
        assert properties["height"] == "6"

    def test_the_entity_reports_its_type_name(self, experiment):
        assert experiment.entityType["sonic"]["sonic_1"].entityType == "sonic"


# ---------------------------------------------------------------------------
# Reading data
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestTrialGetData:
    def test_it_defaults_to_the_trials_own_start_and_end_properties(self, experiment, deviceData):
        data = experiment.trialSet["Measurements"]["T1"].getData(deviceType="sonic")
        assert data.index.min() == pandas.Timestamp(TRIAL_START, tz=TZ)
        assert data.index.max() == pandas.Timestamp(TRIAL_END, tz=TZ)
        assert len(data) == 10

    def test_an_explicit_time_range_overrides_the_trial_properties(self, experiment, deviceData):
        data = experiment.trialSet["Measurements"]["T1"].getData(
            deviceType="sonic",
            startTime=pandas.Timestamp("2020-01-01 00:02:00", tz=TZ),
            endTime=pandas.Timestamp("2020-01-01 00:04:00", tz=TZ),
        )
        assert len(data) == 3

    def test_a_device_name_filters_the_rows_of_the_shared_file(self, experiment, deviceData):
        data = experiment.trialSet["Measurements"]["T1"].getData(
            deviceType="sonic", deviceName="sonic_1")
        assert set(data.deviceName) == {"sonic_1"}
        assert len(data) == 5

    def test_a_window_with_no_data_raises(self, experiment, deviceData):
        with pytest.raises(ValueError, match="There is no data for sonic"):
            experiment.trialSet["Measurements"]["T1"].getData(
                deviceType="sonic",
                startTime=pandas.Timestamp("2021-01-01 00:00:00", tz=TZ),
                endTime=pandas.Timestamp("2021-01-01 01:00:00", tz=TZ),
            )

    def test_an_unknown_device_type_has_no_documents_and_raises(self, experiment, deviceData):
        with pytest.raises(ValueError, match="There is no data for noSuchDevice"):
            experiment.trialSet["Measurements"]["T1"].getData(deviceType="noSuchDevice")


@pytest.mark.unit
class TestTrialGetDataWithMetadata:
    """B162: the merge key does not exist in argos's trial entities table."""

    def test_the_trial_entities_table_has_no_entity_name_column(self, experiment):
        """Characterisation of B162: this is the table the merge is given."""
        table = experiment.trialSet["Measurements"]["T1"].entitiesTable
        assert "deviceItemName" in table.columns
        assert "entityName" not in table.columns

    @pytest.mark.xfail(
        strict=True,
        reason="B162: TrialWithdata.getData merges on "
               "right_on='entityName', but argos's Trial.entitiesTable "
               "names the device column deviceItemName -- entityName exists "
               "only on Experiment.entitiesTable and "
               "EntityType.entitiesTable. See the consolidated findings issue.",
    )
    def test_it_should_join_the_trial_device_metadata_onto_the_data(self, experiment, deviceData):
        data = experiment.trialSet["Measurements"]["T1"].getData(
            deviceType="sonic", withMetadata=True)
        assert "height" in data.columns

    def test_it_currently_raises_a_key_error_on_the_merge_column(self, experiment, deviceData):
        """Characterisation of B162."""
        with pytest.raises(KeyError, match="entityName"):
            experiment.trialSet["Measurements"]["T1"].getData(
                deviceType="sonic", withMetadata=True)


@pytest.mark.unit
class TestEntityTypeGetData:
    def test_it_returns_every_row_of_the_device_type(self, experiment, deviceData):
        data = experiment.entityType["sonic"].getData()
        assert len(data) == 20

    def test_a_time_range_narrows_the_result(self, experiment, deviceData):
        data = experiment.entityType["sonic"].getData(
            startTime=pandas.Timestamp("2020-01-01 00:05:00", tz=TZ),
            endTime=pandas.Timestamp("2020-01-01 00:07:00", tz=TZ),
        )
        assert len(data) == 3


@pytest.mark.unit
class TestEntityTypeGetDataTrial:
    def test_it_uses_the_trials_start_and_end_properties(self, experiment, deviceData):
        data = experiment.entityType["sonic"].getDataTrial(
            trialSetName="Measurements", trialName="T1")
        assert len(data) == 10
        assert data.index.max() == pandas.Timestamp(TRIAL_END, tz=TZ)

    def test_an_unknown_trial_name_raises(self, experiment, deviceData):
        with pytest.raises(KeyError):
            experiment.entityType["sonic"].getDataTrial(
                trialSetName="Measurements", trialName="noSuchTrial")


@pytest.mark.unit
class TestEntityTypeGetDataTrialPerDevice:
    """B164: the per-device branch is never given a device name."""

    @pytest.fixture()
    def perDeviceExperiment(self, experimentFactory):
        return experimentFactory(_experiment_setup(storePerDevice=True))

    def test_the_per_device_flag_is_read_from_the_entity_type_attributes(self, perDeviceExperiment):
        """Characterisation of B164: the flag the method forwards is True."""
        properties = perDeviceExperiment.entityType["sonic"].properties
        flag = next(attr["defaultValue"] for attr in properties
                    if attr["name"] == "StoreDataPerDevice")
        assert flag is True

    @pytest.mark.xfail(
        strict=True,
        reason="B164: getDataTrial forwards perDevice=StoreDataPerDevice "
               "but never a deviceName, and parquetDataEngineHera.getData "
               "asserts deviceName in its per-device branch, so the call "
               "raises AssertionError for every entity type whose data is "
               "stored per device. See the consolidated findings issue.",
    )
    def test_it_should_read_per_device_data_for_the_trial(self, perDeviceExperiment, unit_project, tmp_path):
        _registerDeviceData(unit_project, tmp_path, _deviceFrame(), deviceName="sonic_1")
        data = perDeviceExperiment.entityType["sonic"].getDataTrial(
            trialSetName="Measurements", trialName="T1")
        assert len(data) > 0

    def test_it_currently_raises_an_assertion_error(self, perDeviceExperiment, unit_project, tmp_path):
        """Characterisation of B164."""
        _registerDeviceData(unit_project, tmp_path, _deviceFrame(), deviceName="sonic_1")
        with pytest.raises(AssertionError, match="deviceName should be defined"):
            perDeviceExperiment.entityType["sonic"].getDataTrial(
                trialSetName="Measurements", trialName="T1")

    def test_the_entity_level_read_does_pass_the_device_name(self, perDeviceExperiment, unit_project, tmp_path):
        """EntityWithData.getData forwards deviceName, so the same flag works there."""
        _registerDeviceData(unit_project, tmp_path, _deviceFrame(), deviceName="sonic_1")
        data = perDeviceExperiment.entityType["sonic"]["sonic_1"].getData()
        assert len(data) == 20


@pytest.mark.unit
class TestEntityGetData:
    def test_it_filters_the_shared_file_by_the_entity_name(self, experiment, deviceData):
        data = experiment.entityType["sonic"]["sonic_2"].getData()
        assert set(data.deviceName) == {"sonic_2"}
        assert len(data) == 10

    def test_a_time_range_narrows_the_result(self, experiment, deviceData):
        data = experiment.entityType["sonic"]["sonic_1"].getData(
            startTime=pandas.Timestamp("2020-01-01 00:00:00", tz=TZ),
            endTime=pandas.Timestamp("2020-01-01 00:05:00", tz=TZ),
        )
        assert len(data) == 3

    def test_an_entity_type_without_the_store_flag_raises(self, experimentFactory, deviceData):
        setup = _experiment_setup()
        setup["deviceTypes"][0]["attributeTypes"] = [
            attr for attr in setup["deviceTypes"][0]["attributeTypes"]
            if attr["name"] != "StoreDataPerDevice"
        ]
        exp = experimentFactory(setup)
        with pytest.raises(KeyError, match="StoreDataPerDevice"):
            exp.entityType["sonic"]["sonic_1"].getData()


@pytest.mark.unit
class TestGetDataFromDateRange:
    def test_it_returns_the_rows_in_the_range(self, experiment, deviceData):
        data = experiment.getDataFromDateRange(
            deviceType="sonic",
            startTime=pandas.Timestamp("2020-01-01 00:00:00", tz=TZ),
            endTime=pandas.Timestamp("2020-01-01 00:03:00", tz=TZ),
            withMetadata=False,
        )
        assert len(data) == 4

    def test_a_device_name_filters_the_rows(self, experiment, deviceData):
        data = experiment.getDataFromDateRange(
            deviceType="sonic",
            startTime=pandas.Timestamp("2020-01-01 00:00:00", tz=TZ),
            endTime=pandas.Timestamp("2020-01-01 00:19:00", tz=TZ),
            deviceName="sonic_1",
            withMetadata=False,
        )
        assert set(data.deviceName) == {"sonic_1"}

    def test_an_empty_range_raises(self, experiment, deviceData):
        with pytest.raises(ValueError, match="There is no data for sonic"):
            experiment.getDataFromDateRange(
                deviceType="sonic",
                startTime=pandas.Timestamp("2019-01-01 00:00:00", tz=TZ),
                endTime=pandas.Timestamp("2019-01-01 01:00:00", tz=TZ),
                withMetadata=False,
            )

    def test_with_metadata_joins_the_experiment_entities_table(self, experiment, deviceData):
        data = experiment.getDataFromDateRange(
            deviceType="sonic",
            startTime=pandas.Timestamp("2020-01-01 00:00:00", tz=TZ),
            endTime=pandas.Timestamp("2020-01-01 00:01:00", tz=TZ),
        )
        # Experiment.entitiesTable is long-form (one row per property), so
        # the merge fans each data row out over that entity's properties.
        assert set(data.entityName) == {"sonic_1", "sonic_2"}
        assert data.index.name == "timestamp"
        assert set(data["name"]) == {"StoreDataPerDevice", "height"}


# ---------------------------------------------------------------------------
# Coordinate conversion
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestProcessRow:
    """B160: convertCRS returns a list, not a GeoSeries."""

    @pytest.mark.xfail(
        strict=True,
        reason="B160: _process_row reads pp.x[0]/pp.y[0], but "
               "hera.measurements.GIS.utils.convertCRS returns "
               "list(gdf.to_crs(...).geometry) -- a plain list of shapely "
               "Points, which has no .x. The working access is pp[0].x. "
               "See the consolidated findings issue.",
    )
    def test_it_should_return_the_itm_coordinates_of_the_row(self, experiment):
        row = pandas.Series({"Longitude": 34.0, "Latitude": 32.0})
        converted = experiment._process_row(row)
        assert len(converted) == 2
        assert converted[0] > 100000

    def test_it_currently_raises_an_attribute_error(self, experiment):
        """Characterisation of B160."""
        row = pandas.Series({"Longitude": 34.0, "Latitude": 32.0})
        with pytest.raises(AttributeError, match="'list' object has no attribute 'x'"):
            experiment._process_row(row)


@pytest.mark.unit
class TestGetDevicesImageCoordinates:
    def test_it_returns_the_bounding_box_of_capitalised_coordinate_attributes(
            self, coordinateExperiment):
        bbox = coordinateExperiment.get_devices_image_coordinates(
            trialSetName="Measurements", trialName="T1",
            deviceType="sonic", outputCRS=WSG84,
        )
        assert bbox == (32.0, 34.0, 32.5, 34.5)

    def test_an_unmatched_device_type_yields_an_empty_selection_and_raises(
            self, coordinateExperiment):
        with pytest.raises(ValueError):
            coordinateExperiment.get_devices_image_coordinates(
                trialSetName="Measurements", trialName="T1",
                deviceType="noSuchDeviceType", outputCRS=WSG84,
            )


@pytest.mark.unit
class TestGetDevicesImageCoordinatesColumnNames:
    """B161: the method reads Latitude/Longitude; argos writes latitude/longitude."""

    def test_the_trial_entities_table_spells_the_coordinates_in_lower_case(self, experiment):
        """Characterisation of B161: this is the table the method reads."""
        table = experiment.trialSet["Measurements"]["T1"].entitiesTable
        assert {"latitude", "longitude"}.issubset(table.columns)
        assert "Latitude" not in table.columns

    @pytest.mark.xfail(
        strict=True,
        reason="B161: get_devices_image_coordinates reads the columns "
               "Latitude/Longitude, but argos spreads a device-on-trial "
               "location into lowercase latitude/longitude "
               "(fillContained.spread_attributes), so for map-positioned "
               "devices -- what every argos example export does -- the "
               "lookup raises KeyError. See the consolidated findings issue.",
    )
    def test_it_should_use_the_device_locations_from_the_map(self, experiment):
        bbox = experiment.get_devices_image_coordinates(
            trialSetName="Measurements", trialName="T1",
            deviceType="sonic", outputCRS=WSG84,
        )
        assert bbox == (32.0, 34.0, 32.5, 34.5)

    def test_it_currently_raises_a_key_error_on_the_capitalised_column(self, experiment):
        """Characterisation of B161."""
        with pytest.raises(KeyError, match="Latitude"):
            experiment.get_devices_image_coordinates(
                trialSetName="Measurements", trialName="T1",
                deviceType="sonic", outputCRS=WSG84,
            )

    def test_the_itm_branch_is_dead_even_with_capitalised_columns(self, coordinateExperiment):
        """Characterisation of B160, reached through the public method."""
        with pytest.raises(AttributeError, match="'list' object has no attribute 'x'"):
            coordinateExperiment.get_devices_image_coordinates(
                trialSetName="Measurements", trialName="T1",
                deviceType="sonic", outputCRS=ITM,
            )
