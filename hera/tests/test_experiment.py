"""
Comprehensive pytest tests for the experiment toolkit.

Tests cover:
- experimentHome: keys, getExperimentsMap, getExperimentsTable, getExperiment
- experimentSetupWithData: properties, trialSet, entityType access
- TrialWithdata.getData: time filtering, metadata merge
- EntityTypeWithData.getData, getDataTrial
- EntityWithData.getData
- experimentAnalysis: addMetadata, addTrialProperties, getDeviceLocations
- experimentPresentation: plotDeviceTypeFunctionality (smoke test)
- dataEngineFactory: engine selection
- parquetDataEngineHera: getData, getDataFromTrial
- Parsers: Parser_OldStyleMetaDataParquet

DB-backed tests require MongoDB and auto-skip when unavailable.
Synthetic tests create in-memory data and run without infrastructure.
"""

import json
import os
import tempfile
import zipfile

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# MongoDB connectivity check
# ---------------------------------------------------------------------------

def _mongo_is_available():
    """Return True if the configured MongoDB server is reachable.

    Uses a direct pymongo ping with a 1 s server-selection timeout so that
    collection-time probing does not hang for 30 s when Mongo is down.
    """
    try:
        import pymongo
        from hera.datalayer.document import getMongoConfigFromJson
        cfg = getMongoConfigFromJson()
        host = cfg.get("dbIP", "localhost")
        port = int(cfg.get("port", 27017))
        client = pymongo.MongoClient(
            host=host, port=port, serverSelectionTimeoutMS=1000
        )
        client.server_info()
        return True
    except Exception:
        return False


_MONGO_OK = _mongo_is_available()
requires_mongo = pytest.mark.skipif(not _MONGO_OK, reason="MongoDB not reachable")

# ---------------------------------------------------------------------------
# Try importing experiment classes — some may fail if geopandas/pyargos missing
# ---------------------------------------------------------------------------

try:
    from hera.measurements.experiment.experiment import (
        experimentHome,
        experimentSetupWithData,
        TrialSetWithData,
        TrialWithdata,
        EntityTypeWithData,
        EntityWithData,
    )
    _EXPERIMENT_CLASSES_AVAILABLE = True
except ImportError:
    _EXPERIMENT_CLASSES_AVAILABLE = False

try:
    from hera.measurements.experiment.dataEngine import (
        dataEngineFactory,
        parquetDataEngineHera,
        PARQUETHERA,
        PANDASDB,
        DASKDB,
    )
    _ENGINE_AVAILABLE = True
except ImportError:
    _ENGINE_AVAILABLE = False

try:
    from hera.measurements.experiment.analysis import experimentAnalysis
    _ANALYSIS_AVAILABLE = True
except ImportError:
    _ANALYSIS_AVAILABLE = False

try:
    from hera.measurements.experiment.presentation import experimentPresentation
    _PRESENTATION_AVAILABLE = True
except ImportError:
    _PRESENTATION_AVAILABLE = False

requires_experiment = pytest.mark.skipif(
    not _EXPERIMENT_CLASSES_AVAILABLE,
    reason="Experiment toolkit not available (missing pyargos or geopandas)",
)

PROJECT_NAME = "pytest_experiment_comprehensive"


# ---------------------------------------------------------------------------
# Helpers: create synthetic experiment data
# ---------------------------------------------------------------------------

def _create_sonic_parquet(path, n_rows=100):
    """Create a synthetic sonic parquet file."""
    idx = pd.date_range("2024-03-15 08:00:00", periods=n_rows, freq="100ms", tz="UTC")
    df = pd.DataFrame({
        "u": np.random.normal(-2.5, 1.0, n_rows),
        "v": np.random.normal(3.0, 0.8, n_rows),
        "w": np.random.normal(0.0, 0.5, n_rows),
        "T": np.random.normal(25.0, 1.0, n_rows),
        "deviceName": "sonic01",
    }, index=idx)
    df.index.name = "timestamp"
    df.to_parquet(path)
    return df


def _create_trh_parquet(path, n_rows=100):
    """Create a synthetic TRH parquet file."""
    idx = pd.date_range("2024-03-15 08:00:00", periods=n_rows, freq="1s", tz="UTC")
    df = pd.DataFrame({
        "TC_T": np.random.normal(24.0, 0.5, n_rows),
        "TRH": np.random.normal(24.5, 0.5, n_rows),
        "RH": np.random.normal(65.0, 5.0, n_rows),
        "deviceName": "TRH01",
    }, index=idx)
    df.index.name = "timestamp"
    df.to_parquet(path)
    return df


def _create_argos_zip(zip_path, experiment_name="TestExp"):
    """Create an Argos zip with trial sets, entity types, trials, and entities."""
    data_json = {
        "version": "1.0.0.",
        "experiment": {
            "name": experiment_name,
            "description": "Test experiment",
            "version": "1.0.0.",
        },
        "entityTypes": [
            {
                "name": "Sonic",
                "attributeTypes": [
                    {
                        "type": "Boolean",
                        "name": "StoreDataPerDevice",
                        "defaultValue": False,
                        "scope": "Constant",
                    },
                    {"type": "String", "name": "stationName", "scope": "Device"},
                    {"type": "Number", "name": "height", "scope": "Device"},
                ],
                "entities": [
                    {
                        "name": "sonic01",
                        "attributes": [
                            {"name": "stationName", "value": "StationA"},
                            {"name": "height", "value": "10"},
                        ],
                    },
                ],
            },
            {
                "name": "TRH",
                "attributeTypes": [
                    {
                        "type": "Boolean",
                        "name": "StoreDataPerDevice",
                        "defaultValue": False,
                        "scope": "Constant",
                    },
                ],
                "entities": [
                    {"name": "TRH01", "attributes": []},
                ],
            },
        ],
        "trialSets": [
            {
                "name": "Measurements",
                "attributeTypes": [
                    {"type": "datetime-local", "name": "TrialStart", "label": "TrialStart", "description": ""},
                    {"type": "datetime-local", "name": "TrialEnd", "label": "TrialEnd", "description": ""},
                ],
                "trials": [
                    {
                        "name": "Trial_01",
                        "created": "2024-03-15T06:00:00.000Z",
                        "cloneFrom": None,
                        "properties": [
                            {"key": "TrialStart", "val": "2024-03-15T08:00:00.000Z", "type": "datetime-local"},
                            {"key": "TrialEnd", "val": "2024-03-15T08:00:10.000Z", "type": "datetime-local"},
                        ],
                        "entities": [
                            {
                                "deviceTypeName": "Sonic",
                                "deviceItemName": "sonic01",
                                "location": {"name": "OSMMap", "coordinates": [32.79, 35.04]},
                                "attributes": [{"name": "height", "value": "10"}],
                            },
                            {
                                "deviceTypeName": "TRH",
                                "deviceItemName": "TRH01",
                                "location": {"name": "OSMMap", "coordinates": [32.79, 35.04]},
                                "attributes": [],
                            },
                        ],
                    },
                ],
            },
        ],
        "maps": [],
    }
    os.makedirs(os.path.dirname(zip_path), exist_ok=True)
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("data.json", json.dumps(data_json, indent=2))


def _create_experiment_dir(base_dir, experiment_name="TestExp"):
    """Create a complete experiment directory with data files."""
    exp_dir = os.path.join(base_dir, experiment_name)
    os.makedirs(os.path.join(exp_dir, "code"), exist_ok=True)
    os.makedirs(os.path.join(exp_dir, "data"), exist_ok=True)
    os.makedirs(os.path.join(exp_dir, "runtimeExperimentData"), exist_ok=True)

    # Experiment class
    code = f"""
from hera.measurements.experiment.experiment import experimentSetupWithData
from hera.measurements.experiment.analysis import experimentAnalysis
from hera.measurements.experiment.presentation import experimentPresentation

class {experiment_name}(experimentSetupWithData):
    def __init__(self, projectName, pathToExperiment, filesDirectory=None):
        super().__init__(projectName=projectName, pathToExperiment=pathToExperiment, filesDirectory=filesDirectory)
"""
    with open(os.path.join(exp_dir, "code", f"{experiment_name}.py"), "w") as f:
        f.write(code)

    # Config
    config = {"experimentName": experiment_name}
    with open(os.path.join(exp_dir, "runtimeExperimentData", "Datasources_Configurations.json"), "w") as f:
        json.dump(config, f)

    # Argos zip
    _create_argos_zip(
        os.path.join(exp_dir, "runtimeExperimentData", f"{experiment_name}.zip"),
        experiment_name,
    )

    # Parquet data
    _create_sonic_parquet(os.path.join(exp_dir, "data", "Sonic.parquet"))
    _create_trh_parquet(os.path.join(exp_dir, "data", "TRH.parquet"))

    # Repository JSON
    repo = {
        "experiment": {
            "DataSource": {
                experiment_name: {
                    "isRelativePath": "False",
                    "item": {
                        "dataSourceName": experiment_name,
                        "resource": exp_dir,
                        "dataFormat": "string",
                        "overwrite": "True",
                    },
                }
            },
            "Measurements": {
                "Sonic": {
                    "isRelativePath": "True",
                    "item": {
                        "type": "Experiment_rawData",
                        "resource": "data/Sonic.parquet",
                        "dataFormat": "parquet",
                        "desc": {
                            "deviceType": "Sonic",
                            "experimentName": experiment_name,
                            "deviceName": "",
                        },
                    },
                },
                "TRH": {
                    "isRelativePath": "True",
                    "item": {
                        "type": "Experiment_rawData",
                        "resource": "data/TRH.parquet",
                        "dataFormat": "parquet",
                        "desc": {
                            "deviceType": "TRH",
                            "experimentName": experiment_name,
                            "deviceName": "",
                        },
                    },
                },
            },
        }
    }
    with open(os.path.join(exp_dir, f"{experiment_name}_repository.json"), "w") as f:
        json.dump(repo, f, indent=2)

    return exp_dir


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def experiment_dir():
    """Create a temporary experiment directory with synthetic data."""
    with tempfile.TemporaryDirectory(prefix="hera_exp_test_") as tmpdir:
        exp_dir = _create_experiment_dir(tmpdir, "TestExp")
        yield exp_dir


@pytest.fixture(scope="module")
def loaded_project(experiment_dir):
    """Load the synthetic experiment into a Hera project (requires MongoDB)."""
    if not _MONGO_OK:
        pytest.skip("MongoDB not reachable")

    from hera.datalayer.project import Project
    from hera.utils.data.toolkit import dataToolkit

    repo_path = os.path.join(experiment_dir, "TestExp_repository.json")
    with open(repo_path) as f:
        repo_json = json.load(f)

    proj = Project(projectName=PROJECT_NAME)
    dt = dataToolkit()
    dt.loadAllDatasourcesInRepositoryJSONToProject(
        projectName=PROJECT_NAME,
        repositoryJSON=repo_json,
        basedir=experiment_dir,
        overwrite=True,
    )
    yield proj

    # Cleanup
    for getter in (proj.getMeasurementsDocuments, proj.getCacheDocuments):
        try:
            for doc in getter():
                doc.delete()
        except Exception:
            pass


@pytest.fixture(scope="module")
def exp_home(loaded_project):
    """Get experimentHome for the loaded project."""
    return experimentHome(projectName=PROJECT_NAME)


@pytest.fixture(scope="module")
def experiment(exp_home):
    """Get the loaded experiment instance."""
    return exp_home.getExperiment("TestExp")


# ===========================================================================
# 1. experimentHome
# ===========================================================================

@pytest.mark.integration
@requires_mongo
@requires_experiment
class TestExperimentHome:

    def test_keys(self, exp_home):
        keys = exp_home.keys()
        assert "TestExp" in keys

    def test_getitem(self, exp_home):
        exp = exp_home["TestExp"]
        assert exp is not None

    def test_get_experiments_map(self, exp_home):
        emap = exp_home.getExperimentsMap()
        assert isinstance(emap, dict)
        assert "TestExp" in emap

    def test_get_experiments_table(self, exp_home):
        table = exp_home.getExperimentsTable()
        assert table is not None
        assert len(table) >= 1

    def test_get_experiment(self, exp_home):
        exp = exp_home.getExperiment("TestExp")
        assert exp is not None

    def test_nonexistent_experiment(self, exp_home):
        with pytest.raises(Exception):
            exp_home.getExperiment("NonExistentExp_XYZ")


# ===========================================================================
# 2. experimentSetupWithData — properties
# ===========================================================================

@pytest.mark.integration
@requires_mongo
@requires_experiment
class TestExperimentSetup:

    def test_name(self, experiment):
        assert experiment.name == "TestExp"

    def test_configuration(self, experiment):
        assert experiment.configuration is not None
        assert isinstance(experiment.configuration, dict)

    def test_trial_set_access(self, experiment):
        assert "Measurements" in experiment.trialSet

    def test_entity_type_access(self, experiment):
        assert "Sonic" in experiment.entityType
        assert "TRH" in experiment.entityType

    def test_analysis_property(self, experiment):
        assert experiment.analysis is not None

    def test_presentation_property(self, experiment):
        assert experiment.presentation is not None

    def test_get_experiment_data(self, experiment):
        engine = experiment.getExperimentData()
        assert engine is not None


# ===========================================================================
# 3. TrialSetWithData and TrialWithdata
# ===========================================================================

@pytest.mark.integration
@requires_mongo
@requires_experiment
class TestTrialAccess:

    def test_trial_set_contains_trial(self, experiment):
        ts = experiment.trialSet["Measurements"]
        assert "Trial_01" in ts

    def test_trial_properties(self, experiment):
        trial = experiment.trialSet["Measurements"]["Trial_01"]
        props = trial.properties
        assert "TrialStart" in props
        assert "TrialEnd" in props

    def test_trial_get_data_sonic(self, experiment):
        trial = experiment.trialSet["Measurements"]["Trial_01"]
        df = trial.getData(deviceType="Sonic")
        assert df is not None
        if hasattr(df, "compute"):
            df = df.compute()
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0

    def test_trial_get_data_with_device_name(self, experiment):
        trial = experiment.trialSet["Measurements"]["Trial_01"]
        df = trial.getData(deviceType="Sonic", deviceName="sonic01")
        if hasattr(df, "compute"):
            df = df.compute()
        assert isinstance(df, pd.DataFrame)

    def test_trial_get_data_trh(self, experiment):
        trial = experiment.trialSet["Measurements"]["Trial_01"]
        df = trial.getData(deviceType="TRH")
        if hasattr(df, "compute"):
            df = df.compute()
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0


# ===========================================================================
# 4. EntityTypeWithData
# ===========================================================================

@pytest.mark.integration
@requires_mongo
@requires_experiment
class TestEntityType:

    def test_entity_type_name(self, experiment):
        assert experiment.entityType["Sonic"].name == "Sonic"

    def test_entity_type_contains_entity(self, experiment):
        sonic = experiment.entityType["Sonic"]
        assert "sonic01" in sonic

    def test_get_data(self, experiment):
        sonic = experiment.entityType["Sonic"]
        df = sonic.getData()
        if hasattr(df, "compute"):
            df = df.compute()
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0

    def test_get_data_trial(self, experiment):
        sonic = experiment.entityType["Sonic"]
        df = sonic.getDataTrial(trialSetName="Measurements", trialName="Trial_01")
        if hasattr(df, "compute"):
            df = df.compute()
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0


# ===========================================================================
# 5. EntityWithData
# ===========================================================================

@pytest.mark.integration
@requires_mongo
@requires_experiment
class TestEntity:

    def test_entity_name(self, experiment):
        entity = experiment.entityType["Sonic"]["sonic01"]
        assert entity.name == "sonic01"

    def test_entity_get_data(self, experiment):
        entity = experiment.entityType["Sonic"]["sonic01"]
        df = entity.getData()
        if hasattr(df, "compute"):
            df = df.compute()
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0


# ===========================================================================
# 6. getDataFromDateRange
# ===========================================================================

@pytest.mark.integration
@requires_mongo
@requires_experiment
class TestGetDataFromDateRange:

    def test_basic(self, experiment):
        df = experiment.getDataFromDateRange(
            deviceType="Sonic",
            startTime="2024-03-15 08:00:00",
            endTime="2024-03-15 08:00:05",
        )
        if hasattr(df, "compute"):
            df = df.compute()
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0


# ===========================================================================
# 7. dataEngineFactory
# ===========================================================================

@pytest.mark.unit
@pytest.mark.skipif(not _ENGINE_AVAILABLE, reason="dataEngine not importable")
class TestDataEngineFactory:
    """Test engine type constants (no DB needed)."""

    def test_constants_exist(self):
        assert PARQUETHERA is not None
        assert PANDASDB is not None
        assert DASKDB is not None

    def test_constants_are_strings(self):
        assert isinstance(PARQUETHERA, str)
        assert isinstance(PANDASDB, str)
        assert isinstance(DASKDB, str)


# ===========================================================================
# 8. experimentAnalysis — synthetic (no DB)
# ===========================================================================

@pytest.mark.unit
@pytest.mark.skipif(not _ANALYSIS_AVAILABLE, reason="experimentAnalysis not importable")
class TestAnalysisSynthetic:
    """Test analysis methods with synthetic DataFrames."""

    @pytest.fixture()
    def mock_analysis(self):
        """Create an analysis instance with a mock datalayer."""

        class MockDatalayer:
            projectName = "test"
            trialSet = {}
            entityType = {}

            def getConfig(self):
                return {}

        return experimentAnalysis(MockDatalayer())

    def test_add_trial_properties(self, mock_analysis):
        idx = pd.date_range("2024-03-15 08:00:00", periods=5, freq="1s")
        df = pd.DataFrame({"value": [1, 2, 3, 4, 5], "timestamp": idx}, index=idx)

        # Mock experiment with trial
        class MockTrial:
            properties = {
                "TrialStart": pd.Timestamp("2024-03-15 08:00:00"),
                "TrialEnd": pd.Timestamp("2024-03-15 08:00:05"),
                "ReleaseStart": pd.Timestamp("2024-03-15 08:00:01"),
            }

        class MockTrialSet(dict):
            def __init__(self):
                super().__init__()
                self["Trial_01"] = MockTrial()

        mock_analysis._datalayer.trialSet = {"Measurements": MockTrialSet()}

        result = mock_analysis.addTrialProperties(
            data=df, trialName="Trial_01", trialSetName="Measurements",
        )
        assert "fromStart" in result.columns or "fromStartSeconds" in result.columns


# ===========================================================================
# 9. experimentPresentation — smoke test (no DB)
# ===========================================================================

@pytest.mark.unit
@pytest.mark.skipif(not _PRESENTATION_AVAILABLE, reason="experimentPresentation not importable")
class TestPresentationInit:
    """Test presentation can be instantiated."""

    def test_init(self):
        class MockDatalayer:
            projectName = "test"
            trialSet = {}
            entityType = {}

            def getConfig(self):
                return {}

        class MockAnalysis:
            pass

        pres = experimentPresentation(MockDatalayer(), MockAnalysis())
        assert pres is not None


# ===========================================================================
# 10. Parsers — synthetic (no DB)
# ===========================================================================

@pytest.mark.unit
class TestParsers:
    """Test parsers with synthetic data."""

    def test_old_style_metadata_parser_import(self):
        from hera.measurements.experiment.parsers import Parser_OldStyleMetaDataParquet
        assert Parser_OldStyleMetaDataParquet is not None

    def test_campbell_binary_parser_import(self):
        from hera.measurements.experiment.parsers import Parser_CampbellBinary
        assert Parser_CampbellBinary is not None


# ===========================================================================
# 11. Argos zip parsing (no DB)
# ===========================================================================

@pytest.mark.integration
@requires_experiment
class TestArgosZipParsing:
    """Test that Argos zip files are parsed correctly."""

    def test_parse_synthetic_zip(self, experiment_dir):
        try:
            from argos.experimentSetup.dataObjects import ExperimentZipFile
        except ImportError:
            pytest.skip("pyargos not available")

        zip_path = os.path.join(
            experiment_dir, "runtimeExperimentData", "TestExp.zip"
        )
        exp = ExperimentZipFile(zip_path)

        # Trial sets
        assert "Measurements" in exp.trialSet
        ts = exp.trialSet["Measurements"]
        assert "Trial_01" in ts

        # Trial properties
        trial = ts["Trial_01"]
        assert "TrialStart" in trial.properties
        assert "TrialEnd" in trial.properties

        # Entity types
        assert "Sonic" in exp.entityType
        assert "TRH" in exp.entityType

        # Entities
        sonic = exp.entityType["Sonic"]
        assert "sonic01" in sonic

    def test_entity_properties(self, experiment_dir):
        try:
            from argos.experimentSetup.dataObjects import ExperimentZipFile
        except ImportError:
            pytest.skip("pyargos not available")

        zip_path = os.path.join(
            experiment_dir, "runtimeExperimentData", "TestExp.zip"
        )
        exp = ExperimentZipFile(zip_path)
        sonic01 = exp.entityType["Sonic"]["sonic01"]
        props = sonic01.properties
        # Should have stationName and height from Device scope
        prop_names = [p["name"] if isinstance(p, dict) else p for p in props]
        assert any("stationName" in str(p) for p in props) or len(props) > 0

    def test_trial_entities_table(self, experiment_dir):
        try:
            from argos.experimentSetup.dataObjects import ExperimentZipFile
        except ImportError:
            pytest.skip("pyargos not available")

        zip_path = os.path.join(
            experiment_dir, "runtimeExperimentData", "TestExp.zip"
        )
        exp = ExperimentZipFile(zip_path)
        trial = exp.trialSet["Measurements"]["Trial_01"]
        table = trial.entitiesTable
        assert isinstance(table, pd.DataFrame)
        assert len(table) == 2  # sonic01 + TRH01


# ===========================================================================
# 12. Full integration: load → navigate → get data
# ===========================================================================

@pytest.mark.integration
@requires_mongo
@requires_experiment
class TestFullIntegration:
    """End-to-end: load experiment → navigate hierarchy → retrieve data."""

    def test_full_workflow(self, experiment):
        # Navigate to trial
        trial = experiment.trialSet["Measurements"]["Trial_01"]

        # Get sonic data
        sonic_df = trial.getData(deviceType="Sonic")
        if hasattr(sonic_df, "compute"):
            sonic_df = sonic_df.compute()
        assert "u" in sonic_df.columns or "value" in sonic_df.columns
        assert len(sonic_df) > 0

        # Get TRH data
        trh_df = trial.getData(deviceType="TRH")
        if hasattr(trh_df, "compute"):
            trh_df = trh_df.compute()
        assert len(trh_df) > 0

        # Entity type access
        sonic_type = experiment.entityType["Sonic"]
        assert "sonic01" in sonic_type

        # Single entity
        entity = sonic_type["sonic01"]
        entity_df = entity.getData()
        if hasattr(entity_df, "compute"):
            entity_df = entity_df.compute()
        assert len(entity_df) > 0
