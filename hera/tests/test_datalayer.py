"""
Comprehensive pytest tests for hera.datalayer.

Tests hit a real MongoDB instance and clean up after themselves.
Covers: Project CRUD (all 3 collections), counters, configuration,
saveData family, getDocumentByID, getMetadata, getAllDocuments,
addDocumentFromDict, nested queries, Collection direct access,
DataHandlers round-trip, autocache decorator, MetadataFrame,
nonDBMetadataFrame, and export/import.
"""

import json
import os
import tempfile
import shutil

import numpy as np
import pandas as pd
import pytest

from hera.datalayer.project import Project, getProjectList, createProjectDirectory
from hera.datalayer.collection import (
    AbstractCollection,
    Measurements_Collection,
    Simulations_Collection,
    Cache_Collection,
)
from hera.datalayer.datahandler import datatypes, getHandler, guessHandler
from hera.datalayer.document.metadataDocument import nonDBMetadataFrame


# ---------------------------------------------------------------------------
# MongoDB connectivity check
# ---------------------------------------------------------------------------

def _mongo_is_available():
    """Return True if the configured MongoDB server is reachable."""
    try:
        p = Project(projectName="__ping__")
        list(p.getMeasurementsDocuments())
        return True
    except Exception:
        return False


_MONGO_OK = _mongo_is_available()
requires_mongo = pytest.mark.skipif(not _MONGO_OK, reason="MongoDB not reachable")


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PROJECT_NAME = "pytest_datalayer_comprehensive"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def proj():
    """Module-scoped Project with full cleanup."""
    p = Project(projectName=PROJECT_NAME)
    yield p
    # Cleanup: remove all documents from all collections
    for getter in (p.getMeasurementsDocuments, p.getSimulationsDocuments, p.getCacheDocuments):
        try:
            for doc in getter():
                doc.delete()
        except Exception:
            pass
    # Clean config
    try:
        p.setConfig(keep_old_values=False)
    except Exception:
        pass


@pytest.fixture(scope="module")
def tmp_dir():
    """Temporary directory for file-based tests."""
    d = tempfile.mkdtemp(prefix="hera_test_")
    yield d
    shutil.rmtree(d, ignore_errors=True)


@pytest.fixture(scope="function")
def clean_proj():
    """Function-scoped Project that cleans up after each test."""
    name = "pytest_datalayer_func"
    p = Project(projectName=name)
    yield p
    for getter in (p.getMeasurementsDocuments, p.getSimulationsDocuments, p.getCacheDocuments):
        try:
            for doc in getter():
                doc.delete()
        except Exception:
            pass


# ===========================================================================
# 1. Project CRUD — all three collections
# ===========================================================================

@requires_mongo
class TestMeasurementsCRUD:
    """Add, get, get-as-dict, delete for Measurements collection."""

    def test_add_and_get(self, proj):
        doc = proj.addMeasurementsDocument(
            resource="/tmp/test_meas.csv",
            dataFormat="string",
            type="TestType",
            desc={"sensor": "A", "location": "roof"},
        )
        assert doc is not None

        docs = list(proj.getMeasurementsDocuments(type="TestType", sensor="A"))
        assert len(docs) >= 1
        assert any(d.type == "TestType" for d in docs)

    def test_get_as_dict(self, proj):
        proj.addMeasurementsDocument(
            resource="/tmp/dict_test.csv",
            dataFormat="string",
            type="DictTest",
            desc={"key": "value"},
        )
        result = proj.getMeasurementsDocumentsAsDict(type="DictTest")
        assert "documents" in result
        assert len(result["documents"]) >= 1

    def test_get_as_dict_with_id(self, proj):
        proj.addMeasurementsDocument(
            resource="/tmp/dict_id.csv",
            dataFormat="string",
            type="DictIdTest",
            desc={},
        )
        result = proj.getMeasurementsDocumentsAsDict(with_id=True, type="DictIdTest")
        assert len(result["documents"]) >= 1
        assert "_id" in result["documents"][0]

    def test_delete(self, proj):
        proj.addMeasurementsDocument(
            resource="/tmp/to_delete.csv",
            dataFormat="string",
            type="DeleteMe",
            desc={},
        )
        proj.deleteMeasurementsDocuments(type="DeleteMe")
        remaining = list(proj.getMeasurementsDocuments(type="DeleteMe"))
        assert len(remaining) == 0


@requires_mongo
class TestSimulationsCRUD:
    """Add, get, get-as-dict, delete for Simulations collection."""

    def test_add_and_get(self, proj):
        doc = proj.addSimulationsDocument(
            resource="/tmp/sim_output.nc",
            dataFormat="string",
            type="CFD_Run",
            desc={"solver": "OpenFOAM", "case": "windTunnel"},
        )
        assert doc is not None

        docs = list(proj.getSimulationsDocuments(type="CFD_Run", solver="OpenFOAM"))
        assert len(docs) >= 1

    def test_get_as_dict(self, proj):
        proj.addSimulationsDocument(
            resource="/tmp/sim_dict.nc",
            dataFormat="string",
            type="SimDict",
            desc={"run": 1},
        )
        result = proj.getSimulationsDocumentsAsDict(type="SimDict")
        assert "documents" in result
        assert len(result["documents"]) >= 1

    def test_delete(self, proj):
        proj.addSimulationsDocument(
            resource="/tmp/sim_del.nc",
            dataFormat="string",
            type="SimDelete",
            desc={},
        )
        proj.deleteSimulationsDocuments(type="SimDelete")
        remaining = list(proj.getSimulationsDocuments(type="SimDelete"))
        assert len(remaining) == 0


@requires_mongo
class TestCacheCRUD:
    """Add, get, get-as-dict, delete for Cache collection."""

    def test_add_and_get(self, proj):
        doc = proj.addCacheDocument(
            resource="/tmp/cache_result.pkl",
            dataFormat="string",
            type="ProcessedStats",
            desc={"function": "calcMean", "version": 2},
        )
        assert doc is not None

        docs = list(proj.getCacheDocuments(type="ProcessedStats", function="calcMean"))
        assert len(docs) >= 1

    def test_get_as_dict(self, proj):
        proj.addCacheDocument(
            resource="/tmp/cache_dict.pkl",
            dataFormat="string",
            type="CacheDict",
            desc={"step": "final"},
        )
        result = proj.getCacheDocumentsAsDict(type="CacheDict")
        assert "documents" in result
        assert len(result["documents"]) >= 1

    def test_delete(self, proj):
        proj.addCacheDocument(
            resource="/tmp/cache_del.pkl",
            dataFormat="string",
            type="CacheDel",
            desc={},
        )
        proj.deleteCacheDocuments(type="CacheDel")
        remaining = list(proj.getCacheDocuments(type="CacheDel"))
        assert len(remaining) == 0


# ===========================================================================
# 2. Nested queries (dictToMongoQuery)
# ===========================================================================

@requires_mongo
class TestNestedQueries:
    """Test querying with nested desc fields."""

    def test_single_level_query(self, clean_proj):
        clean_proj.addMeasurementsDocument(
            resource="r1", dataFormat="string", type="Nested",
            desc={"station": "Haifa", "elevation": 120},
        )
        clean_proj.addMeasurementsDocument(
            resource="r2", dataFormat="string", type="Nested",
            desc={"station": "TelAviv", "elevation": 5},
        )
        docs = list(clean_proj.getMeasurementsDocuments(type="Nested", station="Haifa"))
        assert len(docs) == 1
        assert docs[0].desc["station"] == "Haifa"

    def test_multi_level_nested_query(self, clean_proj):
        clean_proj.addMeasurementsDocument(
            resource="r3", dataFormat="string", type="DeepNested",
            desc={"config": {"model": "LSM", "resolution": 30}},
        )
        clean_proj.addMeasurementsDocument(
            resource="r4", dataFormat="string", type="DeepNested",
            desc={"config": {"model": "OpenFOAM", "resolution": 10}},
        )
        docs = list(clean_proj.getMeasurementsDocuments(
            type="DeepNested", config={"model": "LSM"}
        ))
        assert len(docs) == 1
        assert docs[0].desc["config"]["model"] == "LSM"


# ===========================================================================
# 3. getAllDocuments — cross-collection query
# ===========================================================================

@requires_mongo
class TestGetAllDocuments:
    """Test getAllDocuments aggregates across collections."""

    def test_returns_from_all_collections(self, clean_proj):
        clean_proj.addMeasurementsDocument(
            resource="all_m", dataFormat="string", type="CrossColl", desc={},
        )
        clean_proj.addSimulationsDocument(
            resource="all_s", dataFormat="string", type="CrossColl", desc={},
        )
        clean_proj.addCacheDocument(
            resource="all_c", dataFormat="string", type="CrossColl", desc={},
        )
        docs = list(clean_proj.getAllDocuments(type="CrossColl"))
        assert len(docs) == 3


# ===========================================================================
# 4. getDocumentByID
# ===========================================================================

@requires_mongo
class TestGetDocumentByID:
    """Test retrieving a document by its MongoDB ID."""

    def test_round_trip(self, clean_proj):
        doc = clean_proj.addMeasurementsDocument(
            resource="id_test", dataFormat="string", type="IDTest", desc={"tag": "findme"},
        )
        doc_id = str(doc.id)
        retrieved = clean_proj.getDocumentByID(doc_id)
        assert retrieved is not None
        assert retrieved.desc["tag"] == "findme"


# ===========================================================================
# 5. getMetadata
# ===========================================================================

@requires_mongo
class TestGetMetadata:
    """Test getMetadata returns DataFrame of document descriptions."""

    def test_returns_dataframe(self, clean_proj):
        clean_proj.addMeasurementsDocument(
            resource="meta1", dataFormat="string", type="MetaTest",
            desc={"station": "A", "height": 10},
        )
        clean_proj.addMeasurementsDocument(
            resource="meta2", dataFormat="string", type="MetaTest",
            desc={"station": "B", "height": 20},
        )
        metadata = clean_proj.getMetadata()
        assert isinstance(metadata, pd.DataFrame)
        assert len(metadata) >= 2


# ===========================================================================
# 6. addDocumentFromDict
# ===========================================================================

@requires_mongo
class TestAddDocumentFromDict:
    """Test adding documents from a dictionary representation."""

    def test_add_measurements_from_dict(self, clean_proj):
        doc_dict = {
            "_cls": "Metadata.Measurements",
            "projectName": clean_proj.projectName,
            "desc": {"fromDict": True, "value": 42},
            "type": "FromDict",
            "resource": "/tmp/from_dict",
            "dataFormat": "string",
        }
        clean_proj.addDocumentFromDict(doc_dict)
        docs = list(clean_proj.getMeasurementsDocuments(type="FromDict"))
        assert len(docs) >= 1
        assert docs[0].desc["fromDict"] is True


# ===========================================================================
# 7. Counters
# ===========================================================================

@requires_mongo
class TestCounters:
    """Test counter operations: define, set, get, getAndAdd."""

    def test_define_counter_new(self, clean_proj):
        result = clean_proj.defineCounter("test_counter_1", defaultValue=0)
        assert result is True

    def test_define_counter_existing(self, clean_proj):
        clean_proj.defineCounter("test_counter_2", defaultValue=5)
        result = clean_proj.defineCounter("test_counter_2", defaultValue=99)
        assert result is False
        assert clean_proj.getCounter("test_counter_2") == 5

    def test_set_and_get_counter(self, clean_proj):
        clean_proj.setCounter("explicit_counter", defaultValue=42)
        val = clean_proj.getCounter("explicit_counter")
        assert val == 42

    def test_get_counter_nonexistent(self, clean_proj):
        val = clean_proj.getCounter("nonexistent_counter_xyz")
        assert val is None

    def test_get_counter_and_add(self, clean_proj):
        # First call initializes to 0, adds 1, returns 0
        val1 = clean_proj.getCounterAndAdd("incr_counter", addition=1)
        # Should be 0 on first access
        assert isinstance(val1, (int, float))

        # Second call should be incremented
        val2 = clean_proj.getCounterAndAdd("incr_counter", addition=5)
        assert val2 > val1

    def test_set_counter_overwrites(self, clean_proj):
        clean_proj.setCounter("overwrite_ctr", defaultValue=10)
        clean_proj.setCounter("overwrite_ctr", defaultValue=99)
        assert clean_proj.getCounter("overwrite_ctr") == 99


# ===========================================================================
# 8. Configuration
# ===========================================================================

@requires_mongo
class TestConfiguration:
    """Test initConfig, setConfig, getConfig."""

    def test_set_and_get_config(self, clean_proj):
        clean_proj.setConfig(myKey="myValue", anotherKey=123)
        cfg = clean_proj.getConfig()
        assert cfg["myKey"] == "myValue"
        assert cfg["anotherKey"] == 123

    def test_set_config_keep_old_values_true(self, clean_proj):
        clean_proj.setConfig(keep_old_values=True, oldKey="old")
        clean_proj.setConfig(keep_old_values=True, newKey="new")
        cfg = clean_proj.getConfig()
        assert cfg.get("oldKey") == "old"
        assert cfg.get("newKey") == "new"

    def test_set_config_keep_old_values_false(self, clean_proj):
        clean_proj.setConfig(keep_old_values=True, first="one", second="two")
        clean_proj.setConfig(keep_old_values=False, only="this")
        cfg = clean_proj.getConfig()
        assert cfg.get("only") == "this"
        # Old keys should be gone
        assert "first" not in cfg
        assert "second" not in cfg

    def test_init_config_does_not_overwrite(self, clean_proj):
        clean_proj.setConfig(preserve="original")
        clean_proj.initConfig(preserve="overwritten", extra="added")
        cfg = clean_proj.getConfig()
        assert cfg["preserve"] == "original"
        assert cfg["extra"] == "added"


# ===========================================================================
# 9. Project properties
# ===========================================================================

@requires_mongo
class TestProjectProperties:
    """Test Project properties."""

    def test_project_name(self, proj):
        assert proj.projectName == PROJECT_NAME

    def test_measurements_property(self, proj):
        assert proj.measurements is not None

    def test_simulations_property(self, proj):
        assert proj.simulations is not None

    def test_cache_property(self, proj):
        assert proj.cache is not None

    def test_all_property(self, proj):
        assert proj.all is not None

    def test_files_directory(self, proj):
        fd = proj.filesDirectory
        # Should return a string or Path
        assert fd is not None


# ===========================================================================
# 10. saveData family — round-trip with real files
# ===========================================================================

@requires_mongo
class TestSaveData:
    """Test saveMeasurementData, saveCacheData, saveSimulationData with real data."""

    def test_save_measurement_dataframe(self, clean_proj, tmp_dir):
        clean_proj._filesDirectory = tmp_dir
        df = pd.DataFrame({"x": [1, 2, 3], "y": [4.0, 5.0, 6.0]})
        doc = clean_proj.saveMeasurementData(
            name="test_meas_df",
            data=df,
            desc={"purpose": "test"},
            dataFormat=datatypes.PARQUET,
        )
        assert doc is not None
        # Retrieve and load
        docs = list(clean_proj.getMeasurementsDocuments(purpose="test"))
        assert len(docs) >= 1
        loaded = docs[0].getData()
        assert isinstance(loaded, pd.DataFrame)
        assert len(loaded) == 3
        assert list(loaded.columns) == ["x", "y"]

    def test_save_cache_json_dict(self, clean_proj, tmp_dir):
        clean_proj._filesDirectory = tmp_dir
        data = {"result": 42, "labels": ["a", "b"]}
        doc = clean_proj.saveCacheData(
            name="test_cache_json",
            data=data,
            desc={"step": "final"},
            dataFormat=datatypes.JSON_DICT,
        )
        assert doc is not None
        docs = list(clean_proj.getCacheDocuments(step="final"))
        assert len(docs) >= 1
        loaded = docs[0].getData()
        assert loaded["result"] == 42
        assert loaded["labels"] == ["a", "b"]

    def test_save_simulation_numpy(self, clean_proj, tmp_dir):
        clean_proj._filesDirectory = tmp_dir
        arr = np.array([1.0, 2.0, 3.0, 4.0])
        doc = clean_proj.saveSimulationData(
            name="test_sim_np",
            data=arr,
            desc={"grid": "10m"},
            dataFormat=datatypes.NUMPY_ARRAY,
        )
        assert doc is not None
        docs = list(clean_proj.getSimulationsDocuments(grid="10m"))
        assert len(docs) >= 1
        loaded = docs[0].getData()
        assert isinstance(loaded, np.ndarray)
        assert np.allclose(loaded, arr)

    def test_save_string_data(self, clean_proj, tmp_dir):
        clean_proj._filesDirectory = tmp_dir
        doc = clean_proj.saveMeasurementData(
            name="test_string",
            data="hello world",
            desc={"format": "text"},
            dataFormat=datatypes.STRING,
        )
        assert doc is not None
        docs = list(clean_proj.getMeasurementsDocuments(format="text"))
        assert len(docs) >= 1
        loaded = docs[0].getData()
        assert "hello" in str(loaded)

    def test_save_csv_pandas(self, clean_proj, tmp_dir):
        clean_proj._filesDirectory = tmp_dir
        df = pd.DataFrame({"name": ["Alice", "Bob"], "score": [95, 87]})
        doc = clean_proj.saveMeasurementData(
            name="test_csv",
            data=df,
            desc={"source": "exam"},
            dataFormat=datatypes.CSV_PANDAS,
        )
        assert doc is not None
        docs = list(clean_proj.getMeasurementsDocuments(source="exam"))
        assert len(docs) >= 1
        loaded = docs[0].getData()
        assert isinstance(loaded, pd.DataFrame)
        assert len(loaded) == 2

    def test_save_pickle(self, clean_proj, tmp_dir):
        clean_proj._filesDirectory = tmp_dir
        data = {"complex": [1, 2, 3], "nested": {"a": True}}
        doc = clean_proj.saveCacheData(
            name="test_pickle",
            data=data,
            desc={"kind": "pickled"},
            dataFormat=datatypes.PICKLE,
        )
        assert doc is not None
        docs = list(clean_proj.getCacheDocuments(kind="pickled"))
        assert len(docs) >= 1
        loaded = docs[0].getData()
        assert loaded["complex"] == [1, 2, 3]
        assert loaded["nested"]["a"] is True


# ===========================================================================
# 11. DataHandler round-trips (unit tests, no DB needed)
# ===========================================================================

class TestDataHandlers:
    """Test DataHandler saveData/getData round-trips for key formats."""

    def test_handler_string(self, tmp_dir):
        handler = getHandler(datatypes.STRING)
        data = "test string content"
        filepath = os.path.join(tmp_dir, "handler_string.txt")
        handler.saveData(data, filepath)
        loaded = handler.getData(filepath)
        # String handler returns the resource path as-is
        assert loaded == filepath or loaded == data

    def test_handler_json_dict(self, tmp_dir):
        handler = getHandler(datatypes.JSON_DICT)
        data = {"key": "value", "number": 42, "list": [1, 2, 3]}
        filepath = os.path.join(tmp_dir, "handler_json.json")
        handler.saveData(data, filepath)
        loaded = handler.getData(filepath)
        assert loaded == data

    def test_handler_parquet(self, tmp_dir):
        handler = getHandler(datatypes.PARQUET)
        df = pd.DataFrame({"a": [1, 2, 3], "b": [4.0, 5.0, 6.0]})
        filepath = os.path.join(tmp_dir, "handler_parquet.parquet")
        handler.saveData(df, filepath)
        # Default returns Dask DataFrame; compute to pandas
        loaded = handler.getData(filepath)
        if hasattr(loaded, "compute"):
            loaded = loaded.compute()
        assert isinstance(loaded, pd.DataFrame)
        assert len(loaded) == 3
        pd.testing.assert_frame_equal(loaded, df)

    def test_handler_csv_pandas(self, tmp_dir):
        handler = getHandler(datatypes.CSV_PANDAS)
        df = pd.DataFrame({"x": [10, 20], "y": ["a", "b"]})
        filepath = os.path.join(tmp_dir, "handler_csv.csv")
        handler.saveData(df, filepath)
        loaded = handler.getData(filepath)
        assert isinstance(loaded, pd.DataFrame)
        assert len(loaded) == 2

    def test_handler_numpy(self, tmp_dir):
        handler = getHandler(datatypes.NUMPY_ARRAY)
        arr = np.array([[1, 2], [3, 4], [5, 6]])
        filepath = os.path.join(tmp_dir, "handler_np.npy")
        handler.saveData(arr, filepath)
        loaded = handler.getData(filepath)
        assert isinstance(loaded, np.ndarray)
        assert np.array_equal(loaded, arr)

    def test_handler_pickle(self, tmp_dir):
        handler = getHandler(datatypes.PICKLE)
        data = {"nested": {"deep": [1, 2, 3]}, "flag": True}
        filepath = os.path.join(tmp_dir, "handler_pickle.pkl")
        handler.saveData(data, filepath)
        loaded = handler.getData(filepath)
        assert loaded == data


# ===========================================================================
# 12. datatypes class methods
# ===========================================================================

class TestDatatypes:
    """Test datatypes introspection methods."""

    def test_get_handler_known_types(self):
        for fmt in [datatypes.STRING, datatypes.PARQUET, datatypes.CSV_PANDAS,
                     datatypes.JSON_DICT, datatypes.PICKLE, datatypes.NUMPY_ARRAY]:
            handler = datatypes.getHandler(fmt)
            assert handler is not None
            assert hasattr(handler, "saveData")
            assert hasattr(handler, "getData")

    def test_get_handler_unknown_raises(self):
        with pytest.raises(Exception):
            datatypes.getHandler("nonexistent_format_xyz")

    def test_guess_handler_dataframe(self):
        df = pd.DataFrame({"a": [1]})
        handler = datatypes.guessHandler(df)
        assert handler is not None

    def test_guess_handler_numpy(self):
        arr = np.array([1, 2, 3])
        handler = datatypes.guessHandler(arr)
        assert handler is not None

    def test_guess_handler_dict(self):
        handler = datatypes.guessHandler({"key": "value"})
        assert handler is not None

    def test_get_data_format_name_dataframe(self):
        df = pd.DataFrame({"a": [1]})
        name = datatypes.getDataFormatName(df)
        assert isinstance(name, str)
        assert len(name) > 0

    def test_get_data_format_extension(self):
        df = pd.DataFrame({"a": [1]})
        ext = datatypes.getDataFormatExtension(df)
        assert isinstance(ext, str)


# ===========================================================================
# 13. Collection direct access
# ===========================================================================

@requires_mongo
class TestCollectionDirect:
    """Test AbstractCollection and subclasses directly."""

    def test_measurements_collection_add_and_get(self):
        coll = Measurements_Collection()
        doc = coll.addDocument(
            projectName="pytest_coll_direct",
            resource="coll_test",
            dataFormat="string",
            type="CollTest",
            desc={"direct": True},
        )
        assert doc is not None

        docs = list(coll.getDocuments(projectName="pytest_coll_direct", type="CollTest"))
        assert len(docs) >= 1

        # Cleanup
        coll.deleteDocuments(projectName="pytest_coll_direct", type="CollTest")

    def test_simulations_collection_add_and_get(self):
        coll = Simulations_Collection()
        doc = coll.addDocument(
            projectName="pytest_coll_direct",
            resource="sim_coll",
            dataFormat="string",
            type="SimCollTest",
            desc={},
        )
        assert doc is not None

        docs = list(coll.getDocuments(projectName="pytest_coll_direct", type="SimCollTest"))
        assert len(docs) >= 1

        # Cleanup
        coll.deleteDocuments(projectName="pytest_coll_direct", type="SimCollTest")

    def test_cache_collection_add_and_get(self):
        coll = Cache_Collection()
        doc = coll.addDocument(
            projectName="pytest_coll_direct",
            resource="cache_coll",
            dataFormat="string",
            type="CacheCollTest",
            desc={},
        )
        assert doc is not None

        docs = list(coll.getDocuments(projectName="pytest_coll_direct", type="CacheCollTest"))
        assert len(docs) >= 1

        # Cleanup
        coll.deleteDocuments(projectName="pytest_coll_direct", type="CacheCollTest")

    def test_get_documents_as_dict(self):
        coll = Measurements_Collection()
        coll.addDocument(
            projectName="pytest_coll_direct",
            resource="dict_coll",
            dataFormat="string",
            type="DictCollTest",
            desc={"val": 1},
        )
        result = coll.getDocumentsAsDict(projectName="pytest_coll_direct", type="DictCollTest")
        assert "documents" in result
        assert len(result["documents"]) >= 1

        # Cleanup
        coll.deleteDocuments(projectName="pytest_coll_direct", type="DictCollTest")

    def test_delete_document_by_id(self):
        coll = Measurements_Collection()
        doc = coll.addDocument(
            projectName="pytest_coll_direct",
            resource="del_by_id",
            dataFormat="string",
            type="DelByID",
            desc={},
        )
        doc_id = str(doc.id)
        deleted = coll.deleteDocumentByID(doc_id)
        assert deleted is not None

        # Verify deletion
        docs = list(coll.getDocuments(projectName="pytest_coll_direct", type="DelByID"))
        assert len(docs) == 0

    def test_get_document_by_id(self):
        coll = Measurements_Collection()
        doc = coll.addDocument(
            projectName="pytest_coll_direct",
            resource="get_by_id",
            dataFormat="string",
            type="GetByID",
            desc={"marker": "found"},
        )
        doc_id = str(doc.id)
        retrieved = coll.getDocumentByID(doc_id)
        assert retrieved is not None
        assert retrieved.desc["marker"] == "found"

        # Cleanup
        coll.deleteDocuments(projectName="pytest_coll_direct", type="GetByID")

    def test_get_project_list(self):
        coll = Measurements_Collection()
        coll.addDocument(
            projectName="pytest_project_list_test",
            resource="pl_test",
            dataFormat="string",
            type="PLTest",
            desc={},
        )
        projects = coll.getProjectList()
        assert "pytest_project_list_test" in projects

        # Cleanup
        coll.deleteDocuments(projectName="pytest_project_list_test", type="PLTest")


# ===========================================================================
# 14. MetadataFrame
# ===========================================================================

@requires_mongo
class TestMetadataFrame:
    """Test MetadataFrame.asDict and MetadataFrame.getData."""

    def test_as_dict(self, clean_proj):
        doc = clean_proj.addMeasurementsDocument(
            resource="as_dict_test",
            dataFormat="string",
            type="AsDictTest",
            desc={"key": "value", "num": 7},
        )
        d = doc.asDict()
        assert isinstance(d, dict)
        assert d.get("type") == "AsDictTest"
        assert d.get("desc", {}).get("key") == "value"

    def test_as_dict_with_id(self, clean_proj):
        doc = clean_proj.addMeasurementsDocument(
            resource="as_dict_id",
            dataFormat="string",
            type="AsDictID",
            desc={},
        )
        d = doc.asDict(with_id=True)
        assert "_id" in d

    def test_get_data_string(self, clean_proj):
        doc = clean_proj.addMeasurementsDocument(
            resource="hello_resource",
            dataFormat="string",
            type="GetDataStr",
            desc={},
        )
        data = doc.getData()
        # String handler returns the resource string
        assert data is not None


# ===========================================================================
# 15. nonDBMetadataFrame
# ===========================================================================

class TestNonDBMetadataFrame:
    """Test nonDBMetadataFrame wrapping."""

    def test_init_and_get_data(self):
        df = pd.DataFrame({"col": [1, 2, 3]})
        frame = nonDBMetadataFrame(data=df, projectName="test", type="test_type")
        result = frame.getData()
        assert isinstance(result, pd.DataFrame)
        pd.testing.assert_frame_equal(result, df)

    def test_getitem(self):
        frame = nonDBMetadataFrame(
            data="test_data",
            projectName="proj",
            type="mytype",
            resource="/path",
            dataFormat="string",
        )
        assert frame["type"] == "mytype"
        assert frame["resource"] == "/path"
        assert frame["dataFormat"] == "string"

    def test_get_data_with_dict(self):
        data = {"key": "value"}
        frame = nonDBMetadataFrame(data=data)
        assert frame.getData() == data

    def test_get_data_with_numpy(self):
        arr = np.array([1, 2, 3])
        frame = nonDBMetadataFrame(data=arr)
        result = frame.getData()
        assert np.array_equal(result, arr)


# ===========================================================================
# 16. autocache — @cacheFunction decorator
# ===========================================================================

@requires_mongo
class TestAutoCache:
    """Test cacheFunction decorator with real MongoDB caching."""

    def test_cache_function_basic(self, tmp_dir):
        from hera.datalayer.autocache import cacheFunction, clearFunctionCache

        call_count = {"n": 0}

        @cacheFunction(returnFormat=datatypes.JSON_DICT, projectName="pytest_autocache")
        def expensive_computation(x, y):
            call_count["n"] += 1
            return {"result": x + y}

        # First call — should execute and cache
        result1 = expensive_computation(3, 4)
        assert result1["result"] == 7
        first_count = call_count["n"]

        # Second call with same args — should return cached
        result2 = expensive_computation(3, 4)
        assert result2["result"] == 7
        # Function should not have been called again
        assert call_count["n"] == first_count

        # Different args — should execute again
        result3 = expensive_computation(10, 20)
        assert result3["result"] == 30
        assert call_count["n"] == first_count + 1

        # Cleanup
        clearFunctionCache("expensive_computation", projectName="pytest_autocache")

    def test_clear_function_cache(self):
        from hera.datalayer.autocache import cacheFunction, clearFunctionCache

        call_count = {"n": 0}

        @cacheFunction(returnFormat=datatypes.JSON_DICT, projectName="pytest_autocache_clear")
        def cached_add(a, b):
            call_count["n"] += 1
            return {"sum": a + b}

        # Call and cache
        cached_add(1, 2)
        count_after_first = call_count["n"]

        # Clear cache
        clearFunctionCache("cached_add", projectName="pytest_autocache_clear")

        # Call again — should re-execute
        cached_add(1, 2)
        assert call_count["n"] == count_after_first + 1

        # Cleanup
        clearFunctionCache("cached_add", projectName="pytest_autocache_clear")

    def test_clear_all_functions_cache(self):
        from hera.datalayer.autocache import cacheFunction, clearAllFunctionsCache

        @cacheFunction(returnFormat=datatypes.JSON_DICT, projectName="pytest_autocache_all")
        def func_a():
            return {"a": 1}

        @cacheFunction(returnFormat=datatypes.JSON_DICT, projectName="pytest_autocache_all")
        def func_b():
            return {"b": 2}

        func_a()
        func_b()

        # Clear all
        clearAllFunctionsCache(projectName="pytest_autocache_all")

        # Verify cache cleared by checking documents
        p = Project(projectName="pytest_autocache_all")
        docs = list(p.getCacheDocuments(type="functionCacheData"))
        assert len(docs) == 0


# ===========================================================================
# 17. Export / Import round-trip
# ===========================================================================

@requires_mongo
class TestExportImport:
    """Test Project export and load round-trip."""

    def test_export_and_load(self, tmp_dir):
        # Create source project with documents
        src_name = "pytest_export_src"
        src = Project(projectName=src_name)
        src.addMeasurementsDocument(
            resource="export_test_1", dataFormat="string", type="ExportType",
            desc={"station": "A", "value": 100},
        )
        src.addMeasurementsDocument(
            resource="export_test_2", dataFormat="string", type="ExportType",
            desc={"station": "B", "value": 200},
        )
        src.addSimulationsDocument(
            resource="export_sim", dataFormat="string", type="ExportSim",
            desc={"run": 1},
        )

        # Export
        export_path = os.path.join(tmp_dir, "export_test")
        src.export(export_path, show_progressbar=False)

        # Create target project and load
        dst_name = "pytest_export_dst"
        dst = Project(projectName=dst_name)
        Project.load(dst, export_path, is_hard_import=False, show_progressbar=False)

        # Verify
        meas_docs = list(dst.getMeasurementsDocuments(type="ExportType"))
        assert len(meas_docs) == 2

        sim_docs = list(dst.getSimulationsDocuments(type="ExportSim"))
        assert len(sim_docs) == 1

        # Cleanup both projects
        for p in (src, dst):
            for getter in (p.getMeasurementsDocuments, p.getSimulationsDocuments, p.getCacheDocuments):
                try:
                    for doc in getter():
                        doc.delete()
                except Exception:
                    pass


# ===========================================================================
# 18. Module-level functions
# ===========================================================================

@requires_mongo
class TestModuleFunctionsDB:
    """Test module-level functions that require MongoDB."""

    def test_get_project_list(self, proj):
        # proj fixture ensures at least one document exists
        proj.addMeasurementsDocument(
            resource="pl_func", dataFormat="string", type="PLFunc", desc={},
        )
        projects = getProjectList()
        assert isinstance(projects, list)
        assert PROJECT_NAME in projects


class TestModuleFunctionsLocal:
    """Test module-level functions that do NOT require MongoDB."""

    def test_create_project_directory(self, tmp_dir):
        out = os.path.join(tmp_dir, "new_project_dir")
        createProjectDirectory(out, projectName="TestCreatedProject")
        assert os.path.isdir(out)
        # Check config file created
        config_path = os.path.join(out, "caseConfiguration.json")
        assert os.path.isfile(config_path)
        with open(config_path) as f:
            cfg = json.load(f)
        assert cfg.get("projectName") == "TestCreatedProject"
