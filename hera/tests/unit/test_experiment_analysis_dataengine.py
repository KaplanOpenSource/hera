"""experimentAnalysis and the data-engine factory.

Three defects surfaced here:

* B75 (fixed): ``experimentAnalysis.getTurbulenceStatistics`` referenced
  ``kaijoData``/``kaijoHeight`` -- names that do not exist anywhere in the
  method, which only takes ``sonicData``/``height`` -- so every call raised
  ``NameError`` before it could even reach the toolkit it is trying to call.
* B76 (fixed): both ``pandasDataEngineDB.DBconfiguration`` and
  ``daskDataEngineDB.DBconfiguration`` were missing a ``return`` statement
  (``self._dbConfiguration`` as a bare expression), so the property always
  returned ``None`` regardless of what was actually configured.
* B77 (fixed): ``pandasDataEngineDB.dbConnect`` catches
  ``except OperationFailure as e:``, but ``OperationFailure`` was never
  imported in this module. When a real connection attempt failed with any
  other exception, Python's own attempt to match it against
  ``OperationFailure`` raised ``NameError`` -- which the outer
  ``except Exception as e: return e`` then silently swallowed and returned
  in place of the real connection error.

Note: constructing ``pandasDataEngineDB`` for real calls ``dbConnect()``,
which opens an actual (blocked-by-the-unit-harness-eventually, but not
before pymongo's server-selection timeout) socket -- one blind attempt
cost this batch a genuine 30-second wait. Every test below either
monkeypatches ``pymongo.MongoClient`` or bypasses ``__init__`` with
``__new__`` to stay hermetic and fast.
"""
import pytest

from hera.measurements.experiment.dataEngine import (
    DASKDB,
    daskDataEngineDB,
    dataEngineFactory,
    pandasDataEngineDB,
)


@pytest.mark.unit
class TestGetTurbulenceStatistics:
    def test_it_should_return_a_turbulence_statistics_object(self):
        from hera.measurements.experiment.analysis import experimentAnalysis

        class _FakeHighFreqToolkit:
            class analysis:
                @staticmethod
                def singlePointTurbulenceStatistics(**kwargs):
                    return object()

        class _FakeExtension:
            sonicHighFreqToolkit = _FakeHighFreqToolkit()

        class _FakeDataLayer:
            toolkitExtension = _FakeExtension()

        analysis = experimentAnalysis(_FakeDataLayer())
        result = analysis.getTurbulenceStatistics(sonicData=object(), samplingWindow="1min", height=2)
        assert result is not None


    def test_it_forwards_its_own_arguments_to_the_toolkit(self):
        """B75: sonicData/height reach the toolkit, not undefined names."""
        from hera.measurements.experiment.analysis import experimentAnalysis

        seen = {}

        class _FakeHighFreqToolkit:
            class analysis:
                @staticmethod
                def singlePointTurbulenceStatistics(**kwargs):
                    seen.update(kwargs)
                    return object()

        class _FakeExtension:
            sonicHighFreqToolkit = _FakeHighFreqToolkit()

        class _FakeDataLayer:
            toolkitExtension = _FakeExtension()

        sonicData = object()
        analysis = experimentAnalysis(_FakeDataLayer())
        analysis.getTurbulenceStatistics(sonicData=sonicData, samplingWindow="1min", height=2)
        assert seen["sonicData"] is sonicData
        assert seen["height"] == 2
        assert seen["samplingWindow"] == "1min"


@pytest.mark.unit
class TestDataEngineFactory:
    def test_daskdb_dispatches_to_dask_data_engine(self):
        engine = dataEngineFactory().getDataEngine(
            "proj", {"DB": {"login": {"username": "u", "password": "p", "ip": "1.2.3.4", "port": "27017"}}},
            experimentObj=None, dataType=DASKDB,
        )
        assert isinstance(engine, daskDataEngineDB)

    def test_an_unknown_data_type_raises(self):
        with pytest.raises(NotImplementedError):
            dataEngineFactory().getDataEngine("proj", {}, experimentObj=None, dataType="NoSuchEngine")


@pytest.mark.unit
class TestDaskDataEngineDB:
    def test_construction_requires_a_db_key(self):
        with pytest.raises(ValueError, match="DB"):
            daskDataEngineDB("proj", {})

    def test_construction_does_not_touch_the_network(self):
        """No dbConnect() call in __init__ -- this must be instant."""
        engine = daskDataEngineDB("proj", {"DB": {"login": {
            "username": "u", "password": "p", "ip": "1.2.3.4", "port": "27017",
        }}})
        assert engine._dbConfiguration["login"]["ip"] == "1.2.3.4"

    def test_connection_string_embeds_the_login_details(self):
        engine = daskDataEngineDB("proj", {"DB": {"login": {
            "username": "u", "password": "p", "ip": "1.2.3.4", "port": "27017",
        }}})
        assert engine.connectionString == "mongodb://u:p@1.2.3.4:27017/"

    def test_db_configuration_returns_the_stored_dict(self):
        engine = daskDataEngineDB("proj", {"DB": {"login": {
            "username": "u", "password": "p", "ip": "1.2.3.4", "port": "27017",
        }}})
        assert engine.DBconfiguration == engine._dbConfiguration
        assert engine.DBconfiguration["login"]["ip"] == "1.2.3.4"


@pytest.mark.unit
class TestPandasDataEngineDBConstructionValidation:
    def test_construction_requires_a_db_key(self):
        """The DB-key check runs before dbConnect(), so this never touches the network."""
        with pytest.raises(ValueError, match="DB"):
            pandasDataEngineDB("proj", {})

    def test_db_configuration_returns_the_stored_dict(self):
        engine = pandasDataEngineDB.__new__(pandasDataEngineDB)
        engine._dbConfiguration = {"login": {}}
        assert engine.DBconfiguration == {"login": {}}


@pytest.mark.unit
class TestPandasDataEngineDBConnect:
    """dbConnect() opens a real socket via pymongo.MongoClient; every test
    here monkeypatches it so nothing ever touches the network."""

    def test_dbconnect_surfaces_the_original_exception(self, monkeypatch):
        import hera.measurements.experiment.dataEngine as de_mod

        class _FakeClient:
            def server_info(self):
                raise RuntimeError("connection refused")

        monkeypatch.setattr(de_mod.pymongo, "MongoClient", lambda uri: _FakeClient())
        engine = pandasDataEngineDB.__new__(pandasDataEngineDB)
        engine._dbConfiguration = {"login": {
            "username": "u", "password": "p", "ip": "1.2.3.4", "port": "27017",
        }}
        result = engine.dbConnect()
        assert isinstance(result, RuntimeError)

    def test_dbconnect_returns_an_operation_failure_as_is(self, monkeypatch):
        """B77: the OperationFailure branch is now reachable."""
        import hera.measurements.experiment.dataEngine as de_mod

        class _FakeClient:
            def server_info(self):
                raise de_mod.OperationFailure("auth failed")

        monkeypatch.setattr(de_mod.pymongo, "MongoClient", lambda uri: _FakeClient())
        engine = pandasDataEngineDB.__new__(pandasDataEngineDB)
        engine._dbConfiguration = {"login": {
            "username": "u", "password": "p", "ip": "1.2.3.4", "port": "27017",
        }}
        result = engine.dbConnect()
        assert isinstance(result, de_mod.OperationFailure)
