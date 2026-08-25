"""experimentAnalysis and the data-engine factory.

Three defects surfaced here:

* B75: ``experimentAnalysis.getTurbulenceStatistics`` references
  ``kaijoData``/``kaijoHeight`` -- names that do not exist anywhere in the
  method, which only takes ``sonicData``/``height``. Every call raises
  ``NameError`` before it can even reach the toolkit it is trying to call.
* B76: both ``pandasDataEngineDB.DBconfiguration`` and
  ``daskDataEngineDB.DBconfiguration`` are missing a ``return`` statement
  (``self._dbConfiguration`` as a bare expression), so the property always
  returns ``None`` regardless of what was actually configured.
* B77: ``pandasDataEngineDB.dbConnect`` catches
  ``except OperationFailure as e:``, but ``OperationFailure`` is never
  imported in this module. When a real connection attempt fails with any
  other exception, Python's own attempt to match it against
  ``OperationFailure`` raises ``NameError`` -- which the outer
  ``except Exception as e: return e`` then silently swallows and returns
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
    dataEngineFactory,
    daskDataEngineDB,
    pandasDataEngineDB,
)


@pytest.mark.unit
class TestGetTurbulenceStatisticsIsBroken:
    @pytest.mark.xfail(
        strict=True,
        reason="B75: the method reads kaijoData/kaijoHeight, names that "
               "exist nowhere in its own body or parameters (only "
               "sonicData/height do) -- every call raises NameError. "
               "See the consolidated findings issue.",
    )
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
        analysis.getTurbulenceStatistics(sonicData=object(), samplingWindow="1min", height=2)

    def test_it_currently_raises_a_nameerror(self):
        """Characterisation of B75."""
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
        with pytest.raises(NameError, match="kaijo"):
            analysis.getTurbulenceStatistics(sonicData=object(), samplingWindow="1min", height=2)


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

    @pytest.mark.xfail(
        strict=True,
        reason="B76: DBconfiguration evaluates self._dbConfiguration as a "
               "bare expression with no return, so it always yields None. "
               "See the consolidated findings issue.",
    )
    def test_db_configuration_should_return_the_stored_dict(self):
        engine = daskDataEngineDB("proj", {"DB": {"login": {
            "username": "u", "password": "p", "ip": "1.2.3.4", "port": "27017",
        }}})
        assert engine.DBconfiguration is not None

    def test_db_configuration_currently_returns_none(self):
        """Characterisation of B76."""
        engine = daskDataEngineDB("proj", {"DB": {"login": {
            "username": "u", "password": "p", "ip": "1.2.3.4", "port": "27017",
        }}})
        assert engine.DBconfiguration is None


@pytest.mark.unit
class TestPandasDataEngineDBConstructionValidation:
    def test_construction_requires_a_db_key(self):
        """The DB-key check runs before dbConnect(), so this never touches the network."""
        with pytest.raises(ValueError, match="DB"):
            pandasDataEngineDB("proj", {})

    @pytest.mark.xfail(
        strict=True,
        reason="B76: same missing return as daskDataEngineDB.DBconfiguration.",
    )
    def test_db_configuration_should_return_the_stored_dict(self):
        engine = pandasDataEngineDB.__new__(pandasDataEngineDB)
        engine._dbConfiguration = {"login": {}}
        assert engine.DBconfiguration is not None

    def test_db_configuration_currently_returns_none(self):
        """Characterisation of B76, on the pandas engine."""
        engine = pandasDataEngineDB.__new__(pandasDataEngineDB)
        engine._dbConfiguration = {"login": {}}
        assert engine.DBconfiguration is None


@pytest.mark.unit
class TestPandasDataEngineDBConnect:
    """dbConnect() opens a real socket via pymongo.MongoClient; every test
    here monkeypatches it so nothing ever touches the network."""

    @pytest.mark.xfail(
        strict=True,
        reason="B77: OperationFailure is never imported, so matching any "
               "raised exception against it raises NameError, which the "
               "outer except Exception silently swallows and returns "
               "instead of the real connection error. "
               "See the consolidated findings issue.",
    )
    def test_dbconnect_should_surface_the_original_exception(self, monkeypatch):
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

    def test_dbconnect_currently_masks_it_with_a_nameerror(self, monkeypatch):
        """Characterisation of B77."""
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
        assert isinstance(result, NameError)
        assert "OperationFailure" in str(result)
