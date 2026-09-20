"""abstractcalculator.py: AbstractCalculator's properties, save-properties
setter, and the DB-touching compute modes.

B105: every DB-touching code path (``_compute_from_db_and_save``,
``_compute_from_db_and_not_save``, ``_save_to_db`` -- and so
``_compute_not_from_db_and_save``, which calls it) does
``from hera.datalayer import collection as datalayer`` then
``datalayer.Cache.getDocuments(...)``/``datalayer.Cache.addDocument(...)``.
``hera.datalayer.collection`` only defines the ``Cache_Collection`` class,
not a ``Cache`` singleton instance -- that instance lives at
``hera.datalayer.Cache`` (package level), one import path up. Every mode
except the default ('not_from_db_and_not_save') is therefore permanently
dead code: it always raises ``AttributeError`` before doing anything useful.
"""
import pandas
import pytest

from hera.measurements.meteorology.highfreqdata.analysis.abstractcalculator import (
    AbstractCalculator,
)


@pytest.fixture(autouse=True)
def _reset_shared_save_properties():
    """B106: _saveProperties is a class-level mutable dict that __init__
    never copies into the instance, so set_saveProperties on any instance
    mutates it for every instance (existing or future) that hasn't
    shadowed it. Reset the shared dict around each test so this file's
    own tests stay isolated from each other and from any other test
    module that touches this class."""
    AbstractCalculator._saveProperties = {"dataFormat": None}
    yield
    AbstractCalculator._saveProperties = {"dataFormat": None}


@pytest.fixture()
def calc():
    raw = pandas.DataFrame({"u": [1.0, 2.0, 3.0]})
    metadata = dict(projectName="P", start="2020-01-01", end="2020-01-02", samplingWindow="30min")
    return AbstractCalculator(raw, metadata)


@pytest.fixture()
def calc_with_params(calc):
    calc._TemporaryData = pandas.DataFrame({"x": [1.0, 2.0, 3.0]})
    calc._CalculatedParams = [("x", {})]
    calc._AllCalculatedParams = []
    return calc


@pytest.mark.unit
class TestProperties:
    def test_rawdata_returns_the_original_frame(self, calc):
        assert calc.RawData is calc._RawData

    def test_temporarydata_starts_empty(self, calc):
        assert calc.TemporaryData.empty

    def test_samplingwindow_reads_from_metadata(self, calc):
        assert calc.SamplingWindow == "30min"

    def test_karman_is_the_default_von_karman_constant(self, calc):
        assert calc.Karman == pytest.approx(0.4)


@pytest.mark.unit
class TestSetSaveProperties:
    def test_it_stores_the_dataformat_and_extra_kwargs(self, calc):
        calc.set_saveProperties(dataFormat="JSON_pandas", path="/tmp/x")
        assert calc._saveProperties["dataFormat"] == "JSON_pandas"
        assert calc._saveProperties["path"] == "/tmp/x"


@pytest.mark.unit
class TestSaveSharedAcrossInstancesIsBroken:
    """B106: `_saveProperties = {'dataFormat': None}` is declared on the
    class body, and __init__ never assigns a fresh instance-level copy --
    so every instance shares (and mutates) the exact same dict until one
    happens to shadow it. set_saveProperties() on one calculator leaks
    into every other one, including calculators constructed afterwards."""

    @pytest.mark.xfail(
        strict=True,
        reason="B106: _saveProperties is a class-level mutable dict, never "
               "copied per-instance in __init__, so set_saveProperties on "
               "one AbstractCalculator mutates the save properties seen by "
               "every other instance. See the consolidated findings issue.",
    )
    def test_a_second_instance_should_be_unaffected_by_the_first(self, calc):
        calc.set_saveProperties(dataFormat="JSON_pandas")
        other = AbstractCalculator(calc.RawData, calc.metaData)
        assert other._saveProperties["dataFormat"] is None

    def test_a_second_instance_currently_leaks_the_first_ones_save_format(self, calc):
        """Characterisation of B106."""
        calc.set_saveProperties(dataFormat="JSON_pandas")
        other = AbstractCalculator(calc.RawData, calc.metaData)
        assert other._saveProperties["dataFormat"] == "JSON_pandas"
        assert other._saveProperties is calc._saveProperties


@pytest.mark.unit
class TestSaveToDbIsUnreachable:
    def test_with_no_save_properties_set_it_raises_before_touching_the_db(self, calc_with_params):
        """Documented, pre-existing guard: unrelated to B105."""
        with pytest.raises(AttributeError, match="save properties"):
            calc_with_params._save_to_db(["x"])

    def test_compute_not_from_db_and_save_hits_b105(self, calc_with_params):
        calc_with_params.set_saveProperties(dataFormat="JSON_pandas")
        with pytest.raises(AttributeError, match="Cache"):
            calc_with_params._compute_not_from_db_and_save()

    def test_compute_from_db_and_save_hits_b105(self, calc_with_params):
        calc_with_params.set_saveProperties(dataFormat="JSON_pandas")
        with pytest.raises(AttributeError, match="Cache"):
            calc_with_params._compute_from_db_and_save()

    def test_compute_from_db_and_not_save_hits_b105(self, calc_with_params):
        with pytest.raises(AttributeError, match="Cache"):
            calc_with_params._compute_from_db_and_not_save()

    def test_compute_dispatches_to_the_broken_mode_and_fails_the_same_way(self, calc_with_params):
        calc_with_params._AllCalculatedParams = []
        calc_with_params.set_saveProperties(dataFormat="JSON_pandas")
        with pytest.raises(AttributeError, match="Cache"):
            calc_with_params.compute(mode="not_from_db_and_save")

    def test_only_the_default_mode_actually_works(self, calc_with_params):
        result = calc_with_params.compute(mode="not_from_db_and_not_save")
        assert list(result.columns) == ["x"]
