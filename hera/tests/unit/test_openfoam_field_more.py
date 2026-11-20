"""OFField.__getitem__ / OFField.__iter__ -- the two members left uncovered.

``test_openfoam_field.py`` covers the initialization, I/O and boundary
methods; these two dunders are the mapping protocol OFField puts over its
``data`` dict, keyed by processor name (``singleProcessor`` for a case that
is not decomposed, ``processor0``, ``processor1``, ... for one that is).

Both are one-liners that delegate straight to ``self.data``, so what is
asserted is the delegation and the key set -- which OFObjectHome/OFField
decide -- never the contents of a value.  PyFoam is the conftest's
``MagicMock`` stub (see ``_stubs.py``), so each ``WriteParameterFile(...)``
returns the *same* shared mock object; identity of a value is therefore not
a usable assertion, but identity against ``field.data[key]`` is, because
that is exactly what the dunder must return.

No bugs found in either.
"""
import pytest

from hera.simulations.openFoam.preprocessOFObjects.OFObjectHome import OFObjectHome


@pytest.fixture()
def home():
    return OFObjectHome()


@pytest.fixture()
def singleProcessorField(home):
    """getEmptyField's default: a field for a case that is not decomposed."""
    return home.getEmptyField("p", flowType="incompressible")


@pytest.fixture()
def parallelField(home):
    return home.getEmptyField("U", flowType="incompressible", noOfProc=3)


@pytest.mark.unit
class TestGetItem:
    def test_a_single_processor_field_is_keyed_by_singleprocessor(self, singleProcessorField):
        assert singleProcessorField["singleProcessor"] is \
            singleProcessorField.data["singleProcessor"]

    def test_every_processor_of_a_parallel_field_is_addressable(self, parallelField):
        for index in range(3):
            key = f"processor{index}"
            assert parallelField[key] is parallelField.data[key]

    def test_an_unknown_processor_name_raises_keyerror(self, singleProcessorField):
        with pytest.raises(KeyError, match="processor0"):
            singleProcessorField["processor0"]

    def test_the_single_processor_key_is_absent_from_a_parallel_field(self, parallelField):
        with pytest.raises(KeyError, match="singleProcessor"):
            parallelField["singleProcessor"]

    def test_it_reads_through_to_the_dict_rather_than_copying_it(self, singleProcessorField):
        """Replacing an entry must be visible through the subscript."""
        sentinel = object()
        singleProcessorField.data["singleProcessor"] = sentinel
        assert singleProcessorField["singleProcessor"] is sentinel

    def test_a_key_added_after_construction_is_reachable(self, singleProcessorField):
        sentinel = object()
        singleProcessorField.data["processor7"] = sentinel
        assert singleProcessorField["processor7"] is sentinel


@pytest.mark.unit
class TestIter:
    def test_iterating_a_single_processor_field_yields_its_one_key(self, singleProcessorField):
        assert list(singleProcessorField) == ["singleProcessor"]

    def test_iterating_a_parallel_field_yields_one_key_per_processor(self, parallelField):
        assert sorted(parallelField) == ["processor0", "processor1", "processor2"]

    def test_it_yields_the_keys_in_the_order_the_dict_holds_them(self, parallelField):
        assert list(parallelField) == list(parallelField.data)

    def test_it_returns_an_iterator_so_the_field_can_be_iterated_twice(self, parallelField):
        first = list(parallelField)
        second = list(parallelField)
        assert first == second == list(parallelField.processors)

    def test_the_keys_it_yields_are_the_keys_the_subscript_accepts(self, parallelField):
        assert [parallelField[key] for key in parallelField] == list(parallelField.data.values())

    def test_the_membership_test_that_iteration_provides_agrees_with_the_dict(
            self, parallelField, singleProcessorField):
        assert "processor1" in list(parallelField)
        assert "processor1" not in list(singleProcessorField)
