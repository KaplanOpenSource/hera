"""OFField: the five public methods batch 11 didn't reach.

PyFoam's ``DictProxy``/``Field``/``WriteParameterFile`` are stubbed as a bare
``MagicMock`` module (see ``_stubs.py``), so they carry no real storage --
``mock["key"] = value`` records the call but a later ``mock["key"]`` does not
return it, and every call to e.g. ``WriteParameterFile(...)`` returns the
same shared mock object. That makes the field-initialization methods
(``addProcBoundary``, ``readBoundariesFromCase``, ``readFromCase``)
untestable for their actual written values -- the same ceiling batch 11 hit
for ``addBoundaryField`` -- so they get a "runs without raising" smoke test
plus, where the method has a real ``except FileNotFoundError`` branch, a
monkeypatch of the (real) PyFoam name to force it.

``setFieldFromDataFrame`` and ``getDataFrame`` are different: their inputs
and outputs are plain pandas, so they're tested against hand-built frames
shaped like ``ParsedParameterFileToDataFrame`` actually produces (see batch
13). That surfaced two real defects:

* B67: ``_processorToPyFoam`` always does
  ``processordataFrame.query("region=='boundaryField'").groupby("boundary")``,
  which requires a ``boundary`` column to exist even when there are no
  boundary rows at all. A field with no boundary patches -- exactly what
  ``getDataFrame()`` returns for one -- cannot be fed back through
  ``setFieldFromDataFrame`` without that column: it turns a same-shape
  round trip into a guaranteed ``KeyError``.
* B68: the parallel branch of ``setFieldFromDataFrame`` reads
  ``self.data[processorName]['boundaryField']`` to merge in the *existing*
  boundary before writing -- but a freshly-built field (the common case,
  via ``OFObjectHome.getEmptyField``) is initialized single-processor, so
  ``self.data`` has no ``processorN`` keys yet. Feeding it a parallel
  DataFrame straight away raises ``KeyError``; it only works if the field
  was already initialized parallel with ``initialize(noOfProc=...)``.

Note also: ``preprocessOFObjects/__init__.py`` does
``from .OFField import OFField``, which shadows the ``OFField`` *module*
attribute on the package with the *class* of the same name. Reaching the
real module to monkeypatch ``ParsedParameterFile``/``BoundaryDict`` needs
``sys.modules["hera.simulations.openFoam.preprocessOFObjects.OFField"]``,
not attribute access through the package.
"""
import sys

import pandas
import pytest

from hera.simulations.openFoam.preprocessOFObjects.OFObjectHome import OFObjectHome


class _FakeVal:
    def __init__(self, val):
        self.val = val


@pytest.fixture()
def home():
    return OFObjectHome()


@pytest.fixture()
def scalar_field(home):
    return home.getEmptyField("p", flowType="incompressible")


@pytest.fixture()
def of_field_module():
    import hera.simulations.openFoam.preprocessOFObjects.OFField  # noqa: F401

    return sys.modules["hera.simulations.openFoam.preprocessOFObjects.OFField"]


@pytest.mark.unit
class TestFieldInitializationSmoke:
    """Weak by necessity -- see the module docstring on the PyFoam stub."""

    def test_add_proc_boundary_does_not_raise(self, scalar_field):
        scalar_field.addProcBoundary()
        assert scalar_field.data is not None

    def test_read_boundaries_from_case_on_an_empty_directory_falls_back_to_single(
        self, scalar_field, tmp_path
    ):
        scalar_field.readBoundariesFromCase(str(tmp_path), readParallel=True)
        assert "singleProcessor" in scalar_field.data

    def test_read_from_case_on_an_empty_directory_falls_back_to_single(self, scalar_field, tmp_path):
        scalar_field.readFromCase(str(tmp_path), timeStep=0, readParallel=True)
        assert "singleProcessor" in scalar_field.data


@pytest.mark.unit
class TestFieldNotFoundIsWrappedAsValueError:
    def test_read_from_case_wraps_a_missing_file(self, scalar_field, tmp_path, monkeypatch, of_field_module):
        def boom(*a, **kw):
            raise FileNotFoundError("no such file")

        monkeypatch.setattr(of_field_module, "ParsedParameterFile", boom)
        with pytest.raises(ValueError, match="Field not found"):
            scalar_field.readFromCase(str(tmp_path), timeStep=0, readParallel=False)

    def test_read_boundaries_from_case_wraps_a_missing_file(
        self, scalar_field, tmp_path, monkeypatch, of_field_module
    ):
        def boom(*a, **kw):
            raise FileNotFoundError("no such file")

        monkeypatch.setattr(of_field_module, "BoundaryDict", boom)
        with pytest.raises(ValueError, match="Field not found"):
            scalar_field.readBoundariesFromCase(str(tmp_path), readParallel=False)


@pytest.mark.unit
class TestGetDataFrame:
    def test_a_single_processor_field_becomes_one_frame(self, scalar_field):
        scalar_field.data = {
            "singleProcessor": {
                "internalField": _FakeVal([[1.0], [2.0]]),
                "boundaryField": {},
            }
        }
        df = scalar_field.getDataFrame()
        assert list(df["p"]) == [1.0, 2.0]
        assert list(df["region"]) == ["internalField", "internalField"]

    def test_a_multi_processor_field_concatenates_all_processors(self, scalar_field):
        scalar_field.data = {
            "processor0": {"internalField": _FakeVal([[1.0]]), "boundaryField": {}},
            "processor1": {"internalField": _FakeVal([[2.0]]), "boundaryField": {}},
        }
        df = scalar_field.getDataFrame()
        assert sorted(df["p"]) == [1.0, 2.0]
        assert set(df["processor"]) == {0, 1}

    def test_the_processor_number_is_parsed_from_the_key_name(self, scalar_field):
        """'processor7' -> 7, sliced from character index 9 onward."""
        scalar_field.data = {"processor7": {"internalField": _FakeVal([[1.0]]), "boundaryField": {}}}
        df = scalar_field.getDataFrame()
        assert list(df["processor"]) == [7]


@pytest.mark.unit
class TestSetFieldFromDataFrame:
    def test_a_dataframe_without_a_processor_column_is_treated_as_single(self, scalar_field):
        df = pandas.DataFrame({
            "region": ["internalField"],
            "processorIndex": [0],
            "boundary": [None],
            "p": [5.0],
        })
        scalar_field.setFieldFromDataFrame(df)
        assert list(scalar_field.data.keys()) == ["singleProcessor"]

    def test_a_dataframe_with_a_processor_column_splits_by_processor(self, home):
        field = home.getEmptyField("p", flowType="incompressible")
        field.initialize(noOfProc=2)
        df = pandas.DataFrame({
            "region": ["internalField", "internalField"],
            "processorIndex": [0, 0],
            "processor": [0, 1],
            "boundary": [None, None],
            "p": [5.0, 6.0],
        })
        field.setFieldFromDataFrame(df)
        assert sorted(field.data.keys()) == ["processor0", "processor1"]

    @pytest.mark.xfail(
        strict=True,
        reason="B67: _processorToPyFoam always groups the boundary rows by "
               "a 'boundary' column, even when there are none -- a field "
               "with no boundary patches at all (exactly what getDataFrame() "
               "returns for one) has no such column, so feeding it straight "
               "back raises KeyError. See the consolidated findings issue.",
    )
    def test_a_dataframe_with_no_boundary_column_should_still_work(self, scalar_field):
        df = pandas.DataFrame({"region": ["internalField"], "processorIndex": [0], "p": [5.0]})
        scalar_field.setFieldFromDataFrame(df)

    def test_a_dataframe_with_no_boundary_column_currently_raises(self, scalar_field):
        """Characterisation of B67."""
        df = pandas.DataFrame({"region": ["internalField"], "processorIndex": [0], "p": [5.0]})
        with pytest.raises(KeyError, match="boundary"):
            scalar_field.setFieldFromDataFrame(df)

    @pytest.mark.xfail(
        strict=True,
        reason="B68: the parallel branch reads "
               "self.data[processorName]['boundaryField'] to merge in the "
               "EXISTING boundary before writing, but a freshly-built field "
               "(getEmptyField's default) is single-processor, so there is "
               "no such key yet. Setting a parallel DataFrame on it should "
               "be able to initialize the field from scratch. "
               "See the consolidated findings issue.",
    )
    def test_a_fresh_single_processor_field_can_be_set_parallel_from_scratch(self, scalar_field):
        df = pandas.DataFrame({
            "region": ["internalField", "internalField"],
            "processorIndex": [0, 0],
            "processor": [0, 1],
            "boundary": [None, None],
            "p": [5.0, 6.0],
        })
        scalar_field.setFieldFromDataFrame(df)

    def test_a_fresh_single_processor_field_currently_raises_on_parallel_data(self, scalar_field):
        """Characterisation of B68."""
        df = pandas.DataFrame({
            "region": ["internalField", "internalField"],
            "processorIndex": [0, 0],
            "processor": [0, 1],
            "boundary": [None, None],
            "p": [5.0, 6.0],
        })
        with pytest.raises(KeyError, match="processor0"):
            scalar_field.setFieldFromDataFrame(df)
