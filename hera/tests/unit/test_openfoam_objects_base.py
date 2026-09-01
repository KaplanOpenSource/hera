"""The OpenFOAM preprocessing object base class and field-file parser.

``OFObject`` is the pure-Python base of ``OFField``/``OFList`` -- dimension
bookkeeping, component-name derivation and disk I/O, no PyFoam objects
involved.  ``preprocessOFObjects/utils.py`` converts a parsed OpenFOAM field
file into a DataFrame; it is tested here against a hand-built fake that
mimics the shape PyFoam's ``ParsedParameterFile`` actually returns
(``obj['key'].val``), since the stubbed PyFoam itself is a MagicMock and
would not exercise the real branching.

B61, resolved via deletion: ``OFList._writeNew`` (and ``_updateExisting``,
which just calls it) unconditionally called ``self.getHeader()`` as its
first statement, but no class in ``OFList``'s hierarchy defined
``getHeader`` -- nor did anything set the ``self.columnNames`` it read a
few lines later. Every call raised ``AttributeError`` before any of the
scalar/vector logic ran, and nothing elsewhere in the codebase constructed
an ``OFList`` at all -- a separate dead-code-cleanup effort deleted
OFList.py entirely for exactly this reason (unreferenced, internally
broken), independently of this test-expansion effort. TestOFListIsUnusable
is dropped along with the import.
"""
import pandas
import pytest

from hera.simulations.openFoam import FIELDTYPE_SCALAR, FIELDTYPE_TENSOR, FIELDTYPE_VECTOR
from hera.simulations.openFoam.preprocessOFObjects.OFObject import OFObject
from hera.simulations.openFoam.preprocessOFObjects.utils import (
    ParsedParameterFileToDataFrame,
    extractFieldFile,
)


class _FakeVal:
    """Mimics a PyFoam parsed value: the real object of interest is ``.val``."""

    def __init__(self, val):
        self.val = val


@pytest.mark.unit
class TestOFObjectDimensions:
    def test_get_dimensions_is_the_same_format_as_openfoam(self):
        assert OFObject.getDimensions(m=1, s=-1) == "[0 1 -1 0 0 0 0]"

    def test_dimensions_str_reads_the_instance_dimensions_dict(self):
        obj = OFObject("U", "U", FIELDTYPE_VECTOR, dimensions={"m": 1, "s": -1})
        assert obj.dimensionsStr == "[0 1 -1 0 0 0 0]"

    def test_dimensions_list_defaults_missing_keys_to_zero(self):
        obj = OFObject("U", "U", FIELDTYPE_VECTOR, dimensions={"kg": 1})
        assert obj.dimensionsList == [1, 0, 0, 0, 0, 0, 0]

    def test_an_invalid_field_type_raises(self):
        with pytest.raises(ValueError, match="Field type must be"):
            OFObject("x", "x", "notAType")


@pytest.mark.unit
class TestOFObjectComponentNames:
    def test_a_scalar_field_has_one_component_named_after_the_file(self):
        obj = OFObject("p", "p", FIELDTYPE_SCALAR)
        assert obj.componentNames == ["p"]

    def test_a_vector_field_has_xyz_components(self):
        obj = OFObject("U", "U", FIELDTYPE_VECTOR)
        assert obj.componentNames == ["Ux", "Uy", "Uz"]

    def test_a_tensor_field_has_nine_components(self):
        obj = OFObject("tau", "tau", FIELDTYPE_TENSOR)
        assert obj.componentNames == [
            "tauxx", "tauxy", "tauxz",
            "tauyx", "tauyy", "tauyz",
            "tauzx", "tauzy", "tauzz",
        ]


@pytest.mark.unit
class TestOFObjectDataAccess:
    @pytest.fixture()
    def obj(self):
        o = OFObject("p", "p", FIELDTYPE_SCALAR)
        o.data = {
            "proc0": {"internalField": "IF0", "boundaryField": "BF0"},
            "proc1": {"internalField": "IF1", "boundaryField": "BF1"},
        }
        return o

    def test_internal_field_defaults_to_single_processor(self):
        o = OFObject("p", "p", FIELDTYPE_SCALAR)
        o.data = {"singleProcessor": {"internalField": "IF"}}
        assert o.internalField() == "IF"

    def test_internal_field_can_target_a_named_processor(self, obj):
        assert obj.internalField("proc1") == "IF1"

    def test_boundary_field_can_target_a_named_processor(self, obj):
        assert obj.boundaryField("proc1") == "BF1"

    def test_processors_lists_the_processor_names(self, obj):
        assert set(obj.processors) == {"proc0", "proc1"}

    def test_processor_items_pairs_names_with_their_data(self, obj):
        assert dict(obj.processorItems) == obj.data


@pytest.mark.unit
class TestOFObjectWriteToCase:
    def test_single_processor_writes_one_file_under_the_time_directory(self, tmp_path):
        obj = OFObject("p", "p", FIELDTYPE_SCALAR)
        obj.data = {"singleProcessor": "FIELD CONTENT"}
        obj.writeToCase(str(tmp_path), 0)
        written = tmp_path / "0" / "p"
        assert written.read_text() == "FIELD CONTENT"

    def test_multi_processor_writes_one_file_per_processor(self, tmp_path):
        obj = OFObject("p", "p", FIELDTYPE_SCALAR)
        obj.data = {"processor0": "P0 CONTENT", "processor1": "P1 CONTENT"}
        obj.writeToCase(str(tmp_path), "0.5")
        assert (tmp_path / "processor0" / "0.5" / "p").read_text() == "P0 CONTENT"
        assert (tmp_path / "processor1" / "0.5" / "p").read_text() == "P1 CONTENT"

    def test_a_literal_proc_wildcard_is_quoted_on_write(self, tmp_path):
        obj = OFObject("p", "p", FIELDTYPE_SCALAR)
        obj.data = {"singleProcessor": 'boundaryField { "proc.*" { type processor; } }'}
        obj.writeToCase(str(tmp_path), 0)
        written = (tmp_path / "0" / "p").read_text()
        assert '"proc.*"' in written
        assert 'boundaryField { "proc.*" { type processor; } }'.replace(
            "proc.*", '"proc.*"'
        ) in written


@pytest.mark.unit
class TestParsedParameterFileToDataFrame:
    def _fake(self, internal, boundary):
        return {"internalField": _FakeVal(internal), "boundaryField": boundary}

    def test_internal_field_rows_get_a_processor_index_and_no_boundary_columns(self):
        fake = self._fake([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], {})
        df = ParsedParameterFileToDataFrame(fake, columnNames=["x", "y", "z"])
        assert list(df["processorIndex"]) == [0, 1]
        assert list(df["region"]) == ["internalField", "internalField"]
        assert "boundary" not in df.columns

    def test_a_boundary_patch_with_a_value_gets_its_own_row(self):
        fake = self._fake(
            [[1.0]],
            {"inlet": {"value": _FakeVal([[0.0]]), "type": "fixedValue"}},
        )
        df = ParsedParameterFileToDataFrame(fake, columnNames=["p"])
        boundary_rows = df[df["region"] == "boundaryField"]
        assert list(boundary_rows["boundary"]) == ["inlet"]
        assert list(boundary_rows["type"]) == ["fixedValue"]

    def test_a_boundary_patch_without_a_value_key_still_gets_a_row_of_nans(self):
        fake = self._fake([[1.0]], {"outlet": {"type": "zeroGradient"}})
        df = ParsedParameterFileToDataFrame(fake, columnNames=["p"])
        boundary_rows = df[df["region"] == "boundaryField"]
        assert list(boundary_rows["boundary"]) == ["outlet"]
        assert boundary_rows["p"].isna().all()

    def test_a_boundary_patch_with_an_empty_value_is_skipped_entirely(self):
        fake = self._fake([[1.0]], {"empty": {"value": _FakeVal([]), "type": "fixedValue"}})
        df = ParsedParameterFileToDataFrame(fake, columnNames=["p"])
        assert (df["region"] == "boundaryField").sum() == 0

    def test_processor_patches_are_filtered_out_by_default(self):
        fake = self._fake(
            [[1.0]],
            {"procBoundary0to1": {"value": _FakeVal([[9.0]]), "type": "processor"}},
        )
        df = ParsedParameterFileToDataFrame(fake, columnNames=["p"])
        assert (df["region"] == "boundaryField").sum() == 0

    def test_an_explicit_patch_name_list_restricts_which_patches_are_kept(self):
        fake = self._fake(
            [[1.0]],
            {
                "inlet": {"value": _FakeVal([[0.0]]), "type": "fixedValue"},
                "outlet": {"value": _FakeVal([[1.0]]), "type": "fixedValue"},
            },
        )
        df = ParsedParameterFileToDataFrame(fake, columnNames=["p"], patchNameList=["inlet"])
        boundary_rows = df[df["region"] == "boundaryField"]
        assert list(boundary_rows["boundary"]) == ["inlet"]

    def test_extra_kwargs_are_stamped_onto_every_row(self):
        fake = self._fake([[1.0]], {})
        df = ParsedParameterFileToDataFrame(fake, columnNames=["p"], time="10")
        assert list(df["time"]) == ["10"]


@pytest.mark.unit
class TestExtractFieldFile:
    def test_a_parse_failure_is_wrapped_as_a_value_error(self, monkeypatch):
        import hera.simulations.openFoam.preprocessOFObjects.utils as utils_mod

        def boom(path):
            raise OSError("no such case")

        monkeypatch.setattr(utils_mod, "ParsedParameterFile", boom)
        with pytest.raises(ValueError, match="no such case"):
            extractFieldFile("/no/such/path", columnNames=["p"])

    def test_a_successful_parse_is_handed_to_the_dataframe_converter(self, monkeypatch):
        import hera.simulations.openFoam.preprocessOFObjects.utils as utils_mod

        fake = self._fake_module_level()
        monkeypatch.setattr(utils_mod, "ParsedParameterFile", lambda path: fake)
        df = extractFieldFile("/some/case/0/p", columnNames=["p"])
        assert list(df["region"]) == ["internalField"]

    @staticmethod
    def _fake_module_level():
        return {"internalField": _FakeVal([[1.0]]), "boundaryField": {}}
