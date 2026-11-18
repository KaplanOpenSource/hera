"""The OpenFOAM preprocessing object base class and field-file parser.

``OFObject`` is the pure-Python base of ``OFField``/``OFList`` -- dimension
bookkeeping, component-name derivation and disk I/O through foamlib.
``preprocessOFObjects/utils.py`` converts a parsed OpenFOAM field file into
a DataFrame; it is tested here against dicts in the shape foamlib returns:
a float for a uniform scalar, an (N,) array for a scalar list, an (N, 3)
array for a vector list, an empty array for a patch with no faces.

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
import numpy
import pandas
import pytest
from foamlib import FoamFieldFile

from hera.simulations.openFoam import FIELDTYPE_SCALAR, FIELDTYPE_TENSOR, FIELDTYPE_VECTOR
from hera.simulations.openFoam.preprocessOFObjects.OFObject import OFObject
from hera.simulations.openFoam.preprocessOFObjects.utils import (
    ParsedParameterFileToDataFrame,
    extractFieldFile,
)


def _processor(internalValue, boundaryField=None):
    return {"dimensions": [0, 2, -2, 0, 0, 0, 0], "internalField": internalValue,
            "boundaryField": boundaryField or {}}


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
        obj.data = {"singleProcessor": _processor(1.5)}
        obj.writeToCase(str(tmp_path), 0)
        written = FoamFieldFile(tmp_path / "0" / "p")
        assert (written.class_, written.object_) == ("volScalarField", "p")
        assert written.internal_field == 1.5

    def test_multi_processor_writes_one_file_per_processor(self, tmp_path):
        obj = OFObject("p", "p", FIELDTYPE_SCALAR)
        obj.data = {"processor0": _processor(1.0), "processor1": _processor(2.0)}
        obj.writeToCase(str(tmp_path), "0.5")
        assert FoamFieldFile(tmp_path / "processor0" / "0.5" / "p").internal_field == 1.0
        assert FoamFieldFile(tmp_path / "processor1" / "0.5" / "p").internal_field == 2.0

    def test_a_regex_patch_name_is_written_quoted_once(self, tmp_path):
        obj = OFObject("p", "p", FIELDTYPE_SCALAR)
        obj.data = {"singleProcessor": _processor(0.0, {'"proc.*"': {"type": "processor"}})}
        obj.writeToCase(str(tmp_path), 0)
        written = (tmp_path / "0" / "p").read_text()
        assert '"proc.*"' in written
        assert '""' not in written

    @pytest.mark.parametrize("faces", [3, 6, 9])
    def test_a_scalar_list_of_a_vector_or_tensor_length_stays_a_scalar_list(self, tmp_path, faces):
        """foamlib < 1.10 writes any 3/6/9-long array as uniform; OFObject works around it."""
        values = numpy.arange(faces, dtype=float)
        obj = OFObject("p", "p", FIELDTYPE_SCALAR)
        obj.data = {"singleProcessor": _processor(values, {"outlet": {"type": "fixedValue", "value": values}})}
        obj.writeToCase(str(tmp_path), 0)
        text = (tmp_path / "0" / "p").read_text()
        assert text.count(f"nonuniform List<scalar> {faces}(") == 2
        assert numpy.array_equal(FoamFieldFile(tmp_path / "0" / "p").internal_field, values)

    def test_a_vector_list_is_left_to_foamlib(self, tmp_path):
        obj = OFObject("U", "U", FIELDTYPE_VECTOR)
        obj.data = {"singleProcessor": _processor(numpy.eye(3))}
        obj.writeToCase(str(tmp_path), 0)
        assert "nonuniform List<vector> 3(" in (tmp_path / "0" / "U").read_text()


@pytest.mark.unit
class TestParsedParameterFileToDataFrame:
    def _fake(self, internal, boundary):
        return {"internalField": numpy.asarray(internal, dtype=float), "boundaryField": boundary}

    def test_internal_field_rows_get_a_processor_index_and_no_boundary_columns(self):
        fake = self._fake([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], {})
        df = ParsedParameterFileToDataFrame(fake, columnNames=["x", "y", "z"])
        assert list(df["processorIndex"]) == [0, 1]
        assert list(df["region"]) == ["internalField", "internalField"]
        assert "boundary" not in df.columns

    def test_a_boundary_patch_with_a_value_gets_its_own_row(self):
        fake = self._fake(
            [1.0],
            {"inlet": {"value": numpy.array([0.0]), "type": "fixedValue"}},
        )
        df = ParsedParameterFileToDataFrame(fake, columnNames=["p"])
        boundary_rows = df[df["region"] == "boundaryField"]
        assert list(boundary_rows["boundary"]) == ["inlet"]
        assert list(boundary_rows["type"]) == ["fixedValue"]

    def test_a_boundary_patch_without_a_value_key_still_gets_a_row_of_nans(self):
        fake = self._fake([1.0], {"outlet": {"type": "zeroGradient"}})
        df = ParsedParameterFileToDataFrame(fake, columnNames=["p"])
        boundary_rows = df[df["region"] == "boundaryField"]
        assert list(boundary_rows["boundary"]) == ["outlet"]
        assert boundary_rows["p"].isna().all()

    def test_a_boundary_patch_with_an_empty_value_is_skipped_entirely(self):
        fake = self._fake([1.0], {"empty": {"value": numpy.zeros(0), "type": "fixedValue"}})
        df = ParsedParameterFileToDataFrame(fake, columnNames=["p"])
        assert (df["region"] == "boundaryField").sum() == 0

    def test_processor_patches_are_filtered_out_by_default(self):
        fake = self._fake(
            [1.0],
            {"procBoundary0to1": {"value": numpy.array([9.0]), "type": "processor"}},
        )
        df = ParsedParameterFileToDataFrame(fake, columnNames=["p"])
        assert (df["region"] == "boundaryField").sum() == 0

    def test_an_explicit_patch_name_list_restricts_which_patches_are_kept(self):
        fake = self._fake(
            [1.0],
            {
                "inlet": {"value": numpy.array([0.0]), "type": "fixedValue"},
                "outlet": {"value": numpy.array([1.0]), "type": "fixedValue"},
            },
        )
        df = ParsedParameterFileToDataFrame(fake, columnNames=["p"], patchNameList=["inlet"])
        boundary_rows = df[df["region"] == "boundaryField"]
        assert list(boundary_rows["boundary"]) == ["inlet"]

    def test_extra_kwargs_are_stamped_onto_every_row(self):
        fake = self._fake([1.0], {})
        df = ParsedParameterFileToDataFrame(fake, columnNames=["p"], time="10")
        assert list(df["time"]) == ["10"]


@pytest.mark.unit
class TestExtractFieldFile:
    def test_a_missing_file_is_wrapped_as_a_value_error_naming_the_path(self, tmp_path):
        with pytest.raises(ValueError, match="no-such-field"):
            extractFieldFile(str(tmp_path / "no-such-field"), columnNames=["p"])

    def test_an_unparsable_file_is_wrapped_as_a_value_error(self, tmp_path):
        (tmp_path / "p").write_text("internalField nonuniform List<scalar> 2(1 2;\n")
        with pytest.raises(ValueError):
            extractFieldFile(str(tmp_path / "p"), columnNames=["p"])

    def test_a_parsed_file_becomes_a_dataframe(self, tmp_path):
        (tmp_path / "U").write_text(
            "internalField nonuniform List<vector> 2((1 2 3) (4 5 6));\n"
            "boundaryField { inlet { type fixedValue; value uniform (7 8 9); } }\n")
        df = extractFieldFile(str(tmp_path / "U"), columnNames=["Ux", "Uy", "Uz"])
        assert list(df["region"]) == ["internalField", "internalField", "boundaryField"]
        assert df[["Ux", "Uy", "Uz"]].values.tolist() == [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
