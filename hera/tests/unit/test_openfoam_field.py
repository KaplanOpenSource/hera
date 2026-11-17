"""OFField: initialization, reading from a case, and the DataFrame round trip.

The field data is a dict per processor in the form foamlib reads and writes
(see ``preprocessOFObjects/utils.py``), so these tests assert real stored
values.  ``test_openfoam_field_golden.py`` compares the whole read/write path
with what the PyFoam implementation produced before #1169.

Note: ``preprocessOFObjects/__init__.py`` does ``from .OFField import
OFField``, which shadows the ``OFField`` *module* attribute on the package
with the *class* of the same name.  Reaching the module needs
``sys.modules["hera.simulations.openFoam.preprocessOFObjects.OFField"]``.
"""
import numpy
import pandas
import pytest

from hera.simulations.openFoam.preprocessOFObjects.OFObjectHome import OFObjectHome


@pytest.fixture()
def home():
    return OFObjectHome()


@pytest.fixture()
def scalar_field(home):
    return home.getEmptyField("p", flowType="incompressible")


@pytest.mark.unit
class TestFieldInitialization:
    def test_a_new_field_is_uniform_zero_with_a_zero_gradient_catch_all(self, scalar_field):
        data = scalar_field.data["singleProcessor"]
        assert data["internalField"] == 0
        assert data["boundaryField"] == {'".*"': {"type": "zeroGradient"}}
        assert list(data["dimensions"]) == [0, 2, -2, 0, 0, 0, 0]

    def test_add_proc_boundary_adds_a_quoted_processor_regex(self, scalar_field):
        scalar_field.addProcBoundary()
        assert scalar_field.data["singleProcessor"]["boundaryField"]['"proc.*"'] == {"type": "processor"}

    def test_a_boundary_value_given_as_openfoam_text_is_parsed(self, home):
        field = home.getEmptyField("U", flowType="incompressible")
        field.addBoundaryField("inlet", type="fixedValue", value="uniform (5 0 0)")
        inlet = field.data["singleProcessor"]["boundaryField"]["inlet"]
        assert inlet["type"] == "fixedValue"
        assert numpy.allclose(inlet["value"], [5, 0, 0])

    def test_a_uniform_internal_value_given_as_a_list_becomes_an_array(self, home):
        field = home.getEmptyField("U", flowType="incompressible", noOfProc=2)
        field.setInternalUniformFieldValue([5, 0, 0])
        for processor in ("processor0", "processor1"):
            assert numpy.allclose(field.data[processor]["internalField"], [5, 0, 0])

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
    def test_read_from_case_wraps_a_missing_file(self, scalar_field, tmp_path):
        with pytest.raises(ValueError, match="Field not found"):
            scalar_field.readFromCase(str(tmp_path), timeStep=0, readParallel=False)

    def test_read_boundaries_from_case_wraps_a_missing_file(self, scalar_field, tmp_path):
        with pytest.raises(ValueError, match="Field not found"):
            scalar_field.readBoundariesFromCase(str(tmp_path), readParallel=False)


@pytest.mark.unit
class TestGetDataFrame:
    def test_a_single_processor_field_becomes_one_frame(self, scalar_field):
        scalar_field.data = {
            "singleProcessor": {
                "internalField": numpy.array([[1.0], [2.0]]),
                "boundaryField": {},
            }
        }
        df = scalar_field.getDataFrame()
        assert list(df["p"]) == [1.0, 2.0]
        assert list(df["region"]) == ["internalField", "internalField"]

    def test_a_multi_processor_field_concatenates_all_processors(self, scalar_field):
        scalar_field.data = {
            "processor0": {"internalField": numpy.array([[1.0]]), "boundaryField": {}},
            "processor1": {"internalField": numpy.array([[2.0]]), "boundaryField": {}},
        }
        df = scalar_field.getDataFrame()
        assert sorted(df["p"]) == [1.0, 2.0]
        assert set(df["processor"]) == {0, 1}

    def test_the_processor_number_is_parsed_from_the_key_name(self, scalar_field):
        """'processor7' -> 7, sliced from character index 9 onward."""
        scalar_field.data = {"processor7": {"internalField": numpy.array([[1.0]]), "boundaryField": {}}}
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

    def test_a_dataframe_with_no_boundary_column_should_still_work(self, scalar_field):
        df = pandas.DataFrame({"region": ["internalField"], "processorIndex": [0], "p": [5.0]})
        scalar_field.setFieldFromDataFrame(df)
        # Should not raise and result in singleProcessor data
        assert "singleProcessor" in scalar_field.data


    def test_a_fresh_single_processor_field_can_be_set_parallel_from_scratch(self, scalar_field):
        df = pandas.DataFrame({
            "region": ["internalField", "internalField"],
            "processorIndex": [0, 0],
            "processor": [0, 1],
            "boundary": [None, None],
            "p": [5.0, 6.0],
        })
        scalar_field.setFieldFromDataFrame(df)
        # Should initialize parallel data structures
        assert sorted(scalar_field.data.keys()) == ["processor0", "processor1"]

