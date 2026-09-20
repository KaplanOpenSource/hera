"""simulations/utils/coordinateHandler.py: pure numpy/pandas/xarray grid
interpolation helpers, no DB or filesystem access."""
import pandas
import pytest
import xarray

from hera.simulations.utils.coordinateHandler import coordinateHandler


@pytest.fixture()
def ch():
    return coordinateHandler()


@pytest.fixture()
def slice_data():
    return pandas.DataFrame({
        "x": [0, 0, 1, 1, 0.5],
        "z": [0, 1, 0, 1, 0.5],
        "U_x": [1.0, 2.0, 3.0, 4.0, 2.5],
        "U_z": [0.1, 0.2, 0.3, 0.4, 0.25],
    })


@pytest.mark.unit
class TestRegularizeTimeStepRegular:
    def test_it_interpolates_onto_a_regular_grid(self, ch, slice_data):
        result = ch.regularizeTimeStep_Regular(
            slice_data, fieldList=["U_x"], time=None, n=(4, 4), toPandas=True, addSurface=False,
        )
        assert list(result.columns) == ["U_x"]
        assert len(result) == 16

    def test_an_unknown_field_raises_keyerror(self, ch, slice_data):
        with pytest.raises(KeyError, match="not found"):
            ch.regularizeTimeStep_Regular(
                slice_data, fieldList=["NoSuchField"], time=None, n=(3, 3), toPandas=True, addSurface=False,
            )


@pytest.mark.unit
class TestRegularizeTimeStepSigma:
    def test_it_builds_a_terrain_following_grid(self, ch, slice_data):
        data = slice_data.assign(terrain=[0.0, 0.0, 0.1, 0.1, 0.05])
        result = ch.regularizeTimeStep_Sigma(data, fieldList=["U_x"], n=3, time=0, toPandas=True)
        assert "U_x" in result.columns
        assert len(result) == 9

    def test_no_terrain_column_or_argument_raises(self, ch, slice_data):
        with pytest.raises(KeyError, match="terrain"):
            ch.regularizeTimeStep_Sigma(slice_data, fieldList=["U_x"], n=3, time=0, toPandas=True)


@pytest.mark.unit
class TestRegularizeTimeSteps:
    def test_it_dispatches_to_the_regular_method(self, ch, slice_data):
        result = ch.regularizeTimeSteps(
            slice_data, fieldList=["U_x"], coordinateType="Regular", toPandas=True, n=(3, 3), addSurface=False,
        )
        assert "U_x" in result.columns

    def test_it_runs_once_per_distinct_time_value(self, ch, slice_data):
        data = slice_data.assign(time=[0, 0, 1, 1, 1])
        result = ch.regularizeTimeSteps(
            data, fieldList=["U_x"], coordinateType="Regular", toPandas=True, n=(3, 3), addSurface=False,
        )
        assert len(result) == 2 * 9


@pytest.mark.unit
class TestCombineXarrays:
    def test_it_sums_two_arrays_on_the_merged_coords(self, ch):
        ds1 = xarray.Dataset({"v": (("x", "y"), [[1.0, 2.0], [3.0, 4.0]])}, coords={"x": [0, 1], "y": [0, 1]})
        ds2 = xarray.Dataset({"v": (("x", "y"), [[5.0, 6.0], [7.0, 8.0]])}, coords={"x": [0, 1], "y": [0, 1]})
        result = ch.combineXarrays([ds1, ds2], dims=["x", "y"])
        assert result["v"].values.tolist() == [[6.0, 8.0], [10.0, 12.0]]

    def test_more_than_two_arrays_combine_recursively(self, ch):
        ds1 = xarray.Dataset({"v": (("x", "y"), [[1.0, 1.0], [1.0, 1.0]])}, coords={"x": [0, 1], "y": [0, 1]})
        ds2 = xarray.Dataset({"v": (("x", "y"), [[2.0, 2.0], [2.0, 2.0]])}, coords={"x": [0, 1], "y": [0, 1]})
        ds3 = xarray.Dataset({"v": (("x", "y"), [[3.0, 3.0], [3.0, 3.0]])}, coords={"x": [0, 1], "y": [0, 1]})
        result = ch.combineXarrays([ds1, ds2, ds3], dims=["x", "y"])
        assert result["v"].values.tolist() == [[6.0, 6.0], [6.0, 6.0]]
