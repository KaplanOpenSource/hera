"""datasetToOF: everything in it that does NOT need foamlib.

``hera/tests/unit/test_openfoam_dataset2of.py`` opens with a module-level
``pytest.importorskip("foamlib")``, so in this venv (foamlib is not
pip-installed and is not in requirements.txt) the whole of that file is
skipped and *none* of the module was covered -- 6% of its statements.

foamlib is only needed to read and to write OpenFOAM files.  Nine of the
module's twelve functions never touch it, and are tested here for real
against small synthetic xarray datasets:

===============================  ==========================================
``_getFoamlib``                  its ImportError contract, which is exactly
                                 what this venv exercises
``_getTransformer``              a real pyproj transformer, cached, always_xy
``_fieldMapVariables``           pure dict/list validation
``_horizontalGrid``              rectilinear / curvilinear detection
``_horizontalInterpolation``     bilinear and nearest-column, xarray + scipy
``_verticalInterpolation``       pure numpy
``interpolateDatasetToPoints``   the whole interpolation pipeline: it calls
                                 no foamlib at all
``_cellEdges``                   pure numpy
``_fieldValueTokens``            pure token building
===============================  ==========================================

The three that genuinely need it -- ``caseGeometry``,
``datasetToCaseFields`` and ``datasetToSetFieldsDict`` -- are reached here
only through the pre-foamlib part of their bodies:

* every one of them raises the documented ``ImportError`` (with the install
  command in the message) rather than an obscure ``NameError``, and that is
  asserted below because it is the behaviour this environment actually has;
* their argument validation, which happens before any file is touched, is
  asserted behind a per-test ``pytest.importorskip("foamlib")`` -- the
  established convention -- so those tests skip here and start working the
  day foamlib is installed;
* their value-level behaviour (reading the ``C`` field, writing field files,
  emitting a ``setFieldsDict``) stays in the sibling file, which is where
  the round-trip-through-foamlib assertions belong.

Rough edge characterised but not pinned as a bug: ``datasetToSetFieldsDict``
resolves foamlib on its very first line even when ``outputFile is None``,
i.e. when it will not write anything at all and is a pure region builder.
The module docstring promises only that the import is lazy with respect to
*importing hera*, which it is, so there is no single obviously-intended
behaviour to xfail against -- but it does mean the pure-Python half of that
function cannot be used without foamlib installed.

Assertions are derived from the definitions, not from output: the
interpolation tests use a field that is linear in x, y and z, which any
correct bilinear-plus-linear scheme must reproduce exactly, and the
projection test compares against a transformer built in the test.
"""
import numpy
import pytest
import xarray

from hera.measurements.GIS import ITM, WSG84
from hera.simulations.openFoam.preprocessOFObjects import datasetToOF

_X = numpy.array([180000.0, 180100.0, 180200.0])
_Y = numpy.array([660000.0, 660100.0])
_Z = numpy.array([10.0, 50.0, 100.0])

_FOAMLIB_REASON = "the value-level behaviour of this function is read/written by foamlib"


def _foamlibInstalled():
    try:
        import foamlib  # noqa: F401
    except ImportError:
        return False
    return True


def _linear(x, y, z):
    """A field that is linear in every direction, hence interpolated exactly."""
    return 2.0 * x + 3.0 * y + 5.0 * z


def _dataset():
    """A (z, y, x) dataset in ITM metres carrying three linear fields."""
    zGrid, yGrid, xGrid = numpy.meshgrid(_Z, _Y, _X, indexing="ij")
    return xarray.Dataset(
        data_vars=dict(temperature=(("z", "y", "x"), _linear(xGrid, yGrid, zGrid)),
                       u=(("z", "y", "x"), 0.5 * xGrid),
                       v=(("z", "y", "x"), -0.25 * yGrid)),
        coords=dict(x=_X, y=_Y, z=_Z),
    )


# ---------------------------------------------------------------------------
# _getFoamlib
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetFoamlib:
    def test_a_missing_foamlib_is_reported_with_the_command_that_installs_it(self):
        """This venv has no foamlib, which is precisely the branch under test."""
        if _foamlibInstalled():
            pytest.skip("foamlib is installed, so the ImportError branch cannot be reached")
        with pytest.raises(ImportError, match="pip install foamlib"):
            datasetToOF._getFoamlib()

    def test_the_two_classes_are_returned_in_file_then_fieldfile_order(self):
        pytest.importorskip("foamlib", reason="there is nothing to return without it")
        foamFile, foamFieldFile = datasetToOF._getFoamlib()
        assert foamFile.__name__ == "FoamFile"
        assert foamFieldFile.__name__ == "FoamFieldFile"


# ---------------------------------------------------------------------------
# _getTransformer
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetTransformer:
    def test_the_same_pair_of_systems_is_only_built_once(self):
        """It is asked for once per region and per processor of a case."""
        assert datasetToOF._getTransformer(ITM, WSG84) is datasetToOF._getTransformer(ITM, WSG84)

    def test_a_different_pair_is_a_different_transformer(self):
        assert datasetToOF._getTransformer(ITM, WSG84) is not datasetToOF._getTransformer(WSG84, ITM)

    def test_x_is_the_easting_and_the_longitude_never_the_authority_axis_order(self):
        """always_xy=True: Israel is near 34.8E/32.0N, so a projected ITM
        point must come back longitude-first."""
        longitude, latitude = datasetToOF._getTransformer(ITM, WSG84).transform(180000.0, 660000.0)
        assert 34.0 < longitude < 36.0
        assert 31.0 < latitude < 33.0

    def test_a_crs_given_as_a_string_is_accepted(self):
        """The codes are run through int(), so "2039" is the same as 2039."""
        fromString = datasetToOF._getTransformer(str(ITM), str(WSG84))
        assert fromString.transform(180000.0, 660000.0) == \
            datasetToOF._getTransformer(ITM, WSG84).transform(180000.0, 660000.0)


# ---------------------------------------------------------------------------
# _fieldMapVariables
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestFieldMapVariables:
    def test_a_scalar_mapping_contributes_its_one_variable(self):
        assert datasetToOF._fieldMapVariables(dict(T="temperature")) == ["temperature"]

    def test_the_numeric_components_of_a_vector_are_not_variables(self):
        assert datasetToOF._fieldMapVariables(dict(U=("u", "v", 0))) == ["u", "v"]

    def test_a_variable_used_by_two_fields_is_listed_once_in_order_of_appearance(self):
        assert datasetToOF._fieldMapVariables(dict(U=("u", "v", 0), T="u")) == ["u", "v"]

    def test_a_nine_component_tensor_mapping_is_accepted(self):
        names = tuple(f"c{index}" for index in range(9))
        assert datasetToOF._fieldMapVariables(dict(R=names)) == list(names)

    @pytest.mark.parametrize("componentCount", [2, 4, 8, 10])
    def test_only_one_three_or_nine_components_make_an_openfoam_field(self, componentCount):
        mapping = tuple(f"c{index}" for index in range(componentCount))
        with pytest.raises(ValueError, match="1 \\(scalar\\), 3 \\(vector\\) or 9 \\(tensor\\)"):
            datasetToOF._fieldMapVariables(dict(R=mapping))

    def test_a_component_that_is_neither_a_name_nor_a_number_is_rejected(self):
        with pytest.raises(ValueError, match="neither a dataset variable name nor a number"):
            datasetToOF._fieldMapVariables(dict(U=("u", "v", object())))

    def test_an_empty_field_map_refers_to_nothing(self):
        assert datasetToOF._fieldMapVariables(dict()) == []


# ---------------------------------------------------------------------------
# _horizontalGrid
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestHorizontalGrid:
    def test_two_one_dimensional_coordinates_are_broadcast_into_y_by_x_grids(self):
        grid = datasetToOF._horizontalGrid(_dataset(), "x", "y")
        assert grid["curvilinear"] is False
        assert grid["dimensions"] == ["y", "x"]
        assert grid["x"].shape == (_Y.size, _X.size)
        # every row of x repeats the x axis; every column of y repeats the y axis
        numpy.testing.assert_allclose(grid["x"], numpy.tile(_X, (_Y.size, 1)))
        numpy.testing.assert_allclose(grid["y"], numpy.tile(_Y[:, None], (1, _X.size)))

    def test_two_dimensional_coordinates_are_a_curvilinear_grid_taken_as_they_are(self):
        xGrid, yGrid = numpy.meshgrid(_X, _Y)
        dataset = xarray.Dataset(dict(XLONG=(("row", "column"), xGrid),
                                      XLAT=(("row", "column"), yGrid)))
        grid = datasetToOF._horizontalGrid(dataset, "XLONG", "XLAT")
        assert grid["curvilinear"] is True
        assert grid["dimensions"] == ["row", "column"]
        numpy.testing.assert_allclose(grid["x"], xGrid)

    def test_a_plain_data_variable_works_as_a_coordinate(self):
        """WRF stores XLONG/XLAT as data variables, not as coordinates."""
        xGrid, yGrid = numpy.meshgrid(_X, _Y)
        dataset = xarray.Dataset(dict(XLONG=(("row", "column"), xGrid),
                                      XLAT=(("row", "column"), yGrid)))
        assert "XLONG" not in dataset.coords
        assert datasetToOF._horizontalGrid(dataset, "XLONG", "XLAT")["curvilinear"] is True

    @pytest.mark.parametrize("missing", ["x", "y"])
    def test_a_name_that_is_in_neither_the_coordinates_nor_the_variables_is_named(self, missing):
        arguments = dict(x="x", y="y")
        arguments[missing] = "notThere"
        with pytest.raises(KeyError, match="notThere"):
            datasetToOF._horizontalGrid(_dataset(), arguments["x"], arguments["y"])

    def test_mixing_a_two_dimensional_and_a_one_dimensional_coordinate_is_rejected(self):
        xGrid, _ = numpy.meshgrid(_X, _Y)
        dataset = xarray.Dataset(dict(XLONG=(("row", "column"), xGrid), lat=(("y",), _Y)))
        with pytest.raises(ValueError, match="must both be 1D"):
            datasetToOF._horizontalGrid(dataset, "XLONG", "lat")

    def test_two_curvilinear_coordinates_of_different_shapes_are_rejected(self):
        xGrid, _ = numpy.meshgrid(_X, _Y)
        dataset = xarray.Dataset(dict(XLONG=(("row", "column"), xGrid),
                                      XLAT=(("r2", "c2"), numpy.zeros((_X.size, _Y.size)))))
        with pytest.raises(ValueError, match="different shapes"):
            datasetToOF._horizontalGrid(dataset, "XLONG", "XLAT")


# ---------------------------------------------------------------------------
# _horizontalInterpolation
# ---------------------------------------------------------------------------

def _smallGridded():
    """t[z, y, x] = 4z + 2y + x on a unit square, so every value is known."""
    zGrid, yGrid, xGrid = numpy.meshgrid([0.0, 1.0, 2.0], [0.0, 1.0], [0.0, 1.0], indexing="ij")
    dataset = xarray.Dataset(
        dict(t=(("z", "y", "x"), 4.0 * zGrid + 2.0 * yGrid + xGrid)),
        coords=dict(x=numpy.array([0.0, 1.0]), y=numpy.array([0.0, 1.0]),
                    z=numpy.array([0.0, 1.0, 2.0])),
    )
    return dataset, datasetToOF._horizontalGrid(dataset, "x", "y")


@pytest.mark.unit
class TestHorizontalInterpolation:
    def test_bilinear_returns_one_row_per_level_and_one_column_per_point(self):
        dataset, grid = _smallGridded()
        values = datasetToOF._horizontalInterpolation(
            dataset["t"], numpy.array([0.5, 0.0]), numpy.array([0.5, 1.0]),
            "x", "y", grid, "z", datasetToOF.HORIZONTAL_LINEAR)
        assert values.shape == (3, 2)
        # 4z + 2y + x at (0.5, 0.5) and at (0.0, 1.0)
        numpy.testing.assert_allclose(values[:, 0], [1.5, 5.5, 9.5])
        numpy.testing.assert_allclose(values[:, 1], [2.0, 6.0, 10.0])

    def test_nearest_takes_the_whole_column_of_the_closest_grid_point(self):
        dataset, grid = _smallGridded()
        values = datasetToOF._horizontalInterpolation(
            dataset["t"], numpy.array([0.4]), numpy.array([0.4]),
            "x", "y", grid, "z", datasetToOF.HORIZONTAL_NEAREST)
        # the closest column is (0, 0), whose values are 4z
        numpy.testing.assert_allclose(values[:, 0], [0.0, 4.0, 8.0])

    def test_an_array_without_a_vertical_dimension_still_gets_one_row(self):
        dataset, grid = _smallGridded()
        flat = dataset["t"].isel(z=0)
        for method in (datasetToOF.HORIZONTAL_LINEAR, datasetToOF.HORIZONTAL_NEAREST):
            values = datasetToOF._horizontalInterpolation(
                flat, numpy.array([0.0]), numpy.array([0.0]), "x", "y", grid, None, method)
            assert values.shape == (1, 1), method
            numpy.testing.assert_allclose(values, [[0.0]])

    def test_bilinear_refuses_a_curvilinear_grid_and_says_what_to_do_instead(self):
        xGrid, yGrid = numpy.meshgrid([0.0, 1.0], [0.0, 1.0])
        dataset = xarray.Dataset(dict(t=(("row", "column"), numpy.zeros((2, 2))),
                                      XLONG=(("row", "column"), xGrid),
                                      XLAT=(("row", "column"), yGrid)))
        grid = datasetToOF._horizontalGrid(dataset, "XLONG", "XLAT")
        with pytest.raises(ValueError, match="nearest"):
            datasetToOF._horizontalInterpolation(
                dataset["t"], numpy.array([0.5]), numpy.array([0.5]),
                "XLONG", "XLAT", grid, None, datasetToOF.HORIZONTAL_LINEAR)

    def test_an_unknown_method_lists_the_two_that_exist(self):
        dataset, grid = _smallGridded()
        with pytest.raises(ValueError, match="'linear' or 'nearest'"):
            datasetToOF._horizontalInterpolation(
                dataset["t"], numpy.array([0.0]), numpy.array([0.0]),
                "x", "y", grid, "z", "bilinear")


# ---------------------------------------------------------------------------
# _verticalInterpolation
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestVerticalInterpolation:
    HEIGHTS = numpy.array([[0.0, 0.0], [100.0, 100.0]])
    VALUES = numpy.array([[0.0, 10.0], [100.0, 110.0]])

    def test_a_height_between_two_levels_is_the_weighted_average(self):
        # column 0: 0 -> 100 over 0 -> 100 m, so 50 m gives 50
        # column 1: 10 -> 110 over the same 100 m, so 25 m gives 35
        result = datasetToOF._verticalInterpolation(self.VALUES, self.HEIGHTS,
                                                    numpy.array([50.0, 25.0]))
        numpy.testing.assert_allclose(result, [50.0, 35.0])

    def test_a_single_level_is_returned_as_it_is(self):
        result = datasetToOF._verticalInterpolation(self.VALUES[:1], self.HEIGHTS[:1],
                                                    numpy.array([-9999.0, 9999.0]))
        numpy.testing.assert_allclose(result, self.VALUES[0])

    def test_levels_stored_from_the_top_down_give_the_same_answer(self):
        """Pressure levels arrive in descending height order."""
        descending = datasetToOF._verticalInterpolation(
            self.VALUES[::-1, :], self.HEIGHTS[::-1, :], numpy.array([50.0, 25.0]))
        ascending = datasetToOF._verticalInterpolation(
            self.VALUES, self.HEIGHTS, numpy.array([50.0, 25.0]))
        numpy.testing.assert_allclose(descending, ascending)

    def test_a_profile_is_not_extrapolated_beyond_its_top_or_its_bottom(self):
        result = datasetToOF._verticalInterpolation(self.VALUES, self.HEIGHTS,
                                                    numpy.array([-500.0, 5000.0]))
        numpy.testing.assert_allclose(result, [self.VALUES[0, 0], self.VALUES[-1, 1]])

    def test_values_and_heights_on_different_vertical_grids_are_rejected(self):
        with pytest.raises(ValueError, match="different vertical grids"):
            datasetToOF._verticalInterpolation(self.VALUES, self.HEIGHTS[:, :1],
                                               numpy.array([1.0]))

    def test_a_column_whose_heights_go_up_and_then_down_cannot_be_searched(self):
        with pytest.raises(ValueError, match="not monotonic"):
            datasetToOF._verticalInterpolation(numpy.zeros((3, 1)),
                                               numpy.array([[0.0], [5.0], [1.0]]),
                                               numpy.array([1.0]))


# ---------------------------------------------------------------------------
# interpolateDatasetToPoints -- no foamlib anywhere in it
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestInterpolateDatasetToPoints:
    POINTS = numpy.array([[180050.0, 660050.0, 30.0],
                          [180000.0, 660000.0, 10.0],
                          [180175.0, 660025.0, 75.0]])

    def _interpolate(self, **overrides):
        arguments = dict(dataset=_dataset(), points=self.POINTS, fieldMap=dict(T="temperature"),
                         xCoordinate="x", yCoordinate="y", verticalCoordinate="z",
                         datasetCRS=ITM, caseCRS=ITM)
        arguments.update(overrides)
        return datasetToOF.interpolateDatasetToPoints(**arguments)

    def test_a_field_linear_in_x_y_and_z_is_reproduced_exactly(self):
        numpy.testing.assert_allclose(
            self._interpolate()["T"],
            _linear(self.POINTS[:, 0], self.POINTS[:, 1], self.POINTS[:, 2]))

    def test_a_scalar_field_is_one_value_per_point(self):
        assert self._interpolate()["T"].shape == (self.POINTS.shape[0],)

    def test_a_three_component_mapping_becomes_an_n_by_three_array(self):
        values = self._interpolate(fieldMap=dict(U=("u", "v", 0)))["U"]
        assert values.shape == (self.POINTS.shape[0], 3)
        numpy.testing.assert_allclose(values[:, 0], 0.5 * self.POINTS[:, 0])
        numpy.testing.assert_allclose(values[:, 1], -0.25 * self.POINTS[:, 1])
        numpy.testing.assert_allclose(values[:, 2], 0.0)

    def test_several_fields_are_returned_under_their_openfoam_names(self):
        values = self._interpolate(fieldMap=dict(T="temperature", U=("u", "v", 0)))
        assert sorted(values) == ["T", "U"]

    def test_a_dataset_with_no_vertical_dimension_needs_no_vertical_coordinate(self):
        flat = _dataset().isel(z=0, drop=True)
        values = datasetToOF.interpolateDatasetToPoints(
            dataset=flat, points=self.POINTS, fieldMap=dict(T="temperature"),
            xCoordinate="x", yCoordinate="y", datasetCRS=ITM, caseCRS=ITM)
        # the z of the point is ignored: the level is the one the dataset has
        numpy.testing.assert_allclose(
            values["T"], _linear(self.POINTS[:, 0], self.POINTS[:, 1], _Z[0]))

    def test_a_variable_without_the_vertical_dimension_keeps_its_single_level(self):
        dataset = _dataset().assign(surfaceValue=(("y", "x"), numpy.full((_Y.size, _X.size), 7.0)))
        numpy.testing.assert_allclose(self._interpolate(dataset=dataset,
                                                        fieldMap=dict(S="surfaceValue"))["S"],
                                      numpy.full(self.POINTS.shape[0], 7.0))

    def test_a_horizontal_axis_stored_backwards_is_sorted_before_interpolating(self):
        """xarray's interp needs an ascending coordinate."""
        reversed_ = _dataset().isel(y=slice(None, None, -1), x=slice(None, None, -1))
        numpy.testing.assert_allclose(self._interpolate(dataset=reversed_)["T"],
                                      self._interpolate()["T"])

    def test_the_nearest_time_step_is_selected_before_interpolating(self):
        dataset = _dataset().expand_dims(time=[0.0, 100.0]).copy()
        dataset["temperature"] = dataset["temperature"] + xarray.DataArray([0.0, 1000.0], dims="time")
        numpy.testing.assert_allclose(
            self._interpolate(dataset=dataset, time=90.0)["T"],
            _linear(self.POINTS[:, 0], self.POINTS[:, 1], self.POINTS[:, 2]) + 1000.0)

    def test_asking_for_a_time_from_a_dataset_that_has_none_is_reported(self):
        with pytest.raises(KeyError, match="no 'time' coordinate"):
            self._interpolate(time=5.0)

    def test_the_name_of_the_time_coordinate_is_an_argument(self):
        dataset = _dataset().expand_dims(Times=[0.0, 100.0]).copy()
        dataset["temperature"] = dataset["temperature"] + xarray.DataArray([0.0, 1000.0], dims="Times")
        numpy.testing.assert_allclose(
            self._interpolate(dataset=dataset, time=90.0, timeCoordinate="Times")["T"],
            _linear(self.POINTS[:, 0], self.POINTS[:, 1], self.POINTS[:, 2]) + 1000.0)

    def test_the_points_of_the_case_are_projected_into_the_system_of_the_dataset(self):
        """The dataset variable *is* the longitude, so the value at an ITM
        point must be that point's longitude -- which only holds if the point
        was projected instead of being read as a degree pair."""
        longitudes = numpy.array([34.2, 34.6, 35.0, 35.4])
        latitudes = numpy.array([31.4, 31.8, 32.2, 32.6])
        _, longitudeGrid = numpy.meshgrid(latitudes, longitudes, indexing="ij")
        dataset = xarray.Dataset(dict(longitudeValue=(("latitude", "longitude"), longitudeGrid)),
                                 coords=dict(longitude=longitudes, latitude=latitudes))

        from pyproj import Transformer
        expectedLongitude, expectedLatitude = Transformer.from_crs(
            f"EPSG:{ITM}", f"EPSG:{WSG84}", always_xy=True).transform(180000.0, 660000.0)
        assert longitudes[0] < expectedLongitude < longitudes[-1]
        assert latitudes[0] < expectedLatitude < latitudes[-1]

        values = datasetToOF.interpolateDatasetToPoints(
            dataset=dataset, points=numpy.array([[180000.0, 660000.0, 0.0]]),
            fieldMap=dict(L="longitudeValue"), xCoordinate="longitude", yCoordinate="latitude",
            datasetCRS=WSG84, caseCRS=ITM)
        numpy.testing.assert_allclose(values["L"], [expectedLongitude], rtol=1e-9)

    def test_a_terrain_following_vertical_coordinate_uses_the_height_variable(self):
        """WRF's model levels sit at a different altitude in every column."""
        levels = numpy.arange(3.0)
        heights = numpy.stack([numpy.add.outer(numpy.zeros_like(_Y), _X - _X[0]) + 20.0 * (level + 1)
                               for level in levels])
        dataset = xarray.Dataset(
            dict(height=(("level", "y", "x"), heights),
                 temperature=(("level", "y", "x"), 0.5 * heights + 7.0)),
            coords=dict(x=_X, y=_Y, level=levels))
        points = numpy.array([[_X[0], _Y[0], 30.0], [_X[0], _Y[0], 50.0]])
        values = datasetToOF.interpolateDatasetToPoints(
            dataset=dataset, points=points, fieldMap=dict(T="temperature"),
            xCoordinate="x", yCoordinate="y", verticalCoordinate="level",
            heightVariable="height", datasetCRS=ITM, caseCRS=ITM)
        numpy.testing.assert_allclose(values["T"], 0.5 * points[:, 2] + 7.0)

    def test_a_height_variable_that_is_not_in_the_dataset_is_reported(self):
        with pytest.raises(KeyError, match="notThere"):
            self._interpolate(heightVariable="notThere")

    def test_the_nearest_method_handles_a_curvilinear_grid_that_bilinear_refuses(self):
        xGrid, yGrid = numpy.meshgrid(_X, _Y)
        dataset = xarray.Dataset(
            dict(temperature=(("z", "row", "column"),
                              numpy.stack([_linear(xGrid, yGrid, level) for level in _Z])),
                 XLONG=(("row", "column"), xGrid),
                 XLAT=(("row", "column"), yGrid)),
            coords=dict(z=_Z))
        points = numpy.array([[180010.0, 660010.0, 10.0]])
        nearest = datasetToOF.interpolateDatasetToPoints(
            dataset=dataset, points=points, fieldMap=dict(T="temperature"),
            xCoordinate="XLONG", yCoordinate="XLAT", verticalCoordinate="z",
            datasetCRS=ITM, caseCRS=ITM, horizontalMethod=datasetToOF.HORIZONTAL_NEAREST)
        numpy.testing.assert_allclose(nearest["T"], _linear(_X[0], _Y[0], _Z[0]))

    @pytest.mark.parametrize("shape", [(4, 2), (3,), (2, 4)])
    def test_points_that_are_not_x_y_z_triplets_are_rejected(self, shape):
        with pytest.raises(ValueError, match=r"\(N,3\)"):
            self._interpolate(points=numpy.zeros(shape))

    def test_a_pandas_frame_of_points_is_accepted(self):
        pandas = pytest.importorskip("pandas")
        frame = pandas.DataFrame(self.POINTS, columns=["x", "y", "z"])
        numpy.testing.assert_allclose(self._interpolate(points=frame)["T"],
                                      self._interpolate()["T"])

    def test_a_variable_that_is_not_in_the_dataset_lists_the_ones_that_are(self):
        with pytest.raises(KeyError, match="temperature"):
            self._interpolate(fieldMap=dict(T="notThere"))

    def test_an_ill_formed_field_map_is_rejected_before_anything_is_interpolated(self):
        with pytest.raises(ValueError, match="components"):
            self._interpolate(fieldMap=dict(U=("u", "v")))


# ---------------------------------------------------------------------------
# _cellEdges
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestCellEdges:
    CENTRES = numpy.array([[0.0, 10.0, 30.0], [0.0, 10.0, 30.0]])

    def test_an_interior_edge_sits_halfway_between_two_centres(self):
        lower, upper = datasetToOF._cellEdges(self.CENTRES, axis=1)
        # centres 0, 10, 30 -> interior edges at 5 and 20
        numpy.testing.assert_allclose(upper[0][:2], [5.0, 20.0])
        numpy.testing.assert_allclose(lower[0][1:], [5.0, 20.0])

    def test_the_outermost_edges_mirror_the_spacing_of_the_adjacent_cell(self):
        lower, upper = datasetToOF._cellEdges(self.CENTRES, axis=1)
        # the first cell is 0 -> 5 wide inwards, so its lower edge is at -5;
        # the last is 20 -> 30, so its upper edge is at 40
        assert lower[0][0] == pytest.approx(-5.0)
        assert upper[0][-1] == pytest.approx(40.0)

    def test_an_explicit_limit_replaces_the_mirrored_outer_edge(self):
        lower, upper = datasetToOF._cellEdges(self.CENTRES, axis=1,
                                              lowerLimit=-1.0, upperLimit=99.0)
        assert lower[0][0] == pytest.approx(-1.0)
        assert upper[0][-1] == pytest.approx(99.0)
        # the interior edges are untouched
        numpy.testing.assert_allclose(lower[0][1:], [5.0, 20.0])

    def test_the_edges_have_the_shape_of_the_centres(self):
        lower, upper = datasetToOF._cellEdges(self.CENTRES, axis=1)
        assert lower.shape == upper.shape == self.CENTRES.shape

    def test_every_cell_is_covered_exactly_once_along_the_axis(self):
        lower, upper = datasetToOF._cellEdges(self.CENTRES, axis=1)
        numpy.testing.assert_allclose(lower[0][1:], upper[0][:-1])

    def test_the_other_axis_can_be_asked_for_too(self):
        centres = numpy.array([[0.0, 0.0], [10.0, 10.0], [30.0, 30.0]])
        lower, upper = datasetToOF._cellEdges(centres, axis=0)
        numpy.testing.assert_allclose(lower[:, 0], [-5.0, 5.0, 20.0])
        numpy.testing.assert_allclose(upper[:, 0], [5.0, 20.0, 40.0])

    def test_a_single_cell_axis_has_no_spacing_to_derive_an_extent_from(self):
        with pytest.raises(ValueError, match="single cell"):
            datasetToOF._cellEdges(numpy.array([[5.0]]), axis=1)

    def test_a_single_cell_axis_takes_the_extent_that_was_given(self):
        lower, upper = datasetToOF._cellEdges(numpy.array([[5.0]]), axis=1,
                                              lowerLimit=0.0, upperLimit=9.0)
        numpy.testing.assert_allclose(lower, [[0.0]])
        numpy.testing.assert_allclose(upper, [[9.0]])


# ---------------------------------------------------------------------------
# _fieldValueTokens
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestFieldValueTokens:
    ARRAYS = dict(a=numpy.arange(8.0).reshape(2, 2, 2), b=numpy.full((2, 2, 2), -1.0))

    def test_a_scalar_is_three_tokens_ending_in_the_value_at_the_index(self):
        # a[0, 0, 1] is 1.0
        assert datasetToOF._fieldValueTokens(dict(T="a"), self.ARRAYS, (0, 0, 1)) == \
            ["volScalarFieldValue", "T", 1.0]

    def test_a_vector_carries_its_three_components_as_one_list(self):
        assert datasetToOF._fieldValueTokens(dict(U=("a", "b", 3)), self.ARRAYS, (1, 1, 0)) == \
            ["volVectorFieldValue", "U", [6.0, -1.0, 3.0]]

    def test_a_nine_component_mapping_is_written_as_a_tensor(self):
        tokens = datasetToOF._fieldValueTokens(dict(R=tuple(["b"] * 9)), self.ARRAYS, (0, 0, 0))
        assert tokens[0] == "volTensorFieldValue"
        assert tokens[2] == [-1.0] * 9

    def test_a_fixed_number_does_not_have_to_come_from_the_dataset(self):
        assert datasetToOF._fieldValueTokens(dict(T=300), self.ARRAYS, (0, 0, 0)) == \
            ["volScalarFieldValue", "T", 300.0]

    def test_every_value_is_a_python_float_not_a_numpy_scalar(self):
        tokens = datasetToOF._fieldValueTokens(dict(T="a"), self.ARRAYS, (0, 0, 0))
        assert type(tokens[2]) is float

    def test_the_fields_appear_in_the_order_of_the_field_map(self):
        tokens = datasetToOF._fieldValueTokens(dict(T="a", U=("b", "b", "b")),
                                               self.ARRAYS, (0, 0, 0))
        assert tokens.index("T") < tokens.index("U")


# ---------------------------------------------------------------------------
# The three functions that need foamlib
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestTheFoamlibDependentFunctions:
    """Without foamlib these three must fail with the documented ImportError."""

    def setup_method(self):
        if _foamlibInstalled():
            pytest.skip("foamlib is installed, so the ImportError branch cannot be reached")

    def test_casegeometry_asks_for_foamlib_before_looking_at_the_case(self, tmp_path):
        with pytest.raises(ImportError, match="pip install foamlib"):
            datasetToOF.caseGeometry(str(tmp_path / "noSuchCase"))

    def test_datasettocasefields_asks_for_foamlib_before_looking_at_the_case(self, tmp_path):
        with pytest.raises(ImportError, match="pip install foamlib"):
            datasetToOF.datasetToCaseFields(caseDirectory=str(tmp_path), dataset=_dataset(),
                                            fieldMap=dict(T="temperature"),
                                            xCoordinate="x", yCoordinate="y",
                                            verticalCoordinate="z")

    def test_datasettosetfieldsdict_needs_foamlib_even_with_nothing_to_write(self):
        """Characterisation: the region list is pure Python, but foamlib is
        resolved on the first line, so `outputFile=None` needs it too."""
        with pytest.raises(ImportError, match="pip install foamlib"):
            datasetToOF.datasetToSetFieldsDict(dataset=_dataset(), fieldMap=dict(T="temperature"),
                                               xCoordinate="x", yCoordinate="y",
                                               verticalCoordinate="z", outputFile=None)


@pytest.mark.unit
class TestTheFoamlibDependentValidation:
    """Their argument validation, which runs before any file is touched.

    Skipped in this venv: foamlib is resolved on the first line of each of
    these functions, so the validation cannot be reached without it.  The
    value-level behaviour lives in test_openfoam_dataset2of.py, behind the
    same guard.
    """

    def test_a_case_without_cell_centres_names_the_utility_that_writes_them(self, tmp_path):
        pytest.importorskip("foamlib", reason=_FOAMLIB_REASON)
        with pytest.raises(FileNotFoundError, match="writeCellCentres"):
            datasetToOF.caseGeometry(str(tmp_path / "empty"))

    def test_an_unknown_extent_limit_lists_the_six_that_exist(self):
        pytest.importorskip("foamlib", reason=_FOAMLIB_REASON)
        with pytest.raises(ValueError, match="minX/maxX/minY/maxY/minZ/maxZ"):
            datasetToOF.datasetToSetFieldsDict(dataset=_dataset(), fieldMap=dict(T="temperature"),
                                               xCoordinate="x", yCoordinate="y",
                                               verticalCoordinate="z", extent=dict(minx=0.0))

    def test_a_vertical_coordinate_that_is_not_a_dimension_is_reported(self):
        pytest.importorskip("foamlib", reason=_FOAMLIB_REASON)
        with pytest.raises(KeyError, match="notADimension"):
            datasetToOF.datasetToSetFieldsDict(dataset=_dataset(), fieldMap=dict(T="temperature"),
                                               xCoordinate="x", yCoordinate="y",
                                               verticalCoordinate="notADimension")
