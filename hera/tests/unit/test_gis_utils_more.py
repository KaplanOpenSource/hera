"""hera/measurements/GIS/utils.py::stlFactory.rasterizePandas -- the one
function in that module the existing GIS unit tests leave uncovered.

``rasterizeGeopandas`` (the contour/LineString path) is tested in
``test_gis_utils_stlfactory.py``; ``rasterizePandas`` is its *tabular*
sibling, for topography that arrives as x/y/height columns rather than as
iso-height lines, and it delegates the interpolation to
``coordinateHandler.regularizeTimeSteps``.

How the inputs were produced
----------------------------
``_ramp_table`` builds a regular x/y point table whose height is a plain
linear function of x, so the interpolated field is reproduced exactly by the
linear ``griddata`` inside ``regularizeTimeSteps`` -- every expected value
below has a closed form (``height == x``) rather than being a recorded
number.  A ramp in x only (not in y) also makes an accidental axis transpose
visible: a transposed grid would vary along the wrong index.

Bugs pinned (each has an xfail(strict) for the intended behaviour plus a
passing characterisation of what actually happens)
--------------------------------------------------------------------------
The returned dict promises three co-registered 2-D arrays -- ``x``, ``y`` and
``height`` -- and ``vectorToSTL`` feeds them straight to ``rasterToSTL`` as
``grid_x``/``grid_y``/``grid_z``.  They are not co-registered, because the
coordinate arrays and the height array are built with two different
conventions:

* **B239** -- the *shapes* disagree whenever the span is not an exact multiple
  of ``dxdy``.  ``grid_x`` comes from ``numpy.mgrid[xmin:xmax:dxdy]``, which
  *rounds the point count up*; ``Nx`` is ``int((xmax - xmin) / dxdy)``, which
  *truncates*.  A 100 m span at ``dxdy=30`` gives a 4-point coordinate axis
  against a 3-point height axis, so ``rasterToSTL`` then reads coordinates
  for cells that have no height (or vice versa).
* **B240** -- even when the shapes match, the *coordinates* are wrong.
  ``numpy.mgrid`` with a step excludes the endpoint, while
  ``regularizeTimeSteps`` interpolates onto ``numpy.mgrid[...:complex(n)]``,
  which *includes* both endpoints.  For a 100 m span at ``dxdy=50`` the
  heights are sampled at x = 0 and x = 100 but labelled x = 0 and x = 50, so
  every STL vertex except the first is placed at the wrong location.

Deliberately not tested, and why
--------------------------------
* The ``dask.dataframe.DataFrame`` dispatch in ``vectorToSTL`` that would
  route to this function: ``dask`` is never imported in the module, which is
  B74, already pinned in ``test_gis_utils_stlfactory.py``.  ``rasterizePandas``
  is therefore called directly here rather than through ``vectorToSTL``.
* ``_organizeGrid``'s NaN edge filling: it has its own coverage through
  ``rasterizeGeopandas``; only its interaction with the xarray DataArray this
  function hands it is checked below.
"""
import numpy
import pandas
import pytest

from hera.measurements.GIS.utils import stlFactory


def _ramp_table(xmin=0.0, xmax=100.0, ymin=0.0, ymax=100.0, n=6):
    """A regular point cloud whose height equals its x coordinate."""
    xs = numpy.linspace(xmin, xmax, n)
    ys = numpy.linspace(ymin, ymax, n)
    xx, yy = numpy.meshgrid(xs, ys)
    return pandas.DataFrame(
        {"x": xx.ravel(), "y": yy.ravel(), "height": xx.ravel()}
    )


@pytest.fixture()
def factory():
    return stlFactory()


@pytest.mark.unit
class TestRasterizePandasOutputStructure:
    def test_it_returns_the_three_keys_rastertostl_expects(self, factory):
        result = factory.rasterizePandas(_ramp_table(), dxdy=50.0)
        assert set(result) == {"x", "y", "height"}

    def test_every_value_is_a_two_dimensional_numpy_array(self, factory):
        result = factory.rasterizePandas(_ramp_table(), dxdy=50.0)
        for key in ("x", "y", "height"):
            assert isinstance(result[key], numpy.ndarray)
            assert result[key].ndim == 2

    def test_the_coordinate_axes_start_at_the_data_minimum(self, factory):
        result = factory.rasterizePandas(_ramp_table(), dxdy=50.0)
        assert result["x"][0, 0] == pytest.approx(0.0)
        assert result["y"][0, 0] == pytest.approx(0.0)

    def test_the_coordinate_step_is_dxdy(self, factory):
        result = factory.rasterizePandas(_ramp_table(), dxdy=50.0)
        assert result["x"][1, 0] - result["x"][0, 0] == pytest.approx(50.0)

    def test_a_finer_resolution_gives_a_larger_grid(self, factory):
        coarse = factory.rasterizePandas(_ramp_table(), dxdy=50.0)
        fine = factory.rasterizePandas(_ramp_table(), dxdy=25.0)
        assert fine["height"].size > coarse["height"].size

    def test_the_default_resolution_is_fifty(self, factory):
        default = factory.rasterizePandas(_ramp_table())
        explicit = factory.rasterizePandas(_ramp_table(), dxdy=50.0)
        assert default["height"].shape == explicit["height"].shape


@pytest.mark.unit
class TestRasterizePandasColumnNames:
    def test_the_column_names_are_configurable(self, factory):
        renamed = _ramp_table().rename(
            columns={"x": "easting", "y": "northing", "height": "z"}
        )
        result = factory.rasterizePandas(
            renamed, dxdy=50.0, xColumn="easting", yColumn="northing", heightColumn="z"
        )
        assert result["height"].shape == (2, 2)

    def test_a_missing_height_column_is_reported_by_name(self, factory):
        with pytest.raises(KeyError, match="no_such_column"):
            factory.rasterizePandas(
                _ramp_table(), dxdy=50.0, heightColumn="no_such_column"
            )

    def test_the_height_column_property_is_not_consulted(self, factory):
        """Unlike rasterizeGeopandas, this method takes the column as an
        argument and ignores ``heightColumnsNames`` entirely."""
        factory.heightColumnsNames = "IRRELEVANT"
        result = factory.rasterizePandas(_ramp_table(), dxdy=50.0)
        assert result["height"].shape == (2, 2)


@pytest.mark.unit
class TestRasterizePandasInterpolation:
    def test_the_interpolated_field_reproduces_the_linear_ramp(self, factory):
        """height == x by construction, and linear griddata is exact on it."""
        result = factory.rasterizePandas(_ramp_table(), dxdy=50.0)
        # regularizeTimeSteps samples the *endpoint-inclusive* axis 0..100.
        assert result["height"] == pytest.approx(
            numpy.array([[0.0, 0.0], [100.0, 100.0]])
        )

    def test_the_ramp_varies_along_the_first_index_not_the_second(self, factory):
        """A transposed grid would vary along axis 1 instead."""
        result = factory.rasterizePandas(_ramp_table(), dxdy=50.0)
        assert result["height"][0, 0] != result["height"][1, 0]
        assert result["height"][0, 0] == pytest.approx(result["height"][0, 1])

    def test_a_constant_surface_comes_back_constant(self, factory):
        flat = _ramp_table().assign(height=42.0)
        result = factory.rasterizePandas(flat, dxdy=50.0)
        assert result["height"] == pytest.approx(numpy.full((2, 2), 42.0))

    def test_no_nan_survives_the_edge_filling(self, factory):
        result = factory.rasterizePandas(_ramp_table(), dxdy=25.0)
        assert not numpy.isnan(result["height"]).any()


@pytest.mark.unit
class TestRasterizePandasGridShapeMismatch:
    @pytest.mark.xfail(
        strict=True,
        reason="B239: rasterizePandas builds its coordinate axes with "
               "numpy.mgrid[xmin:xmax:dxdy] (point count rounded UP) but sizes "
               "the interpolated height grid with Nx=int((xmax-xmin)/dxdy) "
               "(truncated).  Any span that is not an exact multiple of dxdy "
               "yields co-ordinate and height arrays of different shape, "
               "which rasterToSTL then indexes as if they matched. "
               "See the consolidated findings issue.",
    )
    def test_the_coordinate_and_height_arrays_always_have_the_same_shape(self, factory):
        result = factory.rasterizePandas(_ramp_table(), dxdy=30.0)
        assert result["x"].shape == result["height"].shape
        assert result["y"].shape == result["height"].shape

    def test_a_non_multiple_resolution_produces_mismatched_shapes(self, factory):
        """Characterisation of B239: 100 / 30 -> 4 coordinate rows, 3 height rows."""
        result = factory.rasterizePandas(_ramp_table(), dxdy=30.0)
        assert result["x"].shape == (4, 4)
        assert result["height"].shape == (3, 3)

    def test_an_exact_multiple_happens_to_agree(self, factory):
        """Characterisation of B239: the defect hides on tidy numbers."""
        result = factory.rasterizePandas(_ramp_table(), dxdy=50.0)
        assert result["x"].shape == result["height"].shape == (2, 2)

    def test_the_mismatch_silently_truncates_the_stl_instead_of_erroring(self, factory):
        """Characterisation of B239: the downstream consequence, stated.

        ``rasterToSTL`` loops over ``grid_z.shape``, so a *smaller* height
        grid quietly drops the outermost coordinate row and column rather
        than raising: at ``dxdy=30`` the coordinate axis reaches x = 90 but no
        vertex past x = 60 is ever emitted.  The failure mode is a silently
        undersized solid, which is worse than an exception.
        """
        result = factory.rasterizePandas(_ramp_table(), dxdy=30.0)
        stl = factory.rasterToSTL(
            grid_x=result["x"], grid_y=result["y"],
            grid_z=result["height"], solidName="Topography",
        )
        vertex_x = [
            float(line.split()[1])
            for line in stl.splitlines()
            if line.strip().startswith("vertex")
        ]
        assert result["x"].max() == pytest.approx(90.0)
        assert max(vertex_x) == pytest.approx(60.0)


@pytest.mark.unit
class TestRasterizePandasCoordinateRegistration:
    @pytest.mark.xfail(
        strict=True,
        reason="B240: the coordinate axes come from numpy.mgrid with a *step*, "
               "which excludes the endpoint, while the heights are "
               "interpolated onto numpy.mgrid with a *count* "
               "(complex(n)), which includes it.  The two grids therefore "
               "sample different locations even when their shapes agree, so "
               "every height after the first is labelled with the wrong "
               "coordinate. See the consolidated findings issue.",
    )
    def test_each_height_is_labelled_with_the_coordinate_it_was_sampled_at(self, factory):
        """height == x, so a correctly registered grid has height == grid_x."""
        result = factory.rasterizePandas(_ramp_table(), dxdy=50.0)
        assert result["height"] == pytest.approx(result["x"])

    def test_the_coordinate_axis_stops_one_step_short_of_the_data(self, factory):
        """Characterisation of B240."""
        result = factory.rasterizePandas(_ramp_table(), dxdy=50.0)
        assert result["x"].max() == pytest.approx(50.0)   # data span is 0..100
        assert result["y"].max() == pytest.approx(50.0)

    def test_the_height_axis_spans_the_full_data_range(self, factory):
        """Characterisation of B240: the other half of the disagreement."""
        result = factory.rasterizePandas(_ramp_table(), dxdy=50.0)
        assert result["height"].max() == pytest.approx(100.0)

    def test_the_last_cell_is_off_by_a_full_span(self, factory):
        """Characterisation of B240: quantified.

        The height at index [1, 0] was interpolated at x = 100 but is labelled
        x = 50, a 50 m error -- one whole dxdy, growing with the grid.
        """
        result = factory.rasterizePandas(_ramp_table(), dxdy=50.0)
        assert result["height"][1, 0] - result["x"][1, 0] == pytest.approx(50.0)
