"""hera/measurements/GIS/raster/topography.py::topographyAnalysis -- the
statistics layer, which ``test_gis_raster_topography_more.py`` explicitly
leaves out of scope.

Covered here
------------
``topographyAnalysis.__init__``, the ``datalayer`` attribute, and
``calculateStastics`` (sic).

How the inputs were produced
----------------------------
``calculateStastics`` indexes its argument with ``['X']``, ``['Y']`` and
``['Elevation']`` and then uses ``.values``, so anything supporting those --
a pandas DataFrame or an xarray Dataset -- fits.  The tabular form is what the
repo's own test factory produces (``hera/tests/unit/_factories.py``
``elevation_grid``, whose columns are exactly ``X``/``Y``/``Elevation``), so
that shape is used here.  ``_grid`` below lays out a 3x2 grid whose elevation
is ``100 + 50*i + 10*j``: the minimum and maximum sit at known, *distinct*
corners so a swapped argmin/argmax or a swapped X/Y readout would show, and
mean and standard deviation are computed independently with ``statistics``
from the same six numbers rather than with numpy.

Bug pinned (xfail(strict) for the intended behaviour plus a passing
characterisation of what actually happens)
--------------------------------------------------------------------------
* **B241** -- ``calculateStastics`` and the toolkit's only elevation producer
  disagree on their column names.  ``TopographyToolkit.getElevation`` returns
  an ``xarray.Dataset`` carrying ``lat``, ``lon`` and an ``elevation``
  coordinate; ``calculateStastics`` reads ``X``, ``Y`` and ``Elevation``.
  Feeding the documented producer's output to the documented consumer raises
  ``KeyError: 'X'``.  The docstring's own worked example cannot be followed
  either: it calls ``tk.getDomainElevation(...)``, a method that does not
  exist anywhere in the codebase.

Deliberately not tested, and why
--------------------------------
* Reaching ``getElevation`` through a real raster: the shape of its output is
  what matters for B241, not its values, and the values are already covered in
  ``test_gis_raster_topography_more.py``.  The dataset used for the pinning
  tests is therefore built by calling the toolkit's *own* ``create_xarray``
  and then attaching an ``elevation`` coordinate exactly as ``getElevation``
  does, which needs no raster and cannot drift from the producer.
"""
import statistics

import numpy
import pandas
import pytest

from hera import toolkitHome
from hera.measurements.GIS.raster.topography import TopographyToolkit, topographyAnalysis
from hera.measurements.GIS.utils import WSG84

# 3 columns of X, 2 rows of Y; elevation = 100 + 50*i + 10*j.
_XS = (1.0, 2.0, 3.0)
_YS = (10.0, 20.0)
_ELEVATIONS = [100.0 + 50.0 * i + 10.0 * j
               for j, _ in enumerate(_YS) for i, _ in enumerate(_XS)]


def _grid():
    rows = []
    for j, y in enumerate(_YS):
        for i, x in enumerate(_XS):
            rows.append((x, y, 100.0 + 50.0 * i + 10.0 * j))
    return pandas.DataFrame(rows, columns=["X", "Y", "Elevation"])


@pytest.fixture()
def topo(unit_toolkit_factory):
    return unit_toolkit_factory(toolkitHome.GIS_RASTER_TOPOGRAPHY)


@pytest.fixture()
def stats(topo):
    return topo.analysis.calculateStastics(_grid())


@pytest.mark.unit
class TestAnalysisLayerWiring:
    def test_the_toolkit_builds_a_statistics_layer(self, topo):
        assert isinstance(topo.analysis, topographyAnalysis)

    def test_the_layer_points_back_at_the_toolkit(self, topo):
        assert topo.analysis.datalayer is topo

    def test_it_can_also_be_constructed_directly(self, topo):
        assert topographyAnalysis(datalayer=topo).datalayer is topo

    def test_the_datalayer_is_a_plain_attribute_with_a_none_class_default(self):
        """Not a property: the class attribute is None until __init__ runs."""
        assert topographyAnalysis.datalayer is None


@pytest.mark.unit
class TestCalculateStatisticsBounds:
    def test_the_x_bounds_come_from_the_x_column(self, stats):
        assert (stats["xmin"], stats["xmax"]) == (min(_XS), max(_XS))

    def test_the_y_bounds_come_from_the_y_column(self, stats):
        assert (stats["ymin"], stats["ymax"]) == (min(_YS), max(_YS))

    def test_x_and_y_are_not_swapped(self, stats):
        """The two axes have deliberately disjoint ranges."""
        assert stats["xmax"] < stats["ymin"]

    def test_the_size_is_the_area_of_the_bounding_box(self, stats):
        expected = (max(_XS) - min(_XS)) * (max(_YS) - min(_YS))
        assert stats["size"] == pytest.approx(expected)
        assert stats["size"] == pytest.approx(2.0 * 10.0)

    def test_a_single_point_domain_has_zero_size(self, topo):
        one = pandas.DataFrame([(5.0, 6.0, 7.0)], columns=["X", "Y", "Elevation"])
        assert topo.analysis.calculateStastics(one)["size"] == 0.0


@pytest.mark.unit
class TestCalculateStatisticsMoments:
    def test_the_mean_is_the_arithmetic_mean_of_the_elevations(self, stats):
        assert stats["mean"] == pytest.approx(statistics.fmean(_ELEVATIONS))

    def test_the_standard_deviation_is_the_population_one(self, stats):
        """numpy's default ddof=0, i.e. pstdev rather than stdev."""
        assert stats["std"] == pytest.approx(statistics.pstdev(_ELEVATIONS))

    def test_the_standard_deviation_is_not_the_sample_one(self, stats):
        assert stats["std"] != pytest.approx(statistics.stdev(_ELEVATIONS))

    def test_a_flat_domain_has_zero_standard_deviation(self, topo):
        flat = _grid().assign(Elevation=42.0)
        result = topo.analysis.calculateStastics(flat)
        assert result["std"] == 0.0
        assert result["mean"] == 42.0


@pytest.mark.unit
class TestCalculateStatisticsExtrema:
    def test_the_maximum_is_the_largest_elevation(self, stats):
        assert stats["domainmax"] == max(_ELEVATIONS)

    def test_the_minimum_is_the_smallest_elevation(self, stats):
        assert stats["domainmin"] == min(_ELEVATIONS)

    def test_the_maximum_is_located_at_its_own_x_and_y(self, stats):
        """100 + 50*i + 10*j peaks at the last x and the last y."""
        assert stats["domainmaxlocation"] == (max(_XS), max(_YS))

    def test_the_minimum_is_located_at_its_own_x_and_y(self, stats):
        assert stats["domainminlocation"] == (min(_XS), min(_YS))

    def test_the_two_locations_are_reported_as_x_y_pairs_not_y_x(self, stats):
        """X spans 1..3 and Y spans 10..20, so a swap is unambiguous."""
        for key in ("domainmaxlocation", "domainminlocation"):
            x, y = stats[key]
            assert min(_XS) <= x <= max(_XS)
            assert min(_YS) <= y <= max(_YS)

    def test_the_extrema_locations_follow_the_data_when_it_is_reshuffled(self, topo):
        """A hand-placed peak must be found where it was put."""
        grid = _grid()
        grid.loc[1, "Elevation"] = 9999.0
        result = topo.analysis.calculateStastics(grid)
        assert result["domainmax"] == 9999.0
        assert result["domainmaxlocation"] == (grid.loc[1, "X"], grid.loc[1, "Y"])


@pytest.mark.unit
class TestCalculateStatisticsContract:
    def test_it_returns_exactly_the_eleven_documented_keys(self, stats):
        assert set(stats) == {
            "xmin", "xmax", "ymin", "ymax", "size",
            "mean", "std", "domainmax", "domainmaxlocation",
            "domainmin", "domainminlocation",
        }

    def test_it_does_not_mutate_the_input(self, topo):
        grid = _grid()
        before = grid.copy()
        topo.analysis.calculateStastics(grid)
        pandas.testing.assert_frame_equal(grid, before)

    def test_a_frame_without_the_elevation_column_is_rejected(self, topo):
        with pytest.raises(KeyError):
            topo.analysis.calculateStastics(_grid().rename(columns={"Elevation": "z"}))


# --------------------------------------------------------------------------
# B241: the producer and the consumer disagree on column names
# --------------------------------------------------------------------------

def _get_elevation_shaped_dataset(topo):
    """What ``getElevation`` returns, built without needing a raster.

    ``getElevation`` is ``create_xarray(...)`` followed by
    ``assign_coords(elevation=(('i','j'), <heights>))``; the heights come from
    ``getPointListElevation``, which is the only part that needs a file.  The
    two structural steps are reproduced verbatim here, so this cannot drift
    from the producer's shape.
    """
    dataset = topo.create_xarray(35.0, 32.0, 35.01, 32.01, dxdy=300, inputCRS=WSG84)
    shape = dataset["lat"].shape
    heights = numpy.arange(numpy.prod(shape), dtype=float).reshape(shape)
    return dataset.assign_coords(elevation=(("i", "j"), heights))


@pytest.mark.unit
class TestStatisticsOfTheToolkitsOwnElevationOutput:
    def test_the_producer_emits_lowercase_lat_lon_elevation(self, topo):
        """Characterisation of B241: the producer's side of the mismatch."""
        dataset = _get_elevation_shaped_dataset(topo)
        assert {"lat", "lon"} <= set(dataset.variables)
        assert "elevation" in dataset.coords
        for name in ("X", "Y", "Elevation"):
            assert name not in dataset.variables

    @pytest.mark.xfail(
        strict=True,
        reason="B241: calculateStastics reads ['X'], ['Y'] and ['Elevation'], "
               "but TopographyToolkit.getElevation -- the toolkit's only "
               "elevation producer -- returns lat/lon plus an 'elevation' "
               "coordinate, so the documented consumer cannot read the "
               "documented producer's output and raises KeyError: 'X'.  The "
               "method's own docstring example is unusable too: it calls "
               "tk.getDomainElevation(), which does not exist. "
               "See the consolidated findings issue.",
    )
    def test_the_statistics_can_be_taken_of_a_getelevation_result(self, topo):
        stats = topo.analysis.calculateStastics(_get_elevation_shaped_dataset(topo))
        assert "mean" in stats

    def test_feeding_the_producers_output_in_raises_a_keyerror(self, topo):
        """Characterisation of B241."""
        with pytest.raises(KeyError) as excinfo:
            topo.analysis.calculateStastics(_get_elevation_shaped_dataset(topo))
        assert "X" in str(excinfo.value)

    def test_the_documented_worked_example_names_a_method_that_does_not_exist(self):
        """Characterisation of B241: the docstring example cannot be run."""
        assert "getDomainElevation" in topographyAnalysis.calculateStastics.__doc__
        assert not hasattr(TopographyToolkit, "getDomainElevation")

    def test_renaming_the_variables_makes_the_same_dataset_work(self, topo):
        """Characterisation of B241: nothing else is wrong -- only the names.

        The same dataset, with lat/lon/elevation renamed to X/Y/Elevation,
        goes straight through, which locates the defect precisely at the
        naming boundary rather than in the statistics themselves.
        """
        dataset = _get_elevation_shaped_dataset(topo)
        renamed = dataset.rename({"lat": "Y", "lon": "X", "elevation": "Elevation"})
        stats = topo.analysis.calculateStastics(renamed)
        assert stats["domainmin"] == 0.0
        assert stats["domainmax"] == float(numpy.prod(dataset["lat"].shape) - 1)
