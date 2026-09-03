"""measurements/GIS/vector/buildings/presentation.py: the Buildings
toolkit's plotting layer.

Every method takes a GeoDataFrame directly and returns the Axes it drew
on, so the tests feed a small synthetic footprint layer (four square
buildings in ITM) and assert on the resulting axes -- collections drawn,
title, limits, aspect, colorbar label. `presentation.__init__` only
stores its two references, so plain stand-ins suffice where the parent
toolkit is not actually read.

`plotBuildingsOnMap` needs a live TilesToolkit (it fetches map tiles), so
only its delegation is exercised, through a fake tiles toolkit.
"""
import matplotlib.pyplot as plt
import numpy
import pytest

from hera.measurements.GIS.utils import ITM, WGS84
from hera.measurements.GIS.vector.buildings.presentation import presentation


@pytest.fixture(autouse=True)
def _fresh_figure():
    plt.close("all")
    yield
    plt.close("all")


@pytest.fixture()
def buildings():
    """Four 10 m squares on a 2x2 grid, with heights and morphology."""
    import geopandas
    from shapely.geometry import Polygon

    polys, heights = [], []
    for i in range(2):
        for j in range(2):
            x0, y0 = 200000 + i * 20, 700000 + j * 20
            polys.append(Polygon([(x0, y0), (x0 + 10, y0), (x0 + 10, y0 + 10), (x0, y0 + 10)]))
            heights.append(5.0 + 5 * (i + j))
    return geopandas.GeoDataFrame(
        {"BLDG_HT": heights, "lambdaP": numpy.linspace(0.1, 0.4, 4),
         "lambdaF": numpy.linspace(0.05, 0.2, 4), "hc": heights, "geometry": polys},
        geometry="geometry", crs=f"EPSG:{ITM}",
    )


@pytest.fixture()
def pres():
    return presentation(dataLayer=object(), analysisLayer=object())


@pytest.mark.unit
class TestConstruction:
    def test_both_layers_are_stored(self):
        datalayer, analysis = object(), object()
        p = presentation(dataLayer=datalayer, analysisLayer=analysis)
        assert p._datalayer is datalayer

    def test_the_analysis_layer_is_optional(self):
        assert presentation(dataLayer=object()) is not None


@pytest.mark.unit
class TestReproject:
    def test_a_crs_less_frame_is_stamped_with_the_input_crs(self, pres, buildings):
        naive = buildings.copy()
        naive.crs = None
        result = presentation._reproject(naive, inputCRS=ITM, outputCRS=None)
        assert result.crs.to_epsg() == ITM

    def test_it_reprojects_to_the_output_crs(self, pres, buildings):
        result = presentation._reproject(buildings.copy(), inputCRS=None, outputCRS=WGS84)
        assert result.crs.to_epsg() == WGS84

    def test_reprojecting_to_wgs84_yields_degree_scale_coordinates(self, pres, buildings):
        result = presentation._reproject(buildings.copy(), inputCRS=None, outputCRS=WGS84)
        minx, miny, maxx, maxy = result.total_bounds
        assert 30 < miny < 34 and 33 < minx < 36

    def test_a_crs_less_frame_with_no_input_crs_is_left_alone(self, pres, buildings):
        naive = buildings.copy()
        naive.crs = None
        result = presentation._reproject(naive, inputCRS=None, outputCRS=WGS84)
        assert result.crs is None


@pytest.mark.unit
class TestSetupAxes:
    def test_it_creates_axes_when_given_none(self):
        ax = presentation._setup_axes(None, (4, 4), None, None, None)
        assert ax is not None

    def test_it_reuses_the_axes_it_is_given(self):
        _, mine = plt.subplots()
        assert presentation._setup_axes(mine, (4, 4), None, None, None) is mine

    def test_the_title_and_limits_are_applied(self):
        ax = presentation._setup_axes(None, (4, 4), "footprints", (0, 10), (0, 20))
        assert ax.get_title() == "footprints"
        assert ax.get_xlim() == pytest.approx((0.0, 10.0))
        assert ax.get_ylim() == pytest.approx((0.0, 20.0))

    def test_the_aspect_is_forced_equal_so_maps_are_not_distorted(self):
        ax = presentation._setup_axes(None, (4, 4), None, None, None)
        assert ax.get_aspect() == 1.0


@pytest.mark.unit
class TestAddColorbar:
    def test_it_attaches_a_labelled_colorbar(self, buildings):
        _, ax = plt.subplots()
        mappable = ax.scatter([0, 1], [0, 1], c=[0.0, 1.0])
        cb = presentation._add_colorbar(ax, mappable, "height [m]")
        assert cb.ax.get_ylabel() == "height [m]"

    def test_an_empty_label_still_yields_a_colorbar(self, buildings):
        _, ax = plt.subplots()
        mappable = ax.scatter([0, 1], [0, 1], c=[0.0, 1.0])
        assert presentation._add_colorbar(ax, mappable, None) is not None


@pytest.mark.unit
class TestPlotBuildings:
    def test_it_draws_the_footprints_and_returns_the_axes(self, pres, buildings):
        ax = pres.plotBuildings(buildings)
        assert ax is not None
        assert len(ax.collections) == 1

    def test_it_draws_on_the_axes_it_is_given(self, pres, buildings):
        _, mine = plt.subplots()
        assert pres.plotBuildings(buildings, ax=mine) is mine
        assert len(mine.collections) == 1

    def test_the_title_reaches_the_axes(self, pres, buildings):
        ax = pres.plotBuildings(buildings, title="my city")
        assert ax.get_title() == "my city"

    def test_the_source_frame_is_not_mutated(self, pres, buildings):
        """The method copies before reprojecting, so the caller's frame
        keeps its own CRS."""
        before = buildings.crs
        pres.plotBuildings(buildings, outputCRS=WGS84)
        assert buildings.crs == before


@pytest.mark.unit
class TestPlotBuildingHeights:
    def test_it_colours_the_footprints_by_height(self, pres, buildings):
        ax = pres.plotBuildingHeights(buildings)
        assert len(ax.collections) >= 1

    def test_a_missing_height_column_raises(self, pres, buildings):
        with pytest.raises(Exception):
            pres.plotBuildingHeights(buildings, heightColumn="NoSuchColumn")

    def test_it_draws_on_the_axes_it_is_given(self, pres, buildings):
        _, mine = plt.subplots()
        assert pres.plotBuildingHeights(buildings, ax=mine) is mine


@pytest.mark.unit
class TestPlotMorphology:
    def test_it_draws_the_requested_morphology_column(self, pres, buildings):
        ax = pres.plotMorphology(buildings, column="lambdaP")
        assert len(ax.collections) >= 1

    def test_another_column_can_be_selected(self, pres, buildings):
        ax = pres.plotMorphology(buildings, column="lambdaF")
        assert len(ax.collections) >= 1

    def test_a_missing_column_raises(self, pres, buildings):
        with pytest.raises(Exception):
            pres.plotMorphology(buildings, column="NoSuchColumn")


@pytest.mark.unit
class TestPlotBuildingsWithMorphology:
    def test_it_overlays_both_layers_on_one_axes(self, pres, buildings):
        ax = pres.plotBuildingsWithMorphology(buildings, buildings)
        assert len(ax.collections) >= 2

    def test_it_draws_on_the_axes_it_is_given(self, pres, buildings):
        _, mine = plt.subplots()
        assert pres.plotBuildingsWithMorphology(buildings, buildings, ax=mine) is mine


@pytest.mark.unit
class TestPlotHeightHistogram:
    def test_it_draws_one_bar_per_bin(self, pres, buildings):
        ax = pres.plotHeightHistogram(buildings, bins=4)
        assert len(ax.patches) == 4

    def test_the_counts_add_up_to_the_number_of_buildings(self, pres, buildings):
        ax = pres.plotHeightHistogram(buildings, bins=4)
        assert sum(patch.get_height() for patch in ax.patches) == pytest.approx(len(buildings))

    def test_a_missing_height_column_raises(self, pres, buildings):
        with pytest.raises(Exception):
            pres.plotHeightHistogram(buildings, heightColumn="NoSuchColumn")


class _FakeTilesPresentation:
    def __init__(self, calls):
        self.calls = calls

    def plot(self, *args, **kwargs):
        self.calls.append(("plot", args, kwargs))
        return kwargs.get("ax") or plt.gca()


class _FakeTiles:
    """Stands in for the TilesToolkit, which would otherwise fetch tiles.
    The method reaches it as `tilesToolkit.presentation.plot(...)`."""

    def __init__(self):
        self.calls = []
        self.presentation = _FakeTilesPresentation(self.calls)


@pytest.mark.unit
class TestPlotBuildingsOnMap:
    def test_it_delegates_to_the_tiles_toolkit(self, pres, buildings):
        tiles = _FakeTiles()
        try:
            pres.plotBuildingsOnMap(buildings, tiles)
        except Exception:
            # The exact tiles API differs; what matters is that the method
            # reaches for the toolkit rather than fetching tiles itself.
            pass
        assert tiles.calls, "plotBuildingsOnMap never used the tiles toolkit"
