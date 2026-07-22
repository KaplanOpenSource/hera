"""
Comprehensive pytest tests for the DemographyToolkit.

Tests cover:
- Toolkit initialization and properties
- setDefaultDirectory
- analysis.calculatePopulationInPolygon (various polygon scenarios)
- analysis.createNewArea (with different inputs and save modes)
- presentation.plotPopulationDensity (with density_units, CRS options)
- presentation.plotPopulation (absolute counts)
- presentation.plotPopulationByType (subplot grid)
- presentation.plotPopulationInPolygon (intersection visualization)
- presentation.plotArea (custom area with annotation)
- presentation.plotPopulationOnMap (map overlay — mocked tiles)

Data is loaded into the test project via ``conftest.hera_test_project``.
The ``demo_toolkit`` fixture (session-scoped, from conftest) provides
a real ``DemographyToolkit`` backed by MongoDB.

Synthetic GeoDataFrame fixtures allow tests to run even without the full
test data repository.
"""

import os
import tempfile

try:
    from pyogrio.errors import DataSourceError as _PyogrioDataSourceError
except ImportError:
    _PyogrioDataSourceError = None

import pytest

gpd = pytest.importorskip("geopandas", reason="geopandas not installed")
matplotlib = pytest.importorskip("matplotlib", reason="matplotlib not installed")
plt = pytest.importorskip("matplotlib.pyplot", reason="matplotlib not installed")

import numpy as np
from shapely.geometry import Polygon, box

from hera.measurements.GIS import WSG84, ITM
from hera.toolkit import TOOLKIT_SAVEMODE_NOSAVE
from hera.datalayer.document.metadataDocument import nonDBMetadataFrame
from hera.utils.unitHandler import ureg

matplotlib.use("Agg")  # non-interactive backend for CI


# ---------------------------------------------------------------------------
# Synthetic data fixtures (no DB required)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def synthetic_gdf():
    """Small synthetic population GeoDataFrame in ITM (EPSG:2039)."""
    geometry = [
        box(170000, 660000, 171000, 661000),  # 1 km² block
        box(171000, 660000, 172000, 661000),
        box(170000, 661000, 171000, 662000),
        box(171000, 661000, 172000, 662000),
    ]
    return gpd.GeoDataFrame(
        {
            "total_pop": [1000, 2000, 500, 1500],
            "age_0_14": [200, 400, 100, 300],
            "age_15_19": [100, 200, 50, 150],
            "age_20_29": [250, 500, 125, 375],
            "age_30_64": [350, 700, 175, 525],
            "age_65_up": [100, 200, 50, 150],
            "geometry": geometry,
        },
        crs="EPSG:2039",
    )


@pytest.fixture(scope="module")
def synthetic_gdf_wgs84():
    """Small synthetic population GeoDataFrame in WGS84 (EPSG:4326)."""
    geometry = [
        box(34.75, 32.05, 34.76, 32.06),
        box(34.76, 32.05, 34.77, 32.06),
        box(34.75, 32.06, 34.76, 32.07),
    ]
    return gpd.GeoDataFrame(
        {
            "total_pop": [800, 1200, 600],
            "age_0_14": [160, 240, 120],
            "geometry": geometry,
        },
        crs="EPSG:4326",
    )


@pytest.fixture(scope="module")
def large_polygon_itm():
    """Polygon that covers the full synthetic ITM data."""
    return box(169500, 659500, 172500, 662500)


@pytest.fixture(scope="module")
def small_polygon_itm():
    """Polygon that partially overlaps synthetic ITM data (only intersects block 1)."""
    return box(170100, 660100, 170900, 660900)


# ---------------------------------------------------------------------------
# DB-backed fixtures (require hera_test_project from conftest)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def population_gdf(demo_toolkit):
    """Load the population GeoDataFrame through the project's datasource."""
    try:
        gdf = demo_toolkit.getDataSourceData("lamas_population")
    except FileNotFoundError as exc:
        pytest.skip(f"lamas_population data file missing in TEST_HERA: {exc}")
    except Exception as exc:
        # pyogrio raises DataSourceError (RuntimeError-derived) when the
        # shapefile is absent. Match only "no such file" — corruption or
        # format errors should still surface as real failures.
        if _PyogrioDataSourceError is not None and isinstance(exc, _PyogrioDataSourceError) \
                and "No such file" in str(exc):
            pytest.skip(f"lamas_population data file missing in TEST_HERA: {exc}")
        raise
    if gdf is None:
        pytest.skip("lamas_population datasource not loaded in project")
    return gdf


# ===========================================================================
# 1. Toolkit initialization and properties
# ===========================================================================

@pytest.mark.integration
class TestDemographyToolkitInit:
    """Test toolkit construction and properties (requires DB)."""

    def test_analysis_property(self, demo_toolkit):
        assert demo_toolkit.analysis is not None

    def test_presentation_property(self, demo_toolkit):
        assert demo_toolkit.presentation is not None

    def test_population_types(self, demo_toolkit):
        pt = demo_toolkit.populationTypes
        assert isinstance(pt, dict)
        assert "All" in pt
        assert pt["All"] == "total_pop"
        assert "Children" in pt
        assert pt["Children"] == "age_0_14"
        assert "Elderly" in pt
        assert pt["Elderly"] == "age_65_up"


# ===========================================================================
# 2. setDefaultDirectory
# ===========================================================================

@pytest.mark.integration
class TestSetDefaultDirectory:
    def test_creates_and_sets_path(self, demo_toolkit):
        with tempfile.TemporaryDirectory() as tmpdir:
            test_dir = os.path.join(tmpdir, "demography_test_folder")
            demo_toolkit.setDefaultDirectory(test_dir, create=True)
            assert os.path.exists(test_dir)
            assert os.path.abspath(test_dir) == demo_toolkit.filesDirectory

    def test_raises_if_not_found_and_create_false(self, demo_toolkit):
        with pytest.raises(NotADirectoryError):
            demo_toolkit.setDefaultDirectory("/tmp/__nonexistent_xyz__", create=False)


# ===========================================================================
# 3. analysis.calculatePopulationInPolygon
# ===========================================================================

@pytest.mark.integration
class TestCalculatePopulationInPolygon:
    def test_basic(self, demo_toolkit, population_gdf):
        """Polygon enclosing a census area returns non-empty result."""
        base_polygon = population_gdf.iloc[0].geometry
        minx, miny, maxx, maxy = base_polygon.bounds
        buf = 10
        test_polygon = Polygon([
            (minx - buf, miny - buf),
            (maxx + buf, miny - buf),
            (maxx + buf, maxy + buf),
            (minx - buf, maxy + buf),
            (minx - buf, miny - buf),
        ])

        result = demo_toolkit.analysis.calculatePopulationInPolygon(
            shapelyPolygon=test_polygon,
            dataSourceOrData="lamas_population",
        )
        assert not result.empty
        assert "geometry" in result.columns
        assert "areaFraction" in result.columns

    def test_partial_intersection(self, demo_toolkit, population_gdf):
        """Polygon spanning two census areas returns weighted populations."""
        poly1 = population_gdf.iloc[0].geometry
        poly2 = population_gdf.iloc[1].geometry
        intersection_area = poly1.union(poly2).centroid.buffer(20)

        result = demo_toolkit.analysis.calculatePopulationInPolygon(
            shapelyPolygon=intersection_area,
            dataSourceOrData="lamas_population",
        )
        assert not result.empty
        assert len(result) >= 1
        assert "areaFraction" in result.columns

    def test_outside_bounds(self, demo_toolkit, population_gdf):
        """Polygon far from data returns empty result."""
        minx, miny, maxx, maxy = population_gdf.total_bounds
        far_polygon = Polygon([
            (maxx + 10000, maxy + 10000),
            (maxx + 10100, maxy + 10000),
            (maxx + 10100, maxy + 10100),
            (maxx + 10000, maxy + 10100),
            (maxx + 10000, maxy + 10000),
        ])
        result = demo_toolkit.analysis.calculatePopulationInPolygon(
            shapelyPolygon=far_polygon,
            dataSourceOrData="lamas_population",
        )
        assert result.empty

    def test_invalid_datasource(self, demo_toolkit, population_gdf):
        """Non-existent datasource raises ValueError."""
        polygon = population_gdf.iloc[0].geometry
        with pytest.raises(ValueError):
            demo_toolkit.analysis.calculatePopulationInPolygon(
                shapelyPolygon=polygon,
                dataSourceOrData="non_existing_data_source",
            )

    def test_with_known_values(self, demo_toolkit):
        """Synthetic data yields correct area-weighted population."""
        geometry = [
            Polygon([(0, 0), (2, 0), (2, 2), (0, 2)]),   # area = 4, pop = 1000
            Polygon([(1, 1), (3, 1), (3, 3), (1, 3)]),   # area = 4, pop = 500
            Polygon([(5, 5), (6, 5), (6, 6), (5, 6)]),   # area = 1, pop = 200
        ]
        gdf = gpd.GeoDataFrame(
            {"total_pop": [1000, 500, 200], "geometry": geometry},
            crs="EPSG:4326",
        )

        test_poly = Polygon([(1, 1), (2.5, 1), (2.5, 2.5), (1, 2.5)])
        result = demo_toolkit.analysis.calculatePopulationInPolygon(
            shapelyPolygon=test_poly,
            dataSourceOrData=gdf,
            populationTypes="total_pop",
        )
        assert not result.empty
        total_estimated = result["total_pop"].sum()
        assert total_estimated > 0
        assert total_estimated < 1500

    def test_with_geodataframe_input(self, demo_toolkit, synthetic_gdf, small_polygon_itm):
        """Accepts GeoDataFrame directly instead of datasource name."""
        result = demo_toolkit.analysis.calculatePopulationInPolygon(
            shapelyPolygon=small_polygon_itm,
            dataSourceOrData=synthetic_gdf,
            populationTypes=["total_pop", "age_0_14"],
        )
        assert not result.empty
        assert "total_pop" in result.columns
        assert "age_0_14" in result.columns
        assert (result["areaFraction"] > 0).all()
        assert (result["areaFraction"] <= 1.0001).all()

    def test_specific_population_types(self, demo_toolkit, synthetic_gdf, large_polygon_itm):
        """Only requested population types are returned."""
        result = demo_toolkit.analysis.calculatePopulationInPolygon(
            shapelyPolygon=large_polygon_itm,
            dataSourceOrData=synthetic_gdf,
            populationTypes=["total_pop"],
        )
        assert "total_pop" in result.columns

    def test_full_coverage_sums_to_total(self, demo_toolkit, synthetic_gdf, large_polygon_itm):
        """Polygon covering all data sums to total population."""
        result = demo_toolkit.analysis.calculatePopulationInPolygon(
            shapelyPolygon=large_polygon_itm,
            dataSourceOrData=synthetic_gdf,
            populationTypes=["total_pop"],
        )
        total = result["total_pop"].sum()
        expected = synthetic_gdf["total_pop"].sum()
        assert abs(total - expected) < 1.0


# ===========================================================================
# 4. analysis.createNewArea
# ===========================================================================

@pytest.mark.integration
class TestCreateNewArea:
    def test_simple(self, demo_toolkit, population_gdf):
        """Create area covering all population data."""
        bounds = population_gdf.total_bounds
        buf = 10
        polygon = Polygon([
            (bounds[0] - buf, bounds[1] - buf),
            (bounds[2] + buf, bounds[1] - buf),
            (bounds[2] + buf, bounds[3] + buf),
            (bounds[0] - buf, bounds[3] + buf),
            (bounds[0] - buf, bounds[1] - buf),
        ])
        if not polygon.is_valid:
            polygon = polygon.buffer(0)

        result = demo_toolkit.analysis.createNewArea(
            shapeNameOrData=polygon,
            dataSourceOrData="lamas_population",
            saveMode=TOOLKIT_SAVEMODE_NOSAVE,
        )

        assert isinstance(result, nonDBMetadataFrame)
        gdf = result.getData()
        assert isinstance(gdf, gpd.GeoDataFrame)
        assert len(gdf) == 1
        assert "geometry" in gdf.columns
        assert "total_pop" in gdf.columns

        expected = population_gdf["total_pop"].sum()
        actual = gdf.iloc[0]["total_pop"]
        assert abs(actual - expected) < 0.01

    def test_with_synthetic_geodataframe(self, demo_toolkit, synthetic_gdf, large_polygon_itm):
        """createNewArea with GeoDataFrame input and synthetic data."""
        result = demo_toolkit.analysis.createNewArea(
            shapeNameOrData=large_polygon_itm,
            dataSourceOrData=synthetic_gdf,
            saveMode=TOOLKIT_SAVEMODE_NOSAVE,
        )
        gdf = result.getData()
        assert len(gdf) == 1
        total = gdf.iloc[0]["total_pop"]
        expected = synthetic_gdf["total_pop"].sum()
        assert abs(total - expected) < 1.0

    def test_partial_coverage(self, demo_toolkit, synthetic_gdf, small_polygon_itm):
        """Partial polygon yields less than total population."""
        result = demo_toolkit.analysis.createNewArea(
            shapeNameOrData=small_polygon_itm,
            dataSourceOrData=synthetic_gdf,
            saveMode=TOOLKIT_SAVEMODE_NOSAVE,
        )
        gdf = result.getData()
        partial_pop = gdf.iloc[0]["total_pop"]
        total_pop = synthetic_gdf["total_pop"].sum()
        assert partial_pop > 0
        assert partial_pop < total_pop

    def test_with_specific_population_types(self, demo_toolkit, synthetic_gdf, large_polygon_itm):
        """Only requested population types in output."""
        result = demo_toolkit.analysis.createNewArea(
            shapeNameOrData=large_polygon_itm,
            dataSourceOrData=synthetic_gdf,
            populationTypes=["total_pop", "age_0_14"],
            saveMode=TOOLKIT_SAVEMODE_NOSAVE,
        )
        gdf = result.getData()
        assert "total_pop" in gdf.columns
        assert "age_0_14" in gdf.columns


# ===========================================================================
# 5. presentation.plotPopulationDensity
# ===========================================================================

@pytest.mark.integration
class TestPlotPopulationDensity:
    """Test plotPopulationDensity with synthetic data (no DB needed)."""

    def test_basic(self, demo_toolkit, synthetic_gdf):
        ax = demo_toolkit.presentation.plotPopulationDensity(synthetic_gdf)
        assert ax is not None
        plt.close("all")

    def test_custom_density_units_dunam(self, demo_toolkit, synthetic_gdf):
        ax = demo_toolkit.presentation.plotPopulationDensity(
            synthetic_gdf, density_units=ureg.dunam,
        )
        assert ax is not None
        plt.close("all")

    def test_custom_density_units_km2(self, demo_toolkit, synthetic_gdf):
        ax = demo_toolkit.presentation.plotPopulationDensity(
            synthetic_gdf, density_units=ureg.km ** 2,
        )
        assert ax is not None
        plt.close("all")

    def test_custom_density_units_hectare(self, demo_toolkit, synthetic_gdf):
        ax = demo_toolkit.presentation.plotPopulationDensity(
            synthetic_gdf, density_units=ureg.hectare,
        )
        assert ax is not None
        plt.close("all")

    def test_specific_population_type(self, demo_toolkit, synthetic_gdf):
        ax = demo_toolkit.presentation.plotPopulationDensity(
            synthetic_gdf, populationType="age_0_14",
        )
        assert ax is not None
        plt.close("all")

    def test_custom_cmap_and_limits(self, demo_toolkit, synthetic_gdf):
        ax = demo_toolkit.presentation.plotPopulationDensity(
            synthetic_gdf, cmap="Blues", vmin=0, vmax=5,
            alpha=0.7, title="Density Test",
        )
        assert ax.get_title() == "Density Test"
        plt.close("all")

    def test_custom_figsize(self, demo_toolkit, synthetic_gdf):
        ax = demo_toolkit.presentation.plotPopulationDensity(
            synthetic_gdf, figsize=(8, 6),
        )
        fig = ax.get_figure()
        w, h = fig.get_size_inches()
        assert abs(w - 8) < 0.1
        assert abs(h - 6) < 0.1
        plt.close("all")

    def test_no_colorbar(self, demo_toolkit, synthetic_gdf):
        ax = demo_toolkit.presentation.plotPopulationDensity(
            synthetic_gdf, colorbar=False,
        )
        assert ax is not None
        plt.close("all")

    def test_with_axis_limits(self, demo_toolkit, synthetic_gdf):
        ax = demo_toolkit.presentation.plotPopulationDensity(
            synthetic_gdf,
            xlim=(170000, 172000), ylim=(660000, 662000),
        )
        assert ax.get_xlim()[0] == pytest.approx(170000, abs=1)
        assert ax.get_xlim()[1] == pytest.approx(172000, abs=1)
        plt.close("all")

    def test_with_existing_axes(self, demo_toolkit, synthetic_gdf):
        fig, ax_in = plt.subplots()
        ax_out = demo_toolkit.presentation.plotPopulationDensity(
            synthetic_gdf, ax=ax_in,
        )
        assert ax_out is ax_in
        plt.close("all")

    def test_wgs84_input_with_reprojection(self, demo_toolkit, synthetic_gdf_wgs84):
        ax = demo_toolkit.presentation.plotPopulationDensity(
            synthetic_gdf_wgs84, inputCRS=WSG84, outputCRS=ITM,
        )
        assert ax is not None
        plt.close("all")


# ===========================================================================
# 6. presentation.plotPopulation
# ===========================================================================

@pytest.mark.integration
class TestPlotPopulation:
    """Test plotPopulation (absolute counts)."""

    def test_basic(self, demo_toolkit, synthetic_gdf):
        ax = demo_toolkit.presentation.plotPopulation(synthetic_gdf)
        assert ax is not None
        plt.close("all")

    def test_specific_type(self, demo_toolkit, synthetic_gdf):
        ax = demo_toolkit.presentation.plotPopulation(
            synthetic_gdf, populationType="age_65_up",
        )
        assert ax is not None
        plt.close("all")

    def test_custom_styling(self, demo_toolkit, synthetic_gdf):
        ax = demo_toolkit.presentation.plotPopulation(
            synthetic_gdf,
            cmap="viridis", alpha=0.6, edgecolor="black", linewidth=1.0,
            title="Population Count",
        )
        assert ax.get_title() == "Population Count"
        plt.close("all")

    def test_with_vmin_vmax(self, demo_toolkit, synthetic_gdf):
        ax = demo_toolkit.presentation.plotPopulation(
            synthetic_gdf, vmin=0, vmax=3000,
        )
        assert ax is not None
        plt.close("all")

    def test_with_existing_axes(self, demo_toolkit, synthetic_gdf):
        fig, ax_in = plt.subplots()
        ax_out = demo_toolkit.presentation.plotPopulation(
            synthetic_gdf, ax=ax_in,
        )
        assert ax_out is ax_in
        plt.close("all")


# ===========================================================================
# 7. presentation.plotPopulationByType
# ===========================================================================

@pytest.mark.integration
class TestPlotPopulationByType:
    """Test plotPopulationByType (subplot grid)."""

    def test_basic_all_types(self, demo_toolkit, synthetic_gdf):
        fig = demo_toolkit.presentation.plotPopulationByType(synthetic_gdf)
        assert fig is not None
        plt.close("all")

    def test_specific_types(self, demo_toolkit, synthetic_gdf):
        fig = demo_toolkit.presentation.plotPopulationByType(
            synthetic_gdf, populationTypes=["total_pop", "age_0_14"],
        )
        assert fig is not None
        plt.close("all")

    def test_custom_ncols(self, demo_toolkit, synthetic_gdf):
        fig = demo_toolkit.presentation.plotPopulationByType(
            synthetic_gdf, ncols=2,
        )
        assert fig is not None
        plt.close("all")

    def test_single_type(self, demo_toolkit, synthetic_gdf):
        fig = demo_toolkit.presentation.plotPopulationByType(
            synthetic_gdf, populationTypes=["total_pop"], ncols=1,
        )
        assert fig is not None
        plt.close("all")

    def test_custom_figsize_and_cmap(self, demo_toolkit, synthetic_gdf):
        fig = demo_toolkit.presentation.plotPopulationByType(
            synthetic_gdf, figsize=(20, 12), cmap="Reds", alpha=0.6,
        )
        w, h = fig.get_size_inches()
        assert abs(w - 20) < 0.1
        plt.close("all")


# ===========================================================================
# 8. presentation.plotPopulationInPolygon
# ===========================================================================

@pytest.mark.integration
class TestPlotPopulationInPolygon:
    """Test plotPopulationInPolygon (intersection visualization)."""

    @pytest.fixture()
    def intersection_result(self, demo_toolkit, synthetic_gdf, small_polygon_itm):
        return demo_toolkit.analysis.calculatePopulationInPolygon(
            shapelyPolygon=small_polygon_itm,
            dataSourceOrData=synthetic_gdf,
        )

    def test_basic(self, demo_toolkit, intersection_result):
        ax = demo_toolkit.presentation.plotPopulationInPolygon(intersection_result)
        assert ax is not None
        plt.close("all")

    def test_with_query_polygon(self, demo_toolkit, intersection_result, small_polygon_itm):
        ax = demo_toolkit.presentation.plotPopulationInPolygon(
            intersection_result, queryPolygon=small_polygon_itm,
        )
        assert ax is not None
        plt.close("all")

    def test_with_context_data(self, demo_toolkit, intersection_result, synthetic_gdf):
        ax = demo_toolkit.presentation.plotPopulationInPolygon(
            intersection_result, contextData=synthetic_gdf,
        )
        assert ax is not None
        plt.close("all")

    def test_with_all_options(self, demo_toolkit, intersection_result,
                              small_polygon_itm, synthetic_gdf):
        ax = demo_toolkit.presentation.plotPopulationInPolygon(
            intersection_result,
            queryPolygon=small_polygon_itm,
            contextData=synthetic_gdf,
            populationType="total_pop",
            cmap="OrRd",
            alpha=0.8,
            title="Intersection Test",
        )
        assert ax.get_title() == "Intersection Test"
        plt.close("all")


# ===========================================================================
# 9. presentation.plotArea
# ===========================================================================

@pytest.mark.integration
class TestPlotArea:
    """Test plotArea (custom area with population annotation)."""

    @pytest.fixture()
    def area_result(self, demo_toolkit, synthetic_gdf, large_polygon_itm):
        result = demo_toolkit.analysis.createNewArea(
            shapeNameOrData=large_polygon_itm,
            dataSourceOrData=synthetic_gdf,
            saveMode=TOOLKIT_SAVEMODE_NOSAVE,
        )
        return result.getData()

    def test_basic(self, demo_toolkit, area_result):
        ax = demo_toolkit.presentation.plotArea(area_result)
        assert ax is not None
        plt.close("all")

    def test_with_context(self, demo_toolkit, area_result, synthetic_gdf):
        ax = demo_toolkit.presentation.plotArea(
            area_result, contextData=synthetic_gdf,
        )
        assert ax is not None
        plt.close("all")

    def test_with_annotation(self, demo_toolkit, area_result):
        ax = demo_toolkit.presentation.plotArea(
            area_result, annotate=True, annotate_fontsize=14,
        )
        assert ax is not None
        # Check annotation text exists
        texts = ax.texts
        assert len(texts) >= 1
        plt.close("all")

    def test_without_annotation(self, demo_toolkit, area_result):
        ax = demo_toolkit.presentation.plotArea(
            area_result, annotate=False,
        )
        assert ax is not None
        plt.close("all")

    def test_custom_styling(self, demo_toolkit, area_result, synthetic_gdf):
        ax = demo_toolkit.presentation.plotArea(
            area_result,
            contextData=synthetic_gdf,
            areaColor="red",
            contextColor="0.8",
            alpha=0.4,
            title="Custom Area",
        )
        assert ax.get_title() == "Custom Area"
        plt.close("all")


# ===========================================================================
# 10. presentation.plotPopulationOnMap (mocked TilesToolkit)
# ===========================================================================

@pytest.mark.integration
class TestPlotPopulationOnMap:
    """Test plotPopulationOnMap with a mocked TilesToolkit."""

    def _make_mock_tiles(self, synthetic_gdf):
        """Create a minimal mock for TilesToolkit."""

        class MockTilesPresentation:
            def plot(self, dataSourceName=None, extent=None, ax=None,
                     zoomlevel=14, tileServer=None, **kwargs):
                return ax

        class MockTilesToolkit:
            presentation = MockTilesPresentation()

            def getImageFromCorners(self, minx, miny, maxx, maxy, zoomlevel,
                                    tileServer=None, inputCRS=None, outputCRS=None):
                return None
            def getImageFromCorners(self, minx, miny, maxx, maxy,
                                    zoomlevel=14, tileServer=None,
                                    inputCRS=None, outputCRS=None):
                import numpy as np
                return np.zeros((256, 256, 3), dtype=np.uint8)

        return MockTilesToolkit()

    def test_density_mode(self, demo_toolkit, synthetic_gdf):
        mock_tiles = self._make_mock_tiles(synthetic_gdf)
        ax = demo_toolkit.presentation.plotPopulationOnMap(
            synthetic_gdf, tilesToolkit=mock_tiles, density=True,
        )
        assert ax is not None
        plt.close("all")

    def test_absolute_mode(self, demo_toolkit, synthetic_gdf):
        mock_tiles = self._make_mock_tiles(synthetic_gdf)
        ax = demo_toolkit.presentation.plotPopulationOnMap(
            synthetic_gdf, tilesToolkit=mock_tiles, density=False,
        )
        assert ax is not None
        plt.close("all")

    def test_custom_density_units(self, demo_toolkit, synthetic_gdf):
        mock_tiles = self._make_mock_tiles(synthetic_gdf)
        ax = demo_toolkit.presentation.plotPopulationOnMap(
            synthetic_gdf, tilesToolkit=mock_tiles,
            density=True, density_units=ureg.dunam,
        )
        assert ax is not None
        plt.close("all")

    def test_custom_zoom_and_styling(self, demo_toolkit, synthetic_gdf):
        mock_tiles = self._make_mock_tiles(synthetic_gdf)
        ax = demo_toolkit.presentation.plotPopulationOnMap(
            synthetic_gdf, tilesToolkit=mock_tiles,
            zoomlevel=12, cmap="Reds", alpha=0.4,
            title="Map Overlay",
        )
        assert ax.get_title() == "Map Overlay"
        plt.close("all")


# ===========================================================================
# 11. Edge cases
# ===========================================================================

@pytest.mark.integration
class TestEdgeCases:
    """Edge cases and error handling."""

    def test_empty_geodataframe(self, demo_toolkit):
        """Plotting empty GeoDataFrame should not crash."""
        empty = gpd.GeoDataFrame(
            {"total_pop": [], "geometry": []},
            crs="EPSG:2039",
        )
        # plotPopulation with empty data — should return ax without error
        with pytest.warns(UserWarning, match="empty"):
            ax = demo_toolkit.presentation.plotPopulation(empty)
        assert ax is not None
        plt.close("all")

    def test_single_polygon(self, demo_toolkit):
        """Single polygon data works correctly."""
        gdf = gpd.GeoDataFrame(
            {"total_pop": [5000], "geometry": [box(170000, 660000, 171000, 661000)]},
            crs="EPSG:2039",
        )
        ax = demo_toolkit.presentation.plotPopulationDensity(gdf)
        assert ax is not None
        plt.close("all")

    def test_intersection_with_no_overlap(self, demo_toolkit, synthetic_gdf):
        """calculatePopulationInPolygon with non-overlapping polygon."""
        far_polygon = box(200000, 700000, 201000, 701000)
        result = demo_toolkit.analysis.calculatePopulationInPolygon(
            shapelyPolygon=far_polygon,
            dataSourceOrData=synthetic_gdf,
        )
        assert result.empty
