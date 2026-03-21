"""
Native pytest tests for the DemographyToolkit.

Data is loaded into the test project via ``conftest.hera_test_project``.
The ``demo_toolkit`` fixture (session-scoped, from conftest) provides
a real ``DemographyToolkit`` backed by MongoDB – no monkey-patching
or direct file-path access.

Replaces:
  - hera/measurements/GIS/vector/test_unit_demography_toolkit.py
  - hera/tests/json_definitions/demography_test_definitions.json
"""

import os
import tempfile

import geopandas as gpd
import pytest
from shapely.geometry import Polygon

from hera.toolkit import TOOLKIT_SAVEMODE_NOSAVE
from hera.datalayer.document.metadataDocument import nonDBMetadataFrame


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def population_gdf(demo_toolkit):
    """Load the population GeoDataFrame through the project's datasource."""
    gdf = demo_toolkit.getDataSourceData("lamas_population")
    if gdf is None:
        pytest.skip("lamas_population datasource not loaded in project")
    return gdf


# ---------------------------------------------------------------------------
# calculatePopulationInPolygon tests
# ---------------------------------------------------------------------------

class TestCalculatePopulationInPolygon:
    def test_basic(self, demo_toolkit, population_gdf):
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
        assert not result.empty, "Expected non-empty result"
        assert "geometry" in result.columns
        assert "areaFraction" in result.columns

    def test_partial_intersection(self, demo_toolkit, population_gdf):
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
        assert result.empty, "Expected empty result for polygon outside data bounds"

    def test_invalid_datasource(self, demo_toolkit, population_gdf):
        polygon = population_gdf.iloc[0].geometry
        with pytest.raises(ValueError):
            demo_toolkit.analysis.calculatePopulationInPolygon(
                shapelyPolygon=polygon,
                dataSourceOrData="non_existing_data_source",
            )

    def test_with_known_values(self, demo_toolkit):
        geometry = [
            Polygon([(0, 0), (2, 0), (2, 2), (0, 2)]),
            Polygon([(1, 1), (3, 1), (3, 3), (1, 3)]),
            Polygon([(5, 5), (6, 5), (6, 6), (5, 6)]),
        ]
        total_pop = [1000, 500, 200]
        gdf = gpd.GeoDataFrame(
            {"total_pop": total_pop, "geometry": geometry}, crs="EPSG:4326"
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


# ---------------------------------------------------------------------------
# createNewArea tests
# ---------------------------------------------------------------------------

class TestCreateNewArea:
    def test_simple(self, demo_toolkit, population_gdf):
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


# ---------------------------------------------------------------------------
# setDefaultDirectory tests
# ---------------------------------------------------------------------------

class TestSetDefaultDirectory:
    def test_creates_and_sets_path(self, demo_toolkit):
        with tempfile.TemporaryDirectory() as tmpdir:
            test_dir = os.path.join(tmpdir, "demography_test_folder")
            demo_toolkit.setDefaultDirectory(test_dir, create=True)
            assert os.path.exists(test_dir)
            assert os.path.abspath(test_dir) == demo_toolkit.filesDirectory
