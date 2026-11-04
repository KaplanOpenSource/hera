"""BuildingsToolkit's two pure static helpers: height extraction and a
bounding-box spatial filter over GeoJSON-shaped feature dicts.

The rest of BuildingsToolkit (getBuildingsFromRectangle,
getBuildingHeightFromRasterTopographyToolkit, and the FreeCAD-based STL
export) needs a real registered datasource or a real FreeCAD install --
left for an integration test.
"""
import geopandas
import pytest
from shapely.geometry import Point, Polygon, mapping

from hera.measurements.GIS.vector.buildings.toolkit import BuildingsToolkit


@pytest.mark.unit
class TestGetBuildingsHeight:
    def test_an_explicit_height_is_kept_as_is(self):
        gdf = geopandas.GeoDataFrame([
            {"name": "A", "geometry": Point(0, 0), "height": 12.0},
        ])
        result = BuildingsToolkit.get_buildings_height(gdf)
        assert result.iloc[0]["height"] == 12.0

    def test_missing_height_is_estimated_from_levels_times_three(self):
        gdf = geopandas.GeoDataFrame([
            {"name": "B", "geometry": Point(0, 0), "building:levels": 4},
        ])
        result = BuildingsToolkit.get_buildings_height(gdf)
        assert result.iloc[0]["height"] == 12

    def test_neither_height_nor_levels_gives_none(self):
        gdf = geopandas.GeoDataFrame([{"name": "C", "geometry": Point(0, 0)}])
        result = BuildingsToolkit.get_buildings_height(gdf)
        assert result.iloc[0]["height"] is None

    def test_a_missing_name_defaults_to_unnamed(self):
        gdf = geopandas.GeoDataFrame([{"geometry": Point(0, 0), "height": 5.0}])
        result = BuildingsToolkit.get_buildings_height(gdf)
        assert result.iloc[0]["name"] == "Unnamed"

    def test_the_geometry_is_preserved_unchanged(self):
        point = Point(1, 2)
        gdf = geopandas.GeoDataFrame([{"geometry": point, "height": 5.0}])
        result = BuildingsToolkit.get_buildings_height(gdf)
        assert result.iloc[0]["geometry"] == point


@pytest.mark.unit
class TestFilterBuildingsInArea:
    def _feature(self, poly, **properties):
        return {"type": "Feature", "geometry": mapping(poly), "properties": dict(properties)}

    def test_a_building_inside_the_box_is_kept(self):
        data = {"features": [self._feature(Polygon([(1, 1), (2, 1), (2, 2), (1, 2)]), id=1)]}
        result = BuildingsToolkit.filter_buildings_in_area(data, 0, 0, 5, 5)
        assert len(result) == 1
        assert result.iloc[0]["id"] == 1

    def test_a_building_entirely_outside_the_box_is_dropped(self):
        data = {"features": [self._feature(Polygon([(100, 100), (101, 100), (101, 101), (100, 101)]), id=2)]}
        result = BuildingsToolkit.filter_buildings_in_area(data, 0, 0, 5, 5)
        assert len(result) == 0

    def test_a_building_straddling_the_boundary_is_kept(self):
        data = {"features": [self._feature(Polygon([(4, 4), (6, 4), (6, 6), (4, 6)]), id=3)]}
        result = BuildingsToolkit.filter_buildings_in_area(data, 0, 0, 5, 5)
        assert len(result) == 1

    def test_the_result_is_a_geodataframe_with_the_geometry_attached(self):
        data = {"features": [self._feature(Polygon([(1, 1), (2, 1), (2, 2), (1, 2)]), id=1)]}
        result = BuildingsToolkit.filter_buildings_in_area(data, 0, 0, 5, 5)
        assert isinstance(result, geopandas.GeoDataFrame)
        assert result.iloc[0]["geometry"].is_valid
