"""TilesToolkit: slippy-map tile math and configuration.

``listImages``/``getImageFromCorners`` need a real tile server or stored
imagery and are left for integration tests; ``presentation.plot`` needs
real image data. What's covered here is the pure Web Mercator tile math
and the toolkit's config/doctype plumbing.
"""
import pytest

from hera import toolkitHome


@pytest.fixture()
def tiles(unit_toolkit_factory):
    return unit_toolkit_factory(toolkitHome.GIS_TILES)


@pytest.mark.unit
class TestDoctype:
    def test_doctype_is_the_toolkit_name_plus_png(self, tiles):
        assert tiles.doctype == f"{tiles.toolkitName}_PNG"


@pytest.mark.unit
class TestTileScaleAtLatLonZoom:
    def test_scale_halves_with_each_zoom_level(self, tiles):
        scale_z1 = tiles.tileScaleAtLatLonZoom(latitude=0, longitude=0, zoomlevel=1)
        scale_z2 = tiles.tileScaleAtLatLonZoom(latitude=0, longitude=0, zoomlevel=2)
        assert scale_z1 / scale_z2 == pytest.approx(2.0)

    def test_scale_shrinks_towards_the_poles(self, tiles):
        equator = tiles.tileScaleAtLatLonZoom(latitude=0, longitude=0, zoomlevel=5)
        mid_lat = tiles.tileScaleAtLatLonZoom(latitude=60, longitude=0, zoomlevel=5)
        assert mid_lat < equator


@pytest.mark.unit
class TestTileDegreeConversion:
    def test_deg2tile_and_tile2deg_are_approximate_inverses(self, tiles):
        xtile, ytile = tiles.deg2tile(lat_deg=32.08, lon_deg=34.78, zoom=10)
        lat, lon = tiles.tile2deg(xtile, ytile, zoom=10)
        # a single tile at zoom 10 spans a few tenths of a degree
        assert abs(lat - 32.08) < 1.0
        assert abs(lon - 34.78) < 1.0

    def test_the_origin_maps_to_the_top_left_tile_at_zoom_zero(self, tiles):
        assert tiles.deg2tile(lat_deg=0.0, lon_deg=-180.0, zoom=0) == (0, 0)

    def test_tile2deg_at_zoom_zero_covers_the_whole_map(self, tiles):
        lat, lon = tiles.tile2deg(xtile=0, ytile=0, zoom=0)
        assert lon == pytest.approx(-180.0)
        assert lat > 0


@pytest.mark.unit
class TestSetDefaultTileServer:
    def test_it_is_stored_in_the_project_config(self, tiles):
        tiles.setDefaultTileServer("https://example.com/{z}/{x}/{y}.png")
        assert tiles.getConfig()["defaultTileServer"] == "https://example.com/{z}/{x}/{y}.png"
