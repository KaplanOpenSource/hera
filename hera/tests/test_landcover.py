"""
Native pytest tests for the LandCoverToolkit.

Data is loaded into the test project via ``conftest.hera_test_project``.
The ``lc_toolkit`` fixture (session-scoped, from conftest) provides
a real ``LandCoverToolkit`` backed by MongoDB – no monkey-patching.

Replaces:
  - hera/measurements/GIS/raster/test_unit_test_gis_raster_landcover.py
  - hera/tests/json_definitions/test_definitions_landcover.json
"""

import os
import numpy as np
import pytest
import rasterio

_LC_TIF = os.path.join(
    os.environ.get("TEST_HERA", ""),
    "measurements", "GIS", "raster", "lc_mcd12q1.tif",
)
if not os.path.isfile(_LC_TIF):
    pytest.skip(
        f"Skipping landcover tests — data file not found: {_LC_TIF}",
        allow_module_level=True,
    )


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def landcover_file(lc_toolkit):
    """Get the landcover file path from the project's datasource.

    Since the datasource uses dataFormat="string", ``getDataSourceData``
    returns the absolute file path directly.
    """
    path = lc_toolkit.getDataSourceData("lc_mcd12q1")
    if path is None or not isinstance(path, str):
        pytest.skip("lc_mcd12q1 datasource not loaded or not a string path")
    return path


@pytest.fixture(scope="module")
def test_coords():
    """Standard test coordinates."""
    return {
        "lon": 35.0,
        "lat": 32.0,
        "minx": 34.9,
        "miny": 31.9,
        "maxx": 35.1,
        "maxy": 32.1,
    }


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestGetLandCoverAtPoint:
    def test_basic(self, lc_toolkit, test_coords):
        value = lc_toolkit.getLandCoverAtPoint(test_coords["lon"], test_coords["lat"])
        assert value is not None

    def test_against_raster(self, lc_toolkit, landcover_file, test_coords):
        lon, lat = test_coords["lon"], test_coords["lat"]
        with rasterio.open(landcover_file) as src:
            row, col = src.index(lon, lat)
            value_from_raster = src.read(1)[row, col]

        value_from_toolkit = lc_toolkit.getLandCoverAtPoint(lon, lat)
        assert value_from_raster == value_from_toolkit


class TestGetLandCover:
    def test_basic(self, lc_toolkit, test_coords):
        lc = lc_toolkit.getLandCover(
            test_coords["minx"], test_coords["miny"],
            test_coords["maxx"], test_coords["maxy"],
            dxdy=500,
        )
        assert lc is not None
        assert "landcover" in lc.coords

    def test_map_vs_raster(self, lc_toolkit, landcover_file, test_coords):
        lc = lc_toolkit.getLandCover(
            test_coords["minx"], test_coords["miny"],
            test_coords["maxx"], test_coords["maxy"],
            dxdy=500,
        )
        lc_values = lc["landcover"].values
        lons = lc["lon"].values
        lats = lc["lat"].values

        with rasterio.open(landcover_file) as src:
            raster_band = src.read(1)
            for i in range(0, lc_values.shape[0], 10):
                for j in range(0, lc_values.shape[1], 10):
                    lon, lat = lons[i, j], lats[i, j]
                    try:
                        row, col = src.index(lon, lat)
                        assert raster_band[row, col] == lc_values[i, j], (
                            f"Mismatch at (lon={lon}, lat={lat})"
                        )
                    except IndexError:
                        pass  # Point outside raster extent


class TestGetRoughness:
    def test_at_point(self, lc_toolkit, test_coords):
        value = lc_toolkit.getRoughnessAtPoint(test_coords["lon"], test_coords["lat"])
        assert value is not None

    def test_area(self, lc_toolkit, test_coords):
        roughness = lc_toolkit.getRoughness(
            test_coords["minx"], test_coords["miny"],
            test_coords["maxx"], test_coords["maxy"],
            dxdy=500,
        )
        assert roughness is not None
        assert "z0" in roughness.coords

    def test_values_in_range(self, lc_toolkit, test_coords):
        roughness = lc_toolkit.getRoughness(
            test_coords["minx"], test_coords["miny"],
            test_coords["maxx"], test_coords["maxy"],
            dxdy=500,
        )
        z0 = roughness["z0"].values
        assert np.all((z0 > 0) & (z0 < 2)), "Some roughness values out of expected range"


class TestRoughnessLength:
    def test_roughnesslength2sandgrainroughness(self, lc_toolkit):
        ks = lc_toolkit.roughnesslength2sandgrainroughness(0.1)
        assert ks > 0

    def test_known_landcover(self, lc_toolkit, test_coords):
        known_lc_value = 5  # Mixed forests
        expected_roughness = 1.0

        original = lc_toolkit.getLandCoverAtPoint
        lc_toolkit.getLandCoverAtPoint = lambda lon, lat, inputCRS=None, dataSourceName=None: known_lc_value
        try:
            roughness = lc_toolkit.getRoughnessAtPoint(test_coords["lon"], test_coords["lat"])
            assert abs(roughness - expected_roughness) < 1e-3
        finally:
            lc_toolkit.getLandCoverAtPoint = original


class TestMiscLandCover:
    def test_out_of_bounds(self, lc_toolkit):
        with pytest.raises(IndexError):
            lc_toolkit.getLandCoverAtPoint(1000, 1000)

    def test_get_coding_map(self, lc_toolkit):
        coding_map = lc_toolkit.getCodingMap("Type-1")
        assert isinstance(coding_map, dict)
        assert 0 in coding_map
        assert coding_map[0] == "Water"
