"""GIS/vector/_io_utils.py::readGeoJSONString and GIS/utils.py::create_xarray."""
import json

import pytest

from hera.measurements.GIS.utils import create_xarray
from hera.measurements.GIS.vector._io_utils import readGeoJSONString


@pytest.mark.unit
class TestReadGeoJSONString:
    def test_a_valid_geojson_string_parses_to_a_geodataframe(self):
        gj = json.dumps({
            "type": "FeatureCollection",
            "features": [{"type": "Feature", "properties": {"name": "A"},
                          "geometry": {"type": "Point", "coordinates": [0, 0]}}],
        })
        result = readGeoJSONString(gj)
        assert list(result["name"]) == ["A"]

    def test_a_non_string_raises(self):
        with pytest.raises(ValueError, match="Expected a GeoJSON string"):
            readGeoJSONString(123)

    def test_content_that_is_neither_a_path_nor_geojson_raises(self):
        with pytest.raises(ValueError, match="not a path on disk"):
            readGeoJSONString("not json")


@pytest.mark.unit
class TestCreateXarray:
    def test_it_builds_a_2d_grid_over_the_bounding_box(self):
        result = create_xarray(34.0, 32.0, 34.1, 32.1, dxdy=1000)
        assert result.dims == ("i", "j")
        assert result.ndim == 2
