"""Reading GeoJSON files as WGS84, including the ITM files that declare no CRS."""
import json

import pytest

import dependencies
import geojson_file
from fastapi import HTTPException
from routers import geojson


def _write(path, coordinates, crs=None):
    doc = {
        "type": "FeatureCollection",
        "features": [{
            "type": "Feature",
            "properties": {"name": "p"},
            "geometry": {"type": "Point", "coordinates": coordinates},
        }],
    }
    if crs is not None:
        doc["crs"] = {"type": "name", "properties": {"name": crs}}
    path.write_text(json.dumps(doc))
    return str(path)


@pytest.fixture
def warmed(monkeypatch):
    monkeypatch.setattr(dependencies.warmup, "_ready", True)


def test_degrees_are_left_alone(tmp_path):
    path = _write(tmp_path / "wgs84.geojson", [35.2, 32.1])

    result = geojson_file.read_geojson(path)

    point = result["features"][0]["geometry"]["coordinates"]
    assert point == pytest.approx([35.2, 32.1])


def test_metres_with_no_crs_are_read_as_itm(tmp_path):
    # Tel Aviv in ITM; far outside the range of degrees.
    path = _write(tmp_path / "itm.geojson", [178000.0, 663000.0])

    result = geojson_file.read_geojson(path)

    longitude, latitude = result["features"][0]["geometry"]["coordinates"]
    assert 34 < longitude < 36
    assert 31 < latitude < 33


def test_the_route_reads_the_file(tmp_path, warmed):
    path = _write(tmp_path / "wgs84.geojson", [35.2, 32.1])

    result = geojson.geojson(path)

    assert result["features"][0]["properties"]["name"] == "p"


def test_the_route_reports_a_missing_file(tmp_path, warmed):
    with pytest.raises(HTTPException) as raised:
        geojson.geojson(str(tmp_path / "nope.geojson"))

    assert raised.value.status_code == 404


def test_the_route_waits_for_warmup(tmp_path, monkeypatch):
    monkeypatch.setattr(dependencies.warmup, "_ready", False)
    path = _write(tmp_path / "wgs84.geojson", [35.2, 32.1])

    with pytest.raises(HTTPException) as raised:
        geojson.geojson(path)

    assert raised.value.status_code == 503
