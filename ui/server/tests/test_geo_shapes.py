"""Reading shapes files as WGS84 FlatGeobuf, including the ITM files that declare no CRS."""
import io
import json
import sys
import types

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import LineString, MultiPolygon, Point, Polygon

import dependencies
import geo_shapes
from fastapi import HTTPException
from routers import geoshapes

ITM = 2039
WGS84 = 4326
# Under a millimetre, in degrees.
TOLERANCE = 1e-8

# Tel Aviv area in ITM metres.
X, Y = 178000.0, 663000.0

SQUARE = Polygon([(X, Y), (X + 100, Y), (X + 100, Y + 100), (X, Y + 100)])
SQUARE_WITH_HOLE = Polygon(
    [(X, Y), (X + 100, Y), (X + 100, Y + 100), (X, Y + 100)],
    [[(X + 40, Y + 40), (X + 60, Y + 40), (X + 60, Y + 60), (X + 40, Y + 60)]],
)
TWO_SQUARES = MultiPolygon([SQUARE, Polygon([(X + 200, Y), (X + 300, Y), (X + 300, Y + 100), (X + 200, Y + 100)])])

SHAPES = {
    "points": [Point(X, Y), Point(X + 50, Y + 70)],
    "lines": [LineString([(X, Y), (X + 10, Y + 20), (X + 30, Y + 5)]), LineString([(X, Y + 1), (X + 1, Y + 2)])],
    "polygons": [SQUARE, SQUARE_WITH_HOLE],
    "multipolygons": [TWO_SQUARES, MultiPolygon([SQUARE_WITH_HOLE])],
}

EXTENSIONS = {"shp": "ESRI Shapefile", "geojson": "GeoJSON", "fgb": "FlatGeobuf"}


@pytest.fixture(autouse=True)
def fake_hera_gis(monkeypatch):
    """The two EPSG codes the reader asks hera for, without installing hera."""
    utils = types.ModuleType("hera.measurements.GIS.utils")
    utils.ITM = ITM
    utils.WGS84 = WGS84
    monkeypatch.setitem(sys.modules, "hera", types.ModuleType("hera"))
    monkeypatch.setitem(sys.modules, "hera.measurements", types.ModuleType("hera.measurements"))
    monkeypatch.setitem(sys.modules, "hera.measurements.GIS", types.ModuleType("hera.measurements.GIS"))
    monkeypatch.setitem(sys.modules, "hera.measurements.GIS.utils", utils)


@pytest.fixture
def warmed(monkeypatch):
    monkeypatch.setattr(dependencies.warmup, "_ready", True)


def _frame(geometries, crs):
    return gpd.GeoDataFrame(
        {"name": [f"s{i}" for i in range(len(geometries))], "height": [10.5 * (i + 1) for i in range(len(geometries))]},
        geometry=geometries,
        crs=crs,
    )


def _write(tmp_path, gdf, extension):
    path = tmp_path / f"shapes.{extension}"
    options = {}
    if extension == "fgb":
        # An index would sort the shapes, so the file order would not be ours.
        options["SPATIAL_INDEX"] = "NO"
    gdf.to_file(path, driver=EXTENSIONS[extension], engine="pyogrio", **options)
    return str(path)


def _decode(data):
    return gpd.read_file(io.BytesIO(data))


def _assert_same_shapes(result, expected):
    assert result.crs.to_epsg() == WGS84
    assert len(result) == len(expected)
    pd.testing.assert_frame_equal(
        pd.DataFrame(result.drop(columns="geometry")),
        pd.DataFrame(expected.drop(columns="geometry")),
        check_dtype=False,
    )
    for got, want in zip(result.geometry, expected.geometry):
        assert got.geom_type == want.geom_type
        # Shapefiles may flip ring order, so compare in a fixed vertex order.
        assert got.normalize().equals_exact(want.normalize(), TOLERANCE)


@pytest.mark.parametrize("extension", EXTENSIONS)
@pytest.mark.parametrize("kind", SHAPES)
def test_exact_shapes_from_itm(tmp_path, extension, kind):
    source = _frame(SHAPES[kind], ITM)
    path = _write(tmp_path, source, extension)

    result = _decode(geo_shapes.read_geo_shapes(path))

    _assert_same_shapes(result, source.to_crs(epsg=WGS84))


@pytest.mark.parametrize("extension", EXTENSIONS)
def test_exact_shapes_from_wgs84(tmp_path, extension):
    source = _frame(SHAPES["polygons"], ITM).to_crs(epsg=WGS84)
    path = _write(tmp_path, source, extension)

    result = _decode(geo_shapes.read_geo_shapes(path))

    _assert_same_shapes(result, source)


def test_geojson_metres_with_no_crs_are_read_as_itm(tmp_path):
    doc = {
        "type": "FeatureCollection",
        "features": [{"type": "Feature", "properties": {"name": "p"}, "geometry": {"type": "Point", "coordinates": [X, Y]}}],
    }
    path = tmp_path / "itm.geojson"
    path.write_text(json.dumps(doc))

    result = _decode(geo_shapes.read_geo_shapes(str(path)))

    expected = gpd.GeoSeries([Point(X, Y)], crs=ITM).to_crs(epsg=WGS84)
    assert result.geometry[0].equals_exact(expected[0], TOLERANCE)


def test_shp_with_no_prj_is_read_as_itm(tmp_path):
    source = _frame(SHAPES["polygons"], ITM)
    path = _write(tmp_path, source, "shp")
    (tmp_path / "shapes.prj").unlink()

    result = _decode(geo_shapes.read_geo_shapes(path))

    _assert_same_shapes(result, source.to_crs(epsg=WGS84))


def test_degrees_with_no_crs_are_left_alone(tmp_path):
    source = _frame([Point(35.2, 32.1)], WGS84)
    path = _write(tmp_path, source, "shp")
    (tmp_path / "shapes.prj").unlink()

    result = _decode(geo_shapes.read_geo_shapes(path))

    _assert_same_shapes(result, source)


def test_the_route_sends_flatgeobuf(tmp_path, warmed):
    source = _frame(SHAPES["points"], ITM)
    path = _write(tmp_path, source, "shp")

    response = geoshapes.geoshapes(path)

    assert response.media_type == geoshapes.FLATGEOBUF_TYPE
    _assert_same_shapes(_decode(response.body), source.to_crs(epsg=WGS84))


def test_the_route_reports_a_missing_file(tmp_path, warmed):
    with pytest.raises(HTTPException) as raised:
        geoshapes.geoshapes(str(tmp_path / "nope.shp"))

    assert raised.value.status_code == 404


def test_the_route_waits_for_warmup(tmp_path, monkeypatch):
    monkeypatch.setattr(dependencies.warmup, "_ready", False)
    path = _write(tmp_path, _frame(SHAPES["points"], ITM), "shp")

    with pytest.raises(HTTPException) as raised:
        geoshapes.geoshapes(path)

    assert raised.value.status_code == 503
