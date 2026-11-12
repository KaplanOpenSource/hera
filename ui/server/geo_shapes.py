"""Reading a shapes file off disk as WGS84 FlatGeobuf, for the map preview."""
from __future__ import annotations

import io
from typing import Optional, Tuple

# Coordinates above this cannot be degrees, so the file is in metres.
MAX_DEGREES = 360
# Above this many shapes, send one dot per shape instead of the full shapes.
MAX_SHAPES = 10_000
# Points along each bbox edge when changing CRS, so curved edges stay inside.
BBOX_DENSIFY_POINTS = 21

Bbox = Tuple[float, float, float, float]


def dots_title() -> str:
    """The FlatGeobuf title that marks dots; the client reads the limit from it."""
    return f"dots max={MAX_SHAPES}"


def _file_crs(path: str):
    """The CRS to read the file in: its own, or ITM when it holds metres but says degrees or nothing.

    GeoJSON reads as WGS84 even when the file holds ITM metres.
    """
    import pyogrio
    from pyproj import CRS
    from hera.measurements.GIS.utils import ITM, WGS84

    info = pyogrio.read_info(path, force_total_bounds=True)
    if info["crs"] and not CRS(info["crs"]).is_geographic:
        return CRS(info["crs"])
    minx, miny = info["total_bounds"][0], info["total_bounds"][1]
    if abs(minx) > MAX_DEGREES or abs(miny) > MAX_DEGREES:
        return CRS.from_epsg(ITM)
    if info["crs"]:
        return CRS(info["crs"])
    return CRS.from_epsg(WGS84)


def read_geo_shapes(path: str, bbox: Optional[Bbox] = None) -> bytes:
    """The file (.shp, .geojson, .fgb) as FlatGeobuf bytes in WGS84.

    With no bbox, the whole file. Over MAX_SHAPES shapes, it sends one dot per
    shape at its center, with no columns, and marks it with dots_title().

    With a bbox (WGS84 minx, miny, maxx, maxy), the full shapes in it only.

    hera's files are sometimes written in ITM with no CRS of their own. Large
    coordinates give that away, so they are read as ITM and converted.

    No shapes give empty bytes: GDAL cannot write an empty FlatGeobuf.
    """
    import geopandas as gpd
    import pyogrio
    import shapely
    from pyproj import Transformer
    from hera.measurements.GIS.utils import WGS84

    crs = _file_crs(path)
    options = {}
    if bbox is not None:
        to_file_crs = Transformer.from_crs(WGS84, crs, always_xy=True)
        gdf = gpd.read_file(path, bbox=to_file_crs.transform_bounds(*bbox, densify_pts=BBOX_DENSIFY_POINTS))
    elif pyogrio.read_info(path, force_feature_count=True)["features"] > MAX_SHAPES:
        shapes = gpd.read_file(path, columns=[])
        # Centers in the file's own CRS, before any change to degrees.
        gdf = gpd.GeoDataFrame(geometry=shapely.centroid(shapes.geometry.values))
        options["TITLE"] = dots_title()
    else:
        gdf = gpd.read_file(path)
    if gdf.empty:
        return b""
    gdf = gdf.set_crs(crs, allow_override=True).to_crs(epsg=WGS84)
    buffer = io.BytesIO()
    # No spatial index: it would sort the shapes, and the client does not read by box.
    gdf.to_file(buffer, driver="FlatGeobuf", engine="pyogrio", SPATIAL_INDEX="NO", **options)
    return buffer.getvalue()
