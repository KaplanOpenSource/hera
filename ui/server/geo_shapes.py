"""Reading a shapes file off disk as WGS84 FlatGeobuf, for the map preview."""
from __future__ import annotations

import io

# Coordinates above this cannot be degrees, so the file is in metres.
MAX_DEGREES = 360


def read_geo_shapes(path: str) -> bytes:
    """The whole file (.shp, .geojson, .fgb) as FlatGeobuf bytes in WGS84.

    hera's files are sometimes written in ITM with no CRS of their own. Large
    coordinates give that away, so they are read as ITM and converted.
    """
    import geopandas as gpd
    from hera.measurements.GIS.utils import ITM, WGS84

    gdf = gpd.read_file(path)
    if gdf.crs and not gdf.crs.is_geographic:
        gdf = gdf.to_crs(epsg=WGS84)
    else:
        # GeoJSON reads as WGS84 even when the file holds ITM metres.
        minx, miny = gdf.total_bounds[0], gdf.total_bounds[1]
        if abs(minx) > MAX_DEGREES or abs(miny) > MAX_DEGREES:
            gdf = gdf.set_crs(epsg=ITM, allow_override=True).to_crs(epsg=WGS84)
        elif gdf.crs is None:
            gdf = gdf.set_crs(epsg=WGS84)
    buffer = io.BytesIO()
    # No spatial index: it would sort the shapes, and phase 1 does not read by box.
    gdf.to_file(buffer, driver="FlatGeobuf", engine="pyogrio", SPATIAL_INDEX="NO")
    return buffer.getvalue()
