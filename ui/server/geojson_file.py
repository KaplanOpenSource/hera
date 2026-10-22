"""Reading a GeoJSON file off disk as WGS84, the only projection Leaflet draws."""
from __future__ import annotations

from typing import Any, Dict

# Coordinates above this cannot be degrees, so the file is in metres.
MAX_DEGREES = 360


def read_geojson(path: str) -> Dict[str, Any]:
    """The file as GeoJSON in WGS84.

    GeoJSON is specified to be WGS84, but hera's files are sometimes written in
    ITM, with no CRS of their own. Large coordinates give that away, so they are
    read as ITM and converted.
    """
    import json

    import geopandas as gpd
    from hera.measurements.GIS.utils import ITM, WGS84

    gdf = gpd.read_file(path)
    if gdf.crs and not gdf.crs.is_geographic:
        gdf = gdf.to_crs(epsg=WGS84)
    else:
        minx, miny = gdf.total_bounds[0], gdf.total_bounds[1]
        if abs(minx) > MAX_DEGREES or abs(miny) > MAX_DEGREES:
            gdf = gdf.set_crs(epsg=ITM, allow_override=True)
            gdf = gdf.to_crs(epsg=WGS84)
    return json.loads(gdf.to_json())
