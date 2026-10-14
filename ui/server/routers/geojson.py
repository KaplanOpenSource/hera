"""A GeoJSON file from disk, for the map preview."""
from pathlib import Path

from fastapi import APIRouter, HTTPException

from dependencies import warmup
from geojson_file import read_geojson

router = APIRouter()


@router.get("/geojson")
def geojson(path: str) -> dict:
    """The file at `path` as GeoJSON in WGS84. See ``geojson_file.read_geojson``."""
    if not warmup.ready:
        raise HTTPException(status_code=503, detail="WARMING_UP")
    if not Path(path).is_file():
        raise HTTPException(status_code=404, detail=f"not a file: {path}")
    return read_geojson(path)
