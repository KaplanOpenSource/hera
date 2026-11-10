"""A shapes file from disk as FlatGeobuf, for the map preview."""
from pathlib import Path

from fastapi import APIRouter, HTTPException, Response

from dependencies import warmup
from geo_shapes import read_geo_shapes

router = APIRouter()

FLATGEOBUF_TYPE = "application/octet-stream"


@router.get("/geoshapes")
def geoshapes(path: str) -> Response:
    """The file at `path` as WGS84 FlatGeobuf. See ``geo_shapes.read_geo_shapes``."""
    if not warmup.ready:
        raise HTTPException(status_code=503, detail="WARMING_UP")
    if not Path(path).is_file():
        raise HTTPException(status_code=404, detail=f"not a file: {path}")
    return Response(content=read_geo_shapes(path), media_type=FLATGEOBUF_TYPE)
