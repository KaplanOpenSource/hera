"""A shapes file from disk as FlatGeobuf, for the map preview."""
from pathlib import Path
from typing import Optional

from fastapi import APIRouter, HTTPException, Response

from dependencies import warmup
from geo_shapes import Bbox, read_geo_shapes

router = APIRouter()

FLATGEOBUF_TYPE = "application/octet-stream"


def parse_bbox(text: str) -> Bbox:
    """`minx,miny,maxx,maxy` in WGS84 degrees."""
    parts = text.split(",")
    if len(parts) != 4:
        raise HTTPException(status_code=400, detail=f"bbox needs 4 numbers: {text}")
    try:
        minx, miny, maxx, maxy = (float(part) for part in parts)
    except ValueError:
        # A bad number is the caller's mistake, so it is a 400, not a server error.
        raise HTTPException(status_code=400, detail=f"bbox needs 4 numbers: {text}")
    if minx > maxx or miny > maxy:
        raise HTTPException(status_code=400, detail=f"bbox min is over max: {text}")
    return minx, miny, maxx, maxy


@router.get("/geoshapes")
def geoshapes(path: str, bbox: Optional[str] = None) -> Response:
    """The file at `path` as WGS84 FlatGeobuf, all of it or only `bbox`. See ``geo_shapes.read_geo_shapes``."""
    if not warmup.ready:
        raise HTTPException(status_code=503, detail="WARMING_UP")
    if not Path(path).is_file():
        raise HTTPException(status_code=404, detail=f"not a file: {path}")
    box = None
    if bbox is not None:
        box = parse_bbox(bbox)
    return Response(content=read_geo_shapes(path, box), media_type=FLATGEOBUF_TYPE)
