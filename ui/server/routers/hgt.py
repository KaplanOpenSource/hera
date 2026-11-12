"""An SRTM .hgt file from disk as a GeoTIFF, for the TIFF preview."""
from pathlib import Path

from fastapi import APIRouter, HTTPException, Response

from dependencies import warmup
from hgt_file import read_hgt

router = APIRouter()

GEOTIFF_TYPE = "application/octet-stream"


@router.get("/hgt")
def hgt(path: str) -> Response:
    """The file at `path` as a GeoTIFF. See ``hgt_file.read_hgt``."""
    if not warmup.ready:
        raise HTTPException(status_code=503, detail="WARMING_UP")
    if not Path(path).is_file():
        raise HTTPException(status_code=404, detail=f"not a file: {path}")
    return Response(content=read_hgt(path), media_type=GEOTIFF_TYPE)
