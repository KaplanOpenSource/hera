"""Reading an SRTM .hgt elevation file off disk as a GeoTIFF, for the TIFF preview."""
from __future__ import annotations

# SRTM marks a missing height with this value.
HGT_NO_DATA = -32768


def read_hgt(path: str) -> bytes:
    """The .hgt file as GeoTIFF bytes: full detail, compressed, missing heights as no-data."""
    import rasterio
    from rasterio.io import MemoryFile

    with rasterio.open(path) as source:
        heights = source.read(1)
        profile = source.profile
    profile.update(driver="GTiff", compress="deflate", predictor=2, nodata=HGT_NO_DATA)
    with MemoryFile() as memory:
        with memory.open(**profile) as target:
            target.write(heights, 1)
        return memory.read()
