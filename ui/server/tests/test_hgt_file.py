"""Reading SRTM .hgt files as GeoTIFF for the preview."""
from pathlib import Path

import numpy as np
import pytest
from rasterio.io import MemoryFile

import dependencies
import hgt_file
from fastapi import HTTPException
from routers import hgt

# The smaller SRTM size (3 arc-seconds); the real files are 3601.
SIZE = 1201
REPO_HGT = Path(__file__).resolve().parents[3] / "hera" / "measurements" / "GIS" / "raster" / "test_data" / "N33E035.hgt"


@pytest.fixture
def warmed(monkeypatch):
    monkeypatch.setattr(dependencies.warmup, "_ready", True)


def _heights():
    heights = (np.arange(SIZE * SIZE, dtype=np.int32) % 3000).astype(np.int16).reshape(SIZE, SIZE)
    heights[10, 20] = hgt_file.HGT_NO_DATA
    return heights


def _write(tmp_path, heights, name="N32E035.hgt"):
    path = tmp_path / name
    # .hgt is raw big-endian 16-bit heights, north row first, no header.
    path.write_bytes(heights.astype(">i2").tobytes())
    return str(path)


def _open(data):
    memory = MemoryFile(data)
    return memory.open()


def test_heights_are_exact_at_full_detail(tmp_path):
    heights = _heights()

    with _open(hgt_file.read_hgt(_write(tmp_path, heights))) as tiff:
        assert tiff.driver == "GTiff"
        assert (tiff.width, tiff.height) == (SIZE, SIZE)
        assert np.array_equal(tiff.read(1), heights)


def test_position_comes_from_the_name(tmp_path):
    with _open(hgt_file.read_hgt(_write(tmp_path, _heights()))) as tiff:
        half = 0.5 / (SIZE - 1)
        assert tiff.crs.to_epsg() == 4326
        assert tiff.bounds.left == pytest.approx(35 - half)
        assert tiff.bounds.bottom == pytest.approx(32 - half)
        assert tiff.bounds.right == pytest.approx(36 + half)
        assert tiff.bounds.top == pytest.approx(33 + half)


def test_missing_heights_stay_no_data(tmp_path):
    with _open(hgt_file.read_hgt(_write(tmp_path, _heights()))) as tiff:
        assert tiff.nodata == hgt_file.HGT_NO_DATA
        assert tiff.read(1, masked=True).mask[10, 20]


def test_it_is_compressed(tmp_path):
    data = hgt_file.read_hgt(_write(tmp_path, _heights()))

    with _open(data) as tiff:
        assert tiff.compression.value == "DEFLATE"
    assert len(data) < SIZE * SIZE * 2


def test_the_repo_file_matches_its_raw_heights():
    raw = np.fromfile(REPO_HGT, dtype=">i2")
    side = int(np.sqrt(raw.size))
    raw = raw.reshape(side, side)

    with _open(hgt_file.read_hgt(str(REPO_HGT))) as tiff:
        heights = tiff.read(1)
    for row, col in [(0, 0), (side // 2, side // 3), (side - 1, side - 1), (123, 3000)]:
        assert heights[row, col] == raw[row, col]


def test_the_route_sends_the_geotiff(tmp_path, warmed):
    heights = _heights()

    response = hgt.hgt(_write(tmp_path, heights))

    assert response.media_type == hgt.GEOTIFF_TYPE
    with _open(response.body) as tiff:
        assert np.array_equal(tiff.read(1), heights)


def test_the_route_reports_a_missing_file(tmp_path, warmed):
    with pytest.raises(HTTPException) as raised:
        hgt.hgt(str(tmp_path / "N00E000.hgt"))

    assert raised.value.status_code == 404


def test_the_route_waits_for_warmup(tmp_path, monkeypatch):
    monkeypatch.setattr(dependencies.warmup, "_ready", False)

    with pytest.raises(HTTPException) as raised:
        hgt.hgt(_write(tmp_path, _heights()))

    assert raised.value.status_code == 503
