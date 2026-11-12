# Plan: elevation (.hgt) preview

## Today

- `.hgt` files get no preview. The preview knows `.tif` / `.tiff` only.
- `.hgt` is SRTM elevation: one 1° x 1° tile, 16-bit heights, no header. Missing heights are -32768.
- Test file: `hera/measurements/GIS/raster/test_data/N33E035.hgt` (25 MB, 3601 x 3601).
- The TIFF preview reads files in the browser with `geotiff.js`, which cannot read `.hgt`.

## Phase 1

1. Server: `GET /hgt?path=...` opens the `.hgt` file with `rasterio` (already in the requirements).
2. Server: keep -32768 as no-data.
3. Server: send it as a GeoTIFF written in memory, full detail, compressed. Nothing is saved.
4. Client: `isHgtResource`: the resource ends in `.hgt`.
5. Client: the TIFF preview opens the `/hgt` bytes with `fromArrayBuffer`.
6. Client: add `.hgt` to `hasPreview` and `PreviewPanel`.
7. Measure on `N33E035.hgt`: server time, bytes, browser time.

## Later

- Shrink before sending, if full detail is too big: long side at most `MAX_PREVIEW_SIZE` (start at 2048), averaging heights.
- On map tiles in the right place, not a plain pixel grid.
- Colors for height, not grey.

## Tests

- Server: right size, position and heights (a few points against the raw file).
- Server: -32768 stays no-data.
- Server: a missing path gives a 404.
- Client: a `.hgt` file gets a preview; the bytes reach the TIFF preview.
- Then the full `TEST_UI.md` checklist.
