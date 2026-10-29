# Plan: show a GeoTIFF in the browser

## Idea

The browser reads the TIFF itself, a piece at a time, with `geotiff.js`.
No PNG, no server change: `/file/<path>` already answers byte ranges (206).

## First file

`hera/doc/jupyter/User/data/landcover/prefixaac.tif` (projects `a111112`, `a2222`):
21600 x 8962, one 8-bit band, deflate, 1024 x 1024 tiles, no overviews, no geo tags.
So step one shows it as a plain image, with no base map.

## Client

1. Add `geotiff` (npm).
2. New `maps/GeoTiffPreview.tsx`:
   - `isGeoTiffResource`: path ends in `.tif` or `.tiff`.
   - `fromUrl(BASEURL + '/file' + path)`: reads only the header first.
   - Leaflet map with `CRS.Simple` (plain pixels), bounds = image width and height.
   - New `maps/GeoTiffGridLayer.ts`: a Leaflet `GridLayer` whose `createTile` calls
     `image.readRasters({ window, width: 256, height: 256 })` for that tile's area.
     geotiff.js fetches only the blocks that area needs.
   - Color: grey from the band's min to max for now.
3. `PreviewPanel`: show `GeoTiffPreview` for a TIFF resource.
4. Tests: `geoTiffPreview.test.tsx` and `geoTiffGridLayer.test.ts` (window math, color), with `geotiff` mocked.

## Known limit

With no overviews, a zoomed-out view reads every block (about 11 MB here). Fine for this file.

## Later

- Files with geo tags: place them on the base map, with the tile server button.
- Overviews / cloud optimized GeoTIFF for zoomed-out views.
- Land cover class colors.
