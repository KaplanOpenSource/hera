# Plan: shapefile preview with FlatGeobuf (#642)

## Today

- `.shp` documents get no preview. The map preview only takes `.geojson` / `.geo.json`.
- `GET /geojson` reads with geopandas, so it could read `.shp`, but GeoJSON is slow and big.
- The client fetches the file once on open; pan / zoom are local in Leaflet.
- Measured on `JERU-BLDG.shp` (156,679 buildings, really Tel Aviv area):
  - Read + convert to WGS84: 1.2 s.
  - GeoJSON: 5.4 s, 124 MB.
  - FlatGeobuf: 0.6 s, 59 MB; browser decode 0.24 s.
- Goal is `.shp` and `.fgb`. `.gpkg` is for later; see "Future".

## One route: `GET /geoshapes`

- `GET /geoshapes?path=...`
- Takes `.shp`, `.geojson`, `.fgb`.
- Always returns FlatGeobuf bytes in WGS84.
- Replaces `GET /geojson`. The client has one way to fetch and decode.

## Phase 1: read and convert the whole file

1. Server: `/geoshapes` reads the whole file with geopandas, converts to WGS84, writes FlatGeobuf to memory.
2. Server: keep the "large numbers mean ITM" guess for files with no CRS.
3. Server: remove `/geojson` and move its tests to `/geoshapes`.
4. Client: add the `flatgeobuf` package; decode the bytes into GeoJSON features.
5. Client: fetch once on open, as today; pan / zoom stay local.
6. Client: one map preview for all these files (tile server, tile button, fit bounds).
7. Client: draw with Leaflet's canvas renderer, so many shapes stay usable.
8. Client: add `.shp` and `.fgb` to the preview file check and to `hasPreview`.

## Phase 2: to plan after phase 1

Ideas only. Plan the details once phase 1 is done.

- Each open reads and converts again (about 2 s on the big buildings file).
- A memory cache of the last few files, so reopen is fast. Must be thread safe.
- For huge files, read only the map bbox, with a cap on shape count and a "zoom in" message.

## Future: `.gpkg` and similar files

- `.gpkg` (GeoPackage) and similar formats can be many GB, so never read or convert them whole.
- Most have their own spatial index, so a map bbox read (phase 2) is fast.
- Needed: add the extension to the preview file check and `hasPreview`, plus a test.
- A `.gpkg` can hold several layers. Pick one (first by default) or let the user choose.

## Tests (phase 1)

- Server: `.shp` / `.geojson` / `.fgb` return FlatGeobuf with the right row count and WGS84 bounds.
- Server: an ITM file with no CRS still gives WGS84 shapes.
- Server: exact shapes. Read the reply back with geopandas and compare to the source in WGS84: same rows in the same order, same geometry (within a tiny tolerance), same column values.
- Server: exact shapes for each type: points, lines, polygons, multi-polygons, and polygons with holes.
- Client: `.shp` / `.fgb` get a preview; decoded features reach the map.
- Then the full `TEST_UI.md` checklist.
