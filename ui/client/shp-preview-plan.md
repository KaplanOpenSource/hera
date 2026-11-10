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

## Phase 1.25: dots above a cap, no clusters

To see how plain dots look before adding Supercluster.

1. Server: read the shape count first (`pyogrio.read_info`, about 2 ms).
2. Server: over a cap (`MAX_SHAPES`, a constant, start at 50,000), send one dot per shape at its center, no columns.
3. Server: find the centers in the file's own CRS (ITM), then convert to WGS84.
4. Server: at or under the cap, send the full shapes, as in phase 1.
5. Client: draw points as small canvas circles (`circleMarker`), not Leaflet's default pin icons.
6. Measured on Jerusalem: 0.77 s and 10 MB as dots, against 2.6 s and 59 MB as shapes.
7. Tests: over the cap gives one dot per shape at the exact center, no columns; at the cap gives full shapes; ITM with no CRS still works.

## Phase 1.4: dots and clusters for big files (Supercluster, client only)

- Big files are slow as full shapes: 2.6 s and 59 MB for the Jerusalem buildings.
- As dots (one per shape, no table columns) they take 0.54 s and 10 MB.
- Leaflet alone may be slow with 157k dots, so merge near dots into clusters with a count.
- Seen: the browser hangs when switching from Jerusalem (157k dots as Leaflet layers) to a lambda geojson. Clusters should fix it; not confirmed.

1. Server: no change. It already sends dots (phase 1.25).
2. Client: add the `supercluster` package (https://github.com/mapbox/supercluster). Checked: no network code, one small dependency (`kdbush`), ISC license.
3. Client: when the answer is many points, cluster them with Supercluster; redraw the clusters on zoom / pan.
4. Client: point files with many points get clusters too, the same way.
5. Zoomed in, clusters split into single dots. Full shapes come in phase 1.7.

## Phase 1.7: full shapes when zoomed in

- A big file shows dots (phase 1.25). Zoomed in far enough, show full shapes instead.
- One limit for both: `MAX_SHAPES` = 10,000 (was 50,000), on the server only.
  - First open: at or under it, full shapes; over it, dots.
  - Zoomed in: the view holds at most that many dots, so ask for full shapes.
- A flag says the answer is dots, in the FlatGeobuf title: `dots max=10000`. Full shapes have no title.
  - Checked: GDAL writes it (`TITLE=...`); the `flatgeobuf` package reads it (header `title`).
  - The client gets the limit from the flag, so it keeps no copy of the number.

1. Server: set `MAX_SHAPES` to 10,000; write the dots flag in the title when sending dots.
2. Client: read the title while decoding; dots only when the flag is there.
3. Client: on zoom / pan, count the dots inside the map view. No server call: the dots are already loaded.
4. Client: at or under the limit, ask for the full shapes in the view's bbox.
5. Server: `GET /geoshapes?path=...&bbox=minx,miny,maxx,maxy`, bbox in WGS84; full shapes in the bbox only, no flag.
6. Server: no bbox works as before (whole file, or dots over the limit).
7. Server: convert the bbox to the file's own CRS before reading (ITM guess as in phase 1).
8. Server: handle a bbox with no shapes. Writing an empty FlatGeobuf to memory fails today.
9. Client: draw the shapes and hide the clusters; over the limit again, back to clusters.
10. Client: drop a reply for an old view if the map moved since.
11. Measured on Jerusalem: 0.05 s for a tiny bbox, 0.23 s for a 5 km bbox (12,621 shapes). No index needed.
12. Tests: the flag is set for dots and not for shapes; a bbox read gives exactly the shapes in the bbox; an empty bbox gives none; a stale reply is dropped.

## Phase 2: to plan after phase 1.7

Ideas only. Plan the details once phase 1.7 is done.

- No cache needed for `.shp`: bbox reads are fast (phase 1.7).
- `.geojson` has no per-shape boxes, so a bbox read parses the whole file. Not timed yet.

## Future: `.gpkg` and similar files

- `.gpkg` (GeoPackage) and similar formats can be many GB, so never read or convert them whole.
- Most have their own spatial index, so a map bbox read (phase 1.7) is fast.
- Needed: add the extension to the preview file check and `hasPreview`, plus a test.
- A `.gpkg` can hold several layers. Pick one (first by default) or let the user choose.

## Tests (phase 1)

- Server: `.shp` / `.geojson` / `.fgb` return FlatGeobuf with the right row count and WGS84 bounds.
- Server: an ITM file with no CRS still gives WGS84 shapes.
- Server: exact shapes. Read the reply back with geopandas and compare to the source in WGS84: same rows in the same order, same geometry (within a tiny tolerance), same column values.
- Server: exact shapes for each type: points, lines, polygons, multi-polygons, and polygons with holes.
- Client: `.shp` / `.fgb` get a preview; decoded features reach the map.
- Then the full `TEST_UI.md` checklist.
