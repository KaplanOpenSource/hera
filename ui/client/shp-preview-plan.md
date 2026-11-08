# Plan: shapefile preview with FlatGeobuf (#642)

## Today

- `.shp` documents get no preview. The map preview only takes `.geojson` / `.geo.json`.
- `GET /geojson` reads with geopandas, so it could read `.shp`, but GeoJSON is slow and big.
- Measured on `JERU-BLDG.shp` (156,679 buildings, really Tel Aviv area):
  - Read + convert to WGS84: 1.2 s.
  - GeoJSON: 5.4 s, 124 MB.
  - FlatGeobuf: 0.6 s, 59 MB; browser decode 0.24 s.

## Change

1. Server: `GET /flatgeobuf?path=...` reads the file with geopandas and returns FlatGeobuf bytes in WGS84.
2. Server: share the "read as WGS84, guess ITM from large numbers" step with `read_geojson`.
3. Client: add the `flatgeobuf` package; decode the bytes into GeoJSON features.
4. Client: a shapefile preview on the same Leaflet map as GeoJSON (tile server, tile button, fit bounds).
5. Client: draw with Leaflet's canvas renderer, so many shapes stay usable.
6. Client: add `.shp` to the preview file check and to `hasPreview`.

## Later, if big files are still slow

- Draw only what is on screen: write the FlatGeobuf with its spatial index and read by map box.
- Or send shapes only, without the table columns (35 MB instead of 59 MB).

## Tests

- Server: the route returns FlatGeobuf that reads back to the same row count and WGS84 bounds.
- Client: `.shp` gets a preview; decoded features reach the map.
- Then the full `TEST_UI.md` checklist.
