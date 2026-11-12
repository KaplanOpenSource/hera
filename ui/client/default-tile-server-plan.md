# Plan: default tile server for the GeoJSON polygon map (#642)

## Today

- `GeoJsonPreview` draws polygons on a fixed OpenStreetMap tile URL (`TILE_URL` in `maps/mapDarkStyles.ts`).
- The project config can hold `defaultTileServer` (Cache doc of type `<projectName>__config__`, key `desc.defaultTileServer`).
- The client store already has that doc. No server change.

## Change

1. `GeoJsonPreview`: read `desc.defaultTileServer` from the project config doc in `useProjectStore`.
2. If it is set, use it as the tile URL. Else use `TILE_URL`.

## Tests

- GeoJSON preview test: config set uses that URL; no config uses OSM.
- Then the full `TEST_UI.md` checklist.
