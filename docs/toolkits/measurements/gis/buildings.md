# Buildings

**Toolkit name:** `GIS_Buildings`

Manages building footprint data and generates 3D meshes for CFD simulations.

**Data source format:**

| Property | Value |
|----------|-------|
| File format | Shapefile, GeoJSON, or GeoPackage (`dataFormat: "geopandas"`) |
| Geometry type | `Polygon` (building footprints) |
| Required columns | `geometry` |
| Height column | `BLDG_HT` (or as specified in `desc.BuildingHeightColumn`) — building height in meters |
| Ground height | `HT_LAND` (or as specified in `desc.LandHeightColumns`) — ground elevation |
| Fallback | If no height column, uses `building:levels × 3` meters |
| CRS | Any (must be defined; methods handle transformation) |
| Config key | `defaultBuildingDataSource` |

```python
buildings = toolkitHome.getToolkit(toolkitHome.GIS_BUILDINGS, projectName="MY_PROJECT")

# Get building footprints in a bounding box
gdf = buildings.getBuildingsFromRectangle(minx=35.0, miny=32.0, maxx=35.1, maxy=32.1)

# Generate 3D STL using raster topography for ground elevation
buildings.regionToSTL(
    minx=35.0, miny=32.0, maxx=35.1, maxy=32.1,
    dxdy=30, inputCRS=4326, outputCRS=2039,
    solidName="Buildings", fileName="buildings.stl"
)
```

For the full API, see the [API Reference](../../developer_guide/api/measurements.md).
