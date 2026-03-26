# Topography (Vector)

**Toolkit name:** `GIS_Vector_Topography`

Works with vector-based topography — contour lines and survey points.

**Data source format:**

| Property | Value |
|----------|-------|
| File format | Shapefile, GeoJSON, or GeoPackage |
| Geometry type | `LineString` or `MultiLineString` (contour lines) |
| Required columns | `geometry`, height column (default: `HEIGHT`) |
| Height column | Elevation value in meters for each contour line |
| CRS | Any valid CRS (toolkit handles conversion) |

```python
vtopo = toolkitHome.getToolkit(toolkitHome.GIS_VECTOR_TOPOGRAPHY, projectName="MY_PROJECT")

# Cut a region from a data source
region = vtopo.cutRegionFromSource(shape=my_polygon, datasource="contours")

# Generate STL from vector data
vtopo.regionToSTL(shape=my_polygon, dxdy=30, datasource="contours")
```

For the full API, see the [API Reference](../../developer_guide/api/measurements.md).
