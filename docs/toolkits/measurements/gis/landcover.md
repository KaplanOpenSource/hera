# Land Cover

**Toolkit name:** `GIS_LandCover`

Land cover classification and surface roughness estimation for atmospheric modeling.

**Data source format:**

| Property | Value |
|----------|-------|
| File format | GeoTIFF (`dataFormat: "geotiff"`) |
| Structure | Single band (band 1), UINT8 values 0–16 (IGBP classification) |
| Classes | 0: Water, 1–5: Forests, 6–9: Shrublands/Savannas, 10–12: Grasslands/Crops, 13: Urban, 14: Mosaic, 15: Snow, 16: Barren |
| CRS | WGS84 (EPSG:4326) |
| Config key | `defaultLandCover` |

The toolkit maps IGBP land cover classes to aerodynamic roughness lengths (z0) for atmospheric simulations.

```python
lc = toolkitHome.getToolkit(toolkitHome.GIS_LANDCOVER, projectName="MY_PROJECT")

# Get land cover class at a point
cover = lc.getLandCoverAtPoint(lat=32.5, lon=35.2)

# Get gridded land cover for a region
ds = lc.getLandCover(minx=35.0, miny=32.0, maxx=35.1, maxy=32.1, dxdy=30)

# Get aerodynamic roughness length (z0) for CFD boundary conditions
z0 = lc.getRoughness(minx=35.0, miny=32.0, maxx=35.1, maxy=32.1, dxdy=30)
```

For the full API, see the [API Reference](../../developer_guide/api/measurements.md).
