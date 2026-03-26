# Topography (Raster)

**Toolkit name:** `GIS_Raster_Topography`

Works with SRTM elevation data (HGT files) to provide terrain information for any location.

**Data source format:**

| Property | Value |
|----------|-------|
| File format | SRTM `.hgt` binary files (30m resolution) |
| File naming | `N{lat}E{lon}.hgt` (e.g., `N33E035.hgt`) |
| Structure | 1201x1201 pixels, 2 bytes per sample, big-endian signed short |
| CRS | WGS84 (EPSG:4326) |
| Data source | Folder path containing `.hgt` files |
| Config key | `defaultSRTM` |

```python
topo = toolkitHome.getToolkit(toolkitHome.GIS_RASTER_TOPOGRAPHY, projectName="MY_PROJECT")

# Get elevation at a single point
elev = topo.getPointElevation(lat=32.5, long=35.2)

# Get a gridded elevation dataset for a bounding box
ds = topo.getElevation(minx=35.0, miny=32.0, maxx=35.1, maxy=32.1, dxdy=30)

# Generate an STL mesh for CFD simulations
topo.createElevationSTL(ds, solidName="Terrain", fileName="terrain.stl")

# Coordinate transformations
points = topo.convertPointsCRS(points, inputCRS=4326, outputCRS=2039)
```

For the full API, see the [API Reference](../../developer_guide/api/measurements.md).
