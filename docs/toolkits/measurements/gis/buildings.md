# Buildings

**Toolkit name:** `GIS_Buildings`

Manages building footprint data — query by bounding box, extract heights, compute urban morphology parameters (plan area fraction, frontal area density, average building height), generate 3D STL meshes for CFD simulations, and group buildings into convex clusters.

```python
from hera import toolkitHome

# Tip: if you created the project with `hera-project project create`, you can omit projectName
buildings = toolkitHome.getToolkit(toolkitHome.GIS_BUILDINGS, projectName="MY_PROJECT")

# Get building footprints in a bounding box
gdf = buildings.getBuildingsFromRectangle(
    minx=35.0, miny=32.0, maxx=35.1, maxy=32.1
)

# Compute urban morphology (lambda_P, lambda_F, hc) per block
morphology = buildings.analysis.LambdaFromBuildingData(
    windMeteorologicalDirection=270,
    resolution=250,
    buildingsData=gdf,
)
```

For the full API, see the [API Reference](../../developer_guide/api/measurements.md). For implementation details (height resolution logic, Blocks class, cache strategy), see the [Developer Guide > GIS > BuildingsToolkit](../../developer_guide/measurements/gis.md#buildingstoolkit).

---

## Data source format

| Property | Value |
|----------|-------|
| File format | Shapefile, GeoJSON, or GeoPackage (`dataFormat: "geopandas"`) |
| Geometry type | `Polygon` (building footprints) |
| Required columns | `geometry` |
| Height column | `BLDG_HT` (or as specified in `desc.BuildingHeightColumn`) — building height in meters |
| Ground height | `HT_LAND` (or as specified in `desc.LandHeightColumns`) — ground elevation |
| Fallback | If no height column, uses `building:levels * 3` meters |
| CRS | Any (must be defined; methods handle transformation) |
| Config key | `defaultBuildingDataSource` |

---

## Querying buildings

### By bounding box

```python
# WGS84 coordinates (default inputCRS)
gdf = buildings.getBuildingsFromRectangle(
    minx=35.0, miny=32.0, maxx=35.1, maxy=32.1,
)

# ITM coordinates
from hera.measurements.GIS import ITM
gdf = buildings.getBuildingsFromRectangle(
    minx=175000, miny=658000, maxx=185000, maxy=668000,
    inputCRS=ITM,
)

# With elevation from raster topography toolkit
gdf = buildings.getBuildingsFromRectangle(
    minx=35.0, miny=32.0, maxx=35.1, maxy=32.1,
    withElevation=True,   # adds 'elevation' column from SRTM data
)

# From a specific data source
gdf = buildings.getBuildingsFromRectangle(
    minx=35.0, miny=32.0, maxx=35.1, maxy=32.1,
    dataSourceName="BNTL",
)
```

### From GeoJSON

```python
# Filter buildings from a GeoJSON-like dict
gdf = BuildingsToolkit.filter_buildings_in_area(
    buildings_data=geojson_dict,
    min_longitude=35.0, min_latitude=32.0,
    max_longitude=35.1, max_latitude=32.1,
)
```

---

## Extracting building heights

```python
# Extract heights from a GeoDataFrame
# Uses BLDG_HT column if available; falls back to building:levels * 3
height_gdf = BuildingsToolkit.get_buildings_height(gdf)
# Returns GeoDataFrame with columns: name, geometry, height
```

### Height with raster topography

```python
# Add ground elevation to each building (at centroid)
gdf_with_elev = buildings.getBuildingHeightFromRasterTopographyToolkit(
    buildingData=gdf,
    topographyDataSource="SRTM_30m",   # optional, uses default if None
)
# Adds 'elevation' column from the TopographyToolkit
```

---

## Analysis: urban morphology

The analysis layer computes block-averaged urban morphology parameters used in boundary-layer meteorology and CFD modeling.

### Lambda parameters

Compute plan area fraction (**lambda_P**), frontal area density (**lambda_F**), and average building height (**hc**) for each block in a domain:

```python
morphology = buildings.analysis.LambdaFromBuildingData(
    windMeteorologicalDirection=270,   # wind direction in degrees
    resolution=250,                     # block size in meters
    buildingsData=gdf,                  # GeoDataFrame with CRS set
    overwrite=False,                    # recompute even if cached
    saveCache=True,                     # cache results in DB
)
```

**Returns** a GeoDataFrame with one row per block:

| Column | Description |
|--------|-------------|
| `geometry` | Block polygon |
| `lambdaP` | Plan area fraction: total building footprint area / block area |
| `lambdaF` | Frontal area density: total building frontal area / block area (wind-direction dependent) |
| `hc` | Area-weighted average building height within the block |
| `i0`, `j0` | Block grid indices |

Results are cached in the project's cache collection (type `Lambda_Buildings`). Subsequent calls with the same bounds, wind direction, resolution, and CRS return the cached result. Use `overwrite=True` to force recomputation.

### How it works

1. The domain (from `buildingsData.total_bounds`) is divided into a grid of square blocks at the given resolution
2. For each block:
   - **lambda_P** = sum of building plan areas / block area (excludes buildings with `FTYPE` 14 or 16, and buildings with zero height)
   - **lambda_F** = sum of building frontal areas / block area (frontal area is computed by rotating each building footprint to the wind direction and measuring the projected width * height)
   - **hc** = sum(area * height) / sum(area) — area-weighted mean building height
3. Buildings with no land height (`HT_LAND=0`) use the nearest building's land height as a fallback

### Convex building clusters

Group nearby buildings into convex hulls (useful for identifying urban blocks):

```python
clusters = buildings.analysis.ConvexPolygons(
    regionNameOrData=gdf,   # GeoDataFrame or region name
    buffer=100,              # buffer distance for grouping (meters)
)
# Returns GeoDataFrame of convex hull polygons, sorted by area
```

The algorithm:
1. Buffers each building footprint by `buffer` meters
2. Groups intersecting buffered footprints
3. Computes the convex hull of each group
4. Recursively merges overlapping hulls

---

## 3D STL generation

Generate 3D building meshes for CFD simulations using the FreeCAD module.

### From a bounding box (recommended)

```python
buildings.regionToSTL(
    minx=35.0, miny=32.0, maxx=35.1, maxy=32.1,
    dxdy=30,                  # grid resolution
    inputCRS=4326,            # input coordinate system
    outputCRS=2039,           # output coordinate system (ITM)
    solidName="Buildings",
    fileName="buildings.stl",
)
```

### From a GeoDataFrame (low-level)

```python
buildings.buildingsGeopandasToSTLRasterTopography(
    buildingData=gdf,
    buildingHeightColumn="BLDG_HT",
    buildingElevationColumn="elevation",
    outputFileName="/output/buildings.stl",
    flatTerrain=False,             # True: use referenceTopography for all
    referenceTopography=0,         # base height if flatTerrain=True
    nonFlatTopographyShift=10,     # shift buildings above topography
)
```

!!! note "FreeCAD required"
    STL generation requires the FreeCAD Python module. Install it and add to `PYTHONPATH` before using these methods.

---

## Presentation layer

The buildings toolkit has a presentation layer for visualising footprints, heights, and morphology parameters.

**Available methods:**

| Method | What it plots | Key parameters |
|--------|-------------|----------------|
| `plotBuildings(data)` | Building footprints (uniform color) | `facecolor`, `alpha` |
| `plotBuildingHeights(data)` | Footprints colored by height (choropleth) | `heightColumn`, `cmap`, `vmin`/`vmax` |
| `plotMorphology(morphData, column)` | Block grid colored by λp, λf, or hc | `column` (`"lambdaP"`, `"lambdaF"`, `"hc"`) |
| `plotBuildingsWithMorphology(buildings, morph)` | Buildings overlaid on morphology grid | `morphColumn`, building + morph styling |
| `plotHeightHistogram(data)` | Histogram of building heights | `bins`, `heightColumn` |
| `plotBuildingsOnMap(data, tilesToolkit)` | Buildings overlaid on tile server map | `heightColumn`, `zoomlevel` |

All map methods support: `ax`, `figsize`, `cmap`, `alpha`, `edgecolor`, `linewidth`, `colorbar`, `title`, `xlim`/`ylim`, `inputCRS`, `outputCRS`.

### Building footprints

```python
# Simple footprint map
ax = buildings.presentation.plotBuildings(gdf, facecolor="steelblue", alpha=0.8)

# Colored by height
ax = buildings.presentation.plotBuildingHeights(
    gdf, heightColumn="BLDG_HT",
    cmap="YlOrBr", vmin=0, vmax=50,
    title="Building Heights",
)
```

### Morphology choropleth

Visualise the output of `LambdaFromBuildingData` as a colored block grid:

```python
morphology = buildings.analysis.LambdaFromBuildingData(270, 250, gdf)

# Plan area fraction
ax = buildings.presentation.plotMorphology(
    morphology, column="lambdaP",
    cmap="RdYlGn_r", title="Plan Area Fraction (λp)",
)

# Frontal area density
ax = buildings.presentation.plotMorphology(
    morphology, column="lambdaF",
    title="Frontal Area Density (λf)",
)

# Average building height
ax = buildings.presentation.plotMorphology(
    morphology, column="hc",
    cmap="YlOrBr", colorbar_label="Average height (m)",
)
```

### Buildings overlaid on morphology

```python
ax = buildings.presentation.plotBuildingsWithMorphology(
    gdf, morphology,
    morphColumn="lambdaP",
    buildingColor="0.2", buildingAlpha=0.6,
    morphAlpha=0.5,
    title="Buildings with λp blocks",
)
```

### Height histogram

```python
ax = buildings.presentation.plotHeightHistogram(
    gdf, heightColumn="BLDG_HT",
    bins=40, title="Building Height Distribution",
)
```

### Buildings on map

```python
from hera import toolkitHome

tiles = toolkitHome.getToolkit(toolkitHome.GIS_TILES, projectName="MY_PROJECT")

# Uniform color on map
ax = buildings.presentation.plotBuildingsOnMap(gdf, tiles, zoomlevel=15)

# Colored by height on map
ax = buildings.presentation.plotBuildingsOnMap(
    gdf, tiles,
    heightColumn="BLDG_HT",
    cmap="YlOrBr", alpha=0.7, zoomlevel=15,
    title="Building Heights on Map",
)
```

---

## Complete example

```python
from hera import toolkitHome
from hera.measurements.GIS import ITM

# Tip: if you created the project with `hera-project project create`, you can omit projectName
buildings = toolkitHome.getToolkit(toolkitHome.GIS_BUILDINGS, projectName="TelAvivStudy")

# 1. Query buildings in Tel Aviv (ITM coordinates)
gdf = buildings.getBuildingsFromRectangle(
    minx=175000, miny=658000, maxx=185000, maxy=668000,
    inputCRS=ITM,
    withElevation=True,
)
print(f"Found {len(gdf)} buildings")

# 2. Visualise footprints colored by height
buildings.presentation.plotBuildingHeights(gdf, title="Tel Aviv Buildings")

# 3. Height distribution
buildings.presentation.plotHeightHistogram(gdf, bins=40)

# 4. Compute urban morphology for westerly wind
morphology = buildings.analysis.LambdaFromBuildingData(
    windMeteorologicalDirection=270,
    resolution=250,
    buildingsData=gdf,
)

# 5. Visualise morphology
buildings.presentation.plotBuildingsWithMorphology(
    gdf, morphology, morphColumn="lambdaP",
    title="Tel Aviv — Plan Area Fraction",
)

# 6. Buildings on satellite map
tiles = toolkitHome.getToolkit(toolkitHome.GIS_TILES, projectName="TelAvivStudy")
buildings.presentation.plotBuildingsOnMap(
    gdf, tiles, heightColumn="BLDG_HT", zoomlevel=15,
)

# 7. Generate 3D STL for CFD
buildings.buildingsGeopandasToSTLRasterTopography(
    buildingData=gdf,
    buildingHeightColumn="BLDG_HT",
    buildingElevationColumn="elevation",
    outputFileName="tel_aviv_buildings.stl",
    flatTerrain=False,
)
```
