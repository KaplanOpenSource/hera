# Measurements Toolkits

Measurement toolkits work with raw observational data — geographic information, meteorological observations, and experimental recordings. Each toolkit provides specialized methods for loading, processing, and visualizing its domain.

All measurement toolkits are accessed via `toolkitHome.getToolkit()` and bound to a project:

```python
from hera import toolkitHome

toolkit = toolkitHome.getToolkit(toolkitHome.GIS_RASTER_TOPOGRAPHY, projectName="MY_PROJECT")
```

---

## GIS Toolkits

GIS toolkits manage geographic data — elevation models, building footprints, land cover classification, and population data.

### Topography (Raster)

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

### Topography (Vector)

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

### Buildings

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

### Demography

**Toolkit name:** `GIS_Demography`

Population data analysis from census shapefiles.

**Data source format:**

| Property | Value |
|----------|-------|
| File format | Shapefile, GeoJSON, or GeoPackage (`dataFormat: "geopandas"`) |
| Geometry type | `Polygon` (census areas / statistical zones) |
| Required columns | `geometry`, `total_pop` (total population) |
| Optional columns | `age_0_14` (children), `age_15_19` (youth), `age_20_29` (young adults), `age_30_64` (adults), `age_65_up` (elderly) |
| CRS | Any (toolkit handles transformation; internally uses ITM EPSG:2039) |

The toolkit maps display names to column names: `"All"` → `total_pop`, `"Children"` → `age_0_14`, etc.

```python
demo = toolkitHome.getToolkit(toolkitHome.GIS_DEMOGRAPHY, projectName="MY_PROJECT")

# Calculate population within a polygon
pop = demo.calculatePopulationInPolygon(polygon=my_area, datasource="census_2020")

# Create a new area GeoDataFrame
area_gdf = demo.createNewArea(polygon=my_area, datasource="census_2020")
```

### Land Cover

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

### Tiles

**Toolkit name:** `GIS_Tiles`

The Tiles toolkit plots raster map images from a tile server (Google Maps, OpenStreetMap, or a custom server). This is useful for creating satellite or street-map backgrounds for your GIS visualizations.

**Data source format:**

| Property | Value |
|----------|-------|
| Data type | URL template string (`dataFormat: "string"`) |
| Format | XYZ tile URL with `{z}`, `{x}`, `{y}` placeholders |
| Example | `http://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}` |
| Tile format | Standard Web Mercator (XYZ) 256×256 PNG tiles |
| CRS | WGS84 (EPSG:4326) or ITM (EPSG:2039) for input coordinates |
| Config key | `defaultTileServer` |

```python
from hera import toolkitHome
import matplotlib.pyplot as plt

tiles = toolkitHome.getToolkit(toolkitHome.GIS_TILES, projectName="MY_PROJECT")
```

#### Getting a map image (WGS84 coordinates)

Specify a bounding box using WGS84 (lat/lon) coordinates:

```python
region = dict(
    minx=34.775, maxx=34.8,
    miny=32.05,  maxy=32.1,
    zoomlevel=17,
    tileServer=None  # use the default server
)
img = tiles.getImageFromCorners(**region)

# Plot the image
fig, ax = plt.subplots(1, 1, figsize=(12, 12))
tiles.presentation.plot(img, ax=ax)
```

#### Getting a map image (ITM coordinates)

You can also use Israeli Transverse Mercator (ITM) coordinates by specifying `inputCRS`:

```python
from hera.measurements.GIS.utils import ITM, WSG84

region = dict(
    minx=178898.481, maxx=181280.365,
    miny=661943.482, maxy=667478.9,
    zoomlevel=17,
    tileServer=None,
    inputCRS=ITM  # input coordinates are in ITM
)
img = tiles.getImageFromCorners(**region)

# Plot with ITM axes (default)
fig, ax = plt.subplots(1, 1, figsize=(12, 12))
tiles.presentation.plot(img, ax=ax)

# Or plot with WGS84 axes
fig, ax = plt.subplots(1, 1, figsize=(12, 12))
tiles.presentation.plot(img, outputCRS=WSG84, ax=ax)
```

#### Tile servers

The default server is Google Maps satellite imagery. You can check or change it:

```python
# Check the default
tiles.getConfig()['defaultTileServer']
# 'http://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}'

# Set a new default (e.g., OpenStreetMap)
tiles.setDefaultTileServer("https://tile.openstreetmap.org/{z}/{x}/{y}.png")

# Or pass a server for a single call without changing the default
region = dict(
    minx=34.775, maxx=34.8,
    miny=32.05, maxy=32.1,
    zoomlevel=17,
    tileServer="https://tile.openstreetmap.org/{z}/{x}/{y}.png"
)
img = tiles.getImageFromCorners(**region)
```

If you have a tile server registered as a data source in your repository, you can refer to it by name:

```python
tiles.getDataSourceList()
# ['localTileServer']

region = dict(
    minx=178898.481, maxx=181280.365,
    miny=661943.482, maxy=667478.9,
    zoomlevel=17,
    tileServer="localTileServer",  # use the registered data source
    inputCRS=ITM
)
img = tiles.getImageFromCorners(**region)
```

> **Jupyter tutorial:** See the [Tile Toolkit notebook](../tutorials/TileToolkitDocumentation.ipynb) for a full interactive walkthrough with map output.

---

## Meteorology Toolkits

Meteorology toolkits handle weather station data at different temporal resolutions.

### MeteoLowFreq

**Toolkit name:** `MeteoLowFreq`

Hourly and daily meteorological station data with analysis and visualization.

```python
lf = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_LOWFREQ, projectName="MY_PROJECT")

# Load station data
df = lf.getDataSourceData("YAVNEEL")

# Analysis: enrich with date columns (year, month, season, etc.)
enriched = lf.analysis.addDatesColumns(df, datecolumn="datetime")

# Analysis: hourly distribution of a variable
hourly = lf.analysis.calcHourlyDist(enriched, field="wind_speed", density=True)

# Presentation: daily scatter plot
lf.presentation.dailyPlots.plotScatter(enriched, plotField="temperature")

# Presentation: seasonal hourly distribution
lf.presentation.seasonalPlots.plotSeasonalHourly(enriched, field="wind_speed")

# Presentation: probability contour by season
lf.presentation.seasonalPlots.plotProbContourf_bySeason(enriched, fields=["T", "RH"])
```

### MeteoHighFreq

**Toolkit name:** `MeteoHighFreq`

High-frequency (10-20 Hz) sonic anemometer and TRH sensor data for turbulence analysis.

```python
hf = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_HIGHFREQ, projectName="MY_PROJECT")

# Load sonic anemometer data
sonic = hf.getDataSourceData("sonic_station_A")

# Analysis: calculate mean statistics
stats = hf.analysis.calculateMeanData(sonic)
```

---

## Experiment Toolkit

**Toolkit name:** `experiment`

Manages experimental workflows — organizing raw data files into structured experiments with metadata tracking.

```python
exp = toolkitHome.getToolkit(toolkitHome.EXPERIMENT, projectName="MY_PROJECT")

# List available experiments
exp.getExperimentsMap()

# Get a specific experiment
my_exp = exp.getExperiment(experimentName="wind_tunnel_march_2024")
```

For the full API details of each toolkit, see the [Toolkit Catalog](overview.md) and the [API Reference](../developer_guide/api/measurements.md).
