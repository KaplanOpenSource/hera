# Measurements Toolkits

Measurement toolkits work with raw observational data — geographic information, meteorological observations, and experimental recordings. Each toolkit provides specialized methods for loading, processing, and visualizing its domain.

All measurement toolkits are accessed via `toolkitHome.getToolkit()` and bound to a project:

```python
from hera import toolkitHome

toolkit = toolkitHome.getToolkit("GIS_Raster_Topography", projectName="MY_PROJECT")
```

---

## GIS Toolkits

GIS toolkits manage geographic data — elevation models, building footprints, land cover classification, and population data.

### Topography (Raster)

**Toolkit name:** `GIS_Raster_Topography`

Works with SRTM elevation data (HGT files) to provide terrain information for any location.

```python
topo = toolkitHome.getToolkit("GIS_Raster_Topography", projectName="MY_PROJECT")

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

```python
vtopo = toolkitHome.getToolkit("GIS_Vector_Topography", projectName="MY_PROJECT")

# Cut a region from a data source
region = vtopo.cutRegionFromSource(shape=my_polygon, datasource="contours")

# Generate STL from vector data
vtopo.regionToSTL(shape=my_polygon, dxdy=30, datasource="contours")
```

### Buildings

**Toolkit name:** `GIS_Buildings`

Manages building footprint data and generates 3D meshes for CFD simulations.

```python
buildings = toolkitHome.getToolkit("GIS_Buildings", projectName="MY_PROJECT")

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

```python
demo = toolkitHome.getToolkit("GIS_Demography", projectName="MY_PROJECT")

# Calculate population within a polygon
pop = demo.calculatePopulationInPolygon(polygon=my_area, datasource="census_2020")

# Create a new area GeoDataFrame
area_gdf = demo.createNewArea(polygon=my_area, datasource="census_2020")
```

### Land Cover

**Toolkit name:** `GIS_LandCover`

Land cover classification and surface roughness estimation for atmospheric modeling.

```python
lc = toolkitHome.getToolkit("GIS_LandCover", projectName="MY_PROJECT")

# Get land cover class at a point
cover = lc.getLandCoverAtPoint(lat=32.5, lon=35.2)

# Get gridded land cover for a region
ds = lc.getLandCover(minx=35.0, miny=32.0, maxx=35.1, maxy=32.1, dxdy=30)

# Get aerodynamic roughness length (z0) for CFD boundary conditions
z0 = lc.getRoughness(minx=35.0, miny=32.0, maxx=35.1, maxy=32.1, dxdy=30)
```

### Tiles

**Toolkit name:** `GIS_Tiles`

Tile-based raster data management for large datasets. Handles tiled storage and retrieval of raster data that is too large to fit in a single file.

```python
tiles = toolkitHome.getToolkit("GIS_Tiles", projectName="MY_PROJECT")

# List available tile data sources
tiles.getDataSourceList()
```

---

## Meteorology Toolkits

Meteorology toolkits handle weather station data at different temporal resolutions.

### MeteoLowFreq

**Toolkit name:** `MeteoLowFreq`

Hourly and daily meteorological station data with analysis and visualization.

```python
lf = toolkitHome.getToolkit("MeteoLowFreq", projectName="MY_PROJECT")

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
hf = toolkitHome.getToolkit("MeteoHighFreq", projectName="MY_PROJECT")

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
exp = toolkitHome.getToolkit("experiment", projectName="MY_PROJECT")

# List available experiments
exp.getExperimentsMap()

# Get a specific experiment
my_exp = exp.getExperiment(experimentName="wind_tunnel_march_2024")
```

For the full API details of each toolkit, see the [Toolkit Catalog](overview.md) and the [API Reference](../developer_guide/api/measurements.md).
