# Demography

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
pop = demo.analysis.calculatePopulationInPolygon(
    shapelyPolygon=my_area,
    dataSourceOrData="census_2020"
)

# Create a new area GeoDataFrame
area_gdf = demo.analysis.createNewArea(
    shapeNameOrData=my_area,
    dataSourceOrData="census_2020"
)
```

## Presentation layer

The demography toolkit has a presentation layer for visualizing population data on maps.

```python
from hera.measurements.GIS import WSG84, ITM
from hera.utils.unitHandler import ureg
```

**Available methods:**

| Method | What it plots | Key parameters |
|--------|-------------|----------------|
| `plotPopulationDensity(data)` | Population per area unit (choropleth) | `density_units` for unit control |
| `plotPopulation(data)` | Absolute population counts (choropleth) | `populationType` for column selection |
| `plotPopulationByType(data)` | Grid of subplots for all age groups | `ncols` for grid layout |
| `plotPopulationInPolygon(result)` | Intersection result with query polygon outline | `queryPolygon`, `contextData` |
| `plotArea(area)` | Custom area with population annotation | `annotate`, `contextData` |
| `plotPopulationOnMap(data, tiles)` | Population overlaid on tile server map | `density`, `density_units`, `alpha`, `zoomlevel` |

All methods support: `ax`, `figsize`, `cmap`, `vmin`/`vmax`, `alpha`, `edgecolor`, `linewidth`, `colorbar`, `title`, `xlim`/`ylim`, `inputCRS`, `outputCRS`.

**Plot population density** (people per area unit):

```python
census = demo.getDataSourceData("census_2020")

# Default: people per km² in ITM coordinates
ax = demo.presentation.plotPopulationDensity(census)

# Full control — density in dunam, custom colormap and domain
ax = demo.presentation.plotPopulationDensity(
    census,
    populationType="total_pop",
    density_units=ureg.dunam,     # people per dunam (1000 m²)
    inputCRS=WSG84,               # data is in WGS84
    outputCRS=ITM,                # plot in ITM (default)
    cmap="YlOrRd",
    vmin=0, vmax=50,
    alpha=0.8,
    xlim=(170000, 190000),
    ylim=(660000, 670000),
    title="Population Density"
)
```

**Density unit options** (using pint):

| Unit | Code | Description |
|------|------|-------------|
| km² | `ureg.km**2` | People per square kilometer (default) |
| dunam | `ureg.dunam` | People per dunam (1000 m²) |
| hectare | `ureg.hectare` | People per hectare (10,000 m²) |
| m² | `ureg.m**2` | People per square meter |

The colorbar label is auto-generated from the unit (e.g., "Population density [people/dunam]").

**Plot absolute population counts**:

```python
ax = demo.presentation.plotPopulation(
    census,
    populationType="total_pop",
    cmap="Blues",
    outputCRS=ITM
)
```

**Plot all age groups** as a grid of subplots:

```python
fig = demo.presentation.plotPopulationByType(
    census,
    ncols=3,
    cmap="YlOrRd",
    alpha=0.8,
    inputCRS=WSG84,
    outputCRS=ITM
)
```

**Visualize a polygon intersection result:**

```python
# Calculate population in a polygon
result = demo.analysis.calculatePopulationInPolygon(
    shapelyPolygon=my_polygon, dataSourceOrData="census_2020"
)

# Plot with query polygon outline and surrounding census context
ax = demo.presentation.plotPopulationInPolygon(
    result,
    queryPolygon=my_polygon,       # dashed blue outline
    contextData=census,            # gray background
    alpha=0.8, outputCRS=ITM
)
```

**Visualize a custom area with population annotation:**

```python
area = demo.analysis.createNewArea(
    shapeNameOrData=my_polygon, dataSourceOrData="census_2020"
)
ax = demo.presentation.plotArea(
    area.getData(),
    contextData=census,
    annotate=True,                 # shows population number on the area
    annotate_fontsize=14,
    outputCRS=ITM
)
```

**Overlay population on a tile server map:**

```python
tiles = toolkitHome.getToolkit(toolkitHome.GIS_TILES, projectName="MY_PROJECT")

# Density overlay on satellite map
ax = demo.presentation.plotPopulationOnMap(
    census,
    tilesToolkit=tiles,
    density=True,                  # True=density, False=absolute counts
    density_units=ureg.km**2,      # people per km²
    zoomlevel=14,                  # tile detail level
    alpha=0.5,                     # see map underneath
    cmap="YlOrRd",
    outputCRS=ITM
)

# Absolute count overlay
ax = demo.presentation.plotPopulationOnMap(
    census,
    tilesToolkit=tiles,
    density=False,                 # absolute population counts
    alpha=0.4,
    cmap="Blues",
    outputCRS=ITM
)
```

**CRS constants:**

| Constant | Value | Description |
|----------|-------|-------------|
| `ITM` | 2039 | Israeli Transverse Mercator (meters) — default for plotting |
| `WSG84` | 4326 | WGS84 (degrees lat/lon) |

For the full API, see the [API Reference](../../developer_guide/api/measurements.md).
