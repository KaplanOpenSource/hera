# Understanding Toolkits

This page explains what toolkits are, how they work, and how to access them — without diving into specific API calls. For hands-on usage, see the [Toolkit Catalog](overview.md) and the individual toolkit pages.

---

## What is a toolkit?

A toolkit is a specialized module that provides three layers for working with one type of scientific data: a **data layer**, an **analysis layer**, and a **presentation layer**.

Without toolkits, you'd work directly with the Project API — adding documents, remembering file paths, writing your own analysis code, and building plots from scratch. A toolkit wraps all of that into a domain-specific interface:

| Layer | What it does |
|-------|-------------|
| **Data layer** | Normalizes the data: loading from various formats, converting between representations, and storing persistently in the database. You ask for data by name (`"YAVNEEL"`) instead of remembering file paths, formats, and query filters. |
| **Analysis** | Provides domain-specific processing methods. A meteorology toolkit knows how to compute hourly distributions; a GIS toolkit knows how to generate STL meshes. |
| **Presentation** | Generates domain-specific visualizations. Seasonal wind roses, elevation contour maps, casualty plots — each toolkit provides the plots its domain needs. |

Every toolkit inherits from `Project`, so it has full access to the project's data layer (measurements, simulations, cache). The toolkit adds its domain-specific layers on top.

---

## Data sources — your external databases

Toolkits organize their data through **data sources**. A data source is an external dataset — a weather station file, a DEM raster, a building shapefile — that has been registered with a name and a version.

Think of data sources as the toolkit's **external databases**:

- Each data source has a **name** (e.g., `"Israel_DEM_30m"`, `"YAVNEEL"`, `"BNTL_Buildings"`)
- Each data source has a **version** (e.g., `(1, 0, 0)`) so you can track updates
- The actual data stays on disk — the toolkit stores metadata about where to find it and how to read it

When you load a [repository](../user_guide/concepts.md#repositories) into a project, all the data sources it defines are registered with the appropriate toolkits. After that, you access data by name:

```python
# The toolkit knows where the data is and how to read it
elevation = topo.getDataSourceData("Israel_DEM_30m")
```

You don't need to know the file path, the data format, or the MongoDB query — the toolkit handles all of that. For the full details on versions and defaults, see [Working with data sources](../user_guide/working_with_data.md#data-sources).

---

## Accessing toolkits with ToolkitHome

You never create toolkit instances directly. Instead, you use **ToolkitHome** — a central registry that knows all available toolkits and creates instances for you.

A global `toolkitHome` singleton is created when you import Hera:

```python
from hera import toolkitHome
```

To get a toolkit, call `getToolkit` with the toolkit constant and a project name:

```python
# Get a topography toolkit bound to your project
topo = toolkitHome.getToolkit(toolkitHome.GIS_RASTER_TOPOGRAPHY, projectName="WindStudy")

# Get a meteorology toolkit
meteo = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_LOWFREQ, projectName="WindStudy")

# Get a risk assessment toolkit
risk = toolkitHome.getToolkit(toolkitHome.RISKASSESSMENT, projectName="WindStudy")
```

The `projectName` parameter is **always required** — either explicitly or via `caseConfiguration.json`. If you pass `projectName=None`, the toolkit loads the project name from `caseConfiguration.json` in the current working directory:

```python
# Explicit project name
meteo = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_LOWFREQ, projectName="WindStudy")

# projectName=None — loads from caseConfiguration.json in the current directory
meteo = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_LOWFREQ, projectName=None)

# Omitting projectName is equivalent to projectName=None
meteo = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_LOWFREQ)
```

!!! tip "Avoid hardcoding projectName"
    If you create your project directory with the CLI, you never need to specify `projectName` in your code:
    ```bash
    hera-project project create WindStudy --directory /data/wind_study
    cd /data/wind_study
    ```
    Now all Python scripts in that directory can omit `projectName`:
    ```python
    # No projectName needed — loaded from caseConfiguration.json
    meteo = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_LOWFREQ)
    ```
    This makes your scripts portable across different projects. See [Projects > Directory-based projects](../user_guide/projects.md#directory-based-projects) for details.

Each toolkit constant maps to a toolkit name string. Using constants avoids typos:

| Constant | String value |
|----------|-------------|
| `toolkitHome.GIS_RASTER_TOPOGRAPHY` | `"GIS_Raster_Topography"` |
| `toolkitHome.METEOROLOGY_LOWFREQ` | `"MeteoLowFreq"` |
| `toolkitHome.SIMULATIONS_OPENFOAM` | `"OpenFOAM"` |
| `toolkitHome.RISKASSESSMENT` | `"RiskAssessment"` |

See the [Toolkit Catalog](overview.md) for the full list of constants and toolkit names.

### Project binding

Every toolkit is **bound to a project** when you create it. This means:

- The toolkit reads data sources from that project
- Any documents the toolkit creates go into that project
- The toolkit's configuration is project-specific

The project name is resolved as follows:

1. If `projectName="WindStudy"` — use that name directly
2. If `projectName=None` (or omitted) — load from `caseConfiguration.json` in the current directory
3. If no `caseConfiguration.json` exists — use the read-only default project

### Listing available toolkits

```python
# Get a table of all available toolkits (built-in + registered)
toolkitHome.getToolkitTable(projectName="WindStudy")
```

From the CLI:

```bash
hera-toolkit list --project WindStudy
```

---

## Built-in vs dynamic toolkits

### Built-in toolkits

Hera ships with ~15 built-in toolkits covering GIS, meteorology, simulations, and risk assessment. These are hardcoded in the `ToolkitHome` registry and available immediately:

| Domain | Toolkits |
|--------|----------|
| GIS | Topography (raster/vector), Buildings, Demography, LandCover, Tiles, Shapes |
| Meteorology | MeteoLowFreq, MeteoHighFreq |
| Simulations | OpenFOAM, LSM, GaussianDispersion, WindProfile, hermesWorkflows |
| Risk | RiskAssessment |
| Data | experiment, ML/DL |

### Dynamic toolkits

You can also register **custom toolkits** at runtime. Dynamic toolkits are stored in the database and loaded on demand.

This is useful when:

- Your organization has domain-specific toolkits not included in Hera
- You want to share a custom toolkit across multiple projects
- You're developing a new toolkit and want to test it without modifying Hera's source code

#### Registering a dynamic toolkit

From the CLI:

```bash
# Register a toolkit by name and path
hera-project addToolkit myCustomToolkit /path/to/toolkit/directory
```

Or from Python:

```python
toolkitHome.registerToolkit(
    toolkit_name="myCustomToolkit",
    toolkit_path="/path/to/toolkit/directory",
    version=(1, 0, 0)
)
```

The toolkit directory should contain a Python module following the convention `<toolkit_name>/<toolkit_name>.py` with a class that inherits from `abstractToolkit`.

#### Using a dynamic toolkit

Once registered, dynamic toolkits are accessed the same way as built-in ones:

```python
# ToolkitHome searches: built-in registry → database → experiments
my_toolkit = toolkitHome.getToolkit("myCustomToolkit", projectName="MY_PROJECT")
```

#### How ToolkitHome resolves toolkit names

When you call `getToolkit(name)`, ToolkitHome searches three sources in order:

1. **Built-in registry** — the hardcoded dict of ~15 toolkits
2. **Database** — dynamically registered toolkits stored as `ToolkitDataSource` documents
3. **Experiments** — experiments exposed via the experiment toolkit

The first match wins. This means you can override a built-in toolkit by registering a dynamic one with the same name.

---

## The three layers

Every toolkit can have up to three layers, each accessed as a property:

### `toolkit.analysis`

The analysis layer provides domain-specific data processing. Examples:

- **Meteorology**: `addDatesColumns()`, `calcHourlyDist()`, `resampleSecondMoments()`
- **Risk assessment**: `getRiskAreas()` — compute risk zones from dispersion data
- **High-frequency meteo**: `calculateMeanData()`, turbulence statistics

### `toolkit.presentation`

The presentation layer generates domain-specific visualizations. Examples:

- **Meteorology**: `dailyPlots.plotScatter()`, `seasonalPlots.plotSeasonalHourly()`
- **Risk assessment**: `plotCasualtiesRose()` — radial casualty distribution
- **GIS Tiles**: `plot()` — render map tile images with CRS axes

### Data layer (data source management)

The data layer normalizes how data is loaded, converted, and stored. All toolkits inherit these methods from `abstractToolkit`:

- `getDataSourceData(name)` — load data by name (handles format detection and conversion)
- `getDataSourceList()` — list available data sources
- `addDataSource(name, resource, dataFormat)` — register a new data source with persistent storage
- `getDataSourceTable()` — view all data sources as a DataFrame

---

## Next steps

- **[Toolkit Catalog](overview.md)** — reference tables for all toolkits with methods and constants
- **[Measurements](measurements.md)** — GIS, meteorology, and experiment toolkit usage with code examples
- **[Simulations](simulations.md)** — OpenFOAM, LSM, Gaussian, wind profile usage
- **[Risk Assessment](risk_assessment.md)** — agents, effects, protection policies
- **[Working with Data Sources](../user_guide/working_with_data.md#data-sources)** — detailed API for versions, defaults, querying
