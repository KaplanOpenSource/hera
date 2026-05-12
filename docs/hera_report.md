# Hera: Full Platform Report

**Version:** 2.16.1  
**Type:** Scientific Data Management & Analysis Platform  
**Language:** Python  
**Source:** [github.com/KaplanOpenSource/hera](https://github.com/KaplanOpenSource/hera)

---

## Table of Contents

1. [What is Hera?](#1-what-is-hera)
2. [Core Architecture](#2-core-architecture)
3. [Data Layer](#3-data-layer)
4. [Toolkit System](#4-toolkit-system)
5. [GIS Toolkits](#5-gis-toolkits)
6. [Meteorology Toolkits](#6-meteorology-toolkits)
7. [Simulation Toolkits](#7-simulation-toolkits)
8. [Risk Assessment Toolkit](#8-risk-assessment-toolkit)
9. [Experiment & Data Toolkits](#9-experiment--data-toolkits)
10. [CLI Reference](#10-cli-reference)
11. [Best Practices](#11-best-practices)
12. [Project Structure](#12-project-structure)
13. [Quick Start](#13-quick-start)

---

## 1. What is Hera?

Hera is a Python-based platform for managing scientific data across measurements, simulations, and cached results. It provides a unified data layer backed by MongoDB and a rich set of domain-specific **Toolkits** covering GIS, meteorology, atmospheric dispersion modeling, CFD simulations, and risk assessment.

The platform is designed around a key principle: **domain-aware data management**. Rather than asking scientists to manage file paths, database queries, and data formats manually, Hera wraps all of that into clean, reusable toolkit interfaces that speak the language of each domain.

### Primary Use Cases

- Managing geospatial (GIS) data: terrain, buildings, population, land cover
- Ingesting and analysing meteorological station data (low- and high-frequency)
- Running and managing CFD simulations (OpenFOAM) and dispersion models (LSM, Gaussian)
- Conducting agent-based risk assessments from hazardous material releases
- Organizing experimental field campaign data with full lifecycle management

---

## 2. Core Architecture

The entire system is built around three core abstractions: **Project**, **ToolkitHome**, and **abstractToolkit**.

### 2.1 Project

`Project` is the central workspace. It represents a named container that groups all data (measurements, simulations, cached results) and configurations together. Every interaction with data goes through a `Project` instance, which manages three MongoDB collections and a local files directory.

```python
from hera import Project
proj = Project(projectName="MY_PROJECT")
```

Key methods on `Project`:

| Method | Description |
|--------|-------------|
| `addMeasurementsDocument(...)` | Store observational data |
| `addSimulationsDocument(...)` | Store simulation output |
| `addCacheDocument(...)` | Store intermediate results |
| `getDocuments(...)` | Query any collection |
| `setConfig(**kwargs)` | Write project-level key-value settings |
| `getConfig()` | Read project-level settings |
| `setCounter(name, defaultValue)` | Define an atomic sequential counter |
| `getCounterAndAdd(name, addition)` | Atomically increment a counter |

### 2.2 ToolkitHome

`ToolkitHome` is a singleton registry that knows all available toolkits and creates bound instances for you. You never instantiate toolkits directly.

```python
from hera import toolkitHome

# Retrieve a toolkit bound to a project
topo = toolkitHome.getToolkit(toolkitHome.GIS_RASTER_TOPOGRAPHY, projectName="MY_PROJECT")

# List all available toolkits
toolkitHome.getToolkitTable(projectName="MY_PROJECT")
```

From the CLI:
```bash
hera-toolkit list --project MY_PROJECT
```

When resolving a toolkit name, `ToolkitHome` searches three sources in order:
1. **Built-in registry** — the hardcoded set of ~15 toolkits
2. **Database** — dynamically registered toolkits stored as `ToolkitDataSource` documents
3. **Experiments** — experiments exposed via the experiment toolkit

The first match wins, so you can override a built-in by registering a dynamic toolkit with the same name.

### 2.3 abstractToolkit

`abstractToolkit` is the base class for all toolkits. It inherits from `Project`, so every toolkit has full data layer access plus domain-specific analysis and presentation capabilities.

Every toolkit exposes three layers:

| Layer | Property | Description |
|-------|----------|-------------|
| Data | (inherited) | Datasource management: load, register, version |
| Analysis | `toolkit.analysis` | Domain-specific data processing and statistics |
| Presentation | `toolkit.presentation` | Domain-specific plots and visualizations |

---

## 3. Data Layer

### 3.1 MongoDB Collections

Hera uses three MongoDB collections, each managed by a corresponding class:

| Collection | Class | Purpose |
|------------|-------|---------|
| **Measurements** | `Measurements_Collection` | Observational data, toolkit datasources |
| **Simulations** | `Simulations_Collection` | Simulation model outputs |
| **Cache** | `Cache_Collection` | Intermediate results, configurations |

All collections provide `addDocument()`, `getDocuments()`, and `deleteDocuments()` methods.

### 3.2 Documents

A Document is a MongoDB record representing a piece of data. Every document has:

| Field | Description |
|-------|-------------|
| `projectName` | The project it belongs to |
| `_cls` | Type discriminator: `Metadata.Measurements`, `Metadata.Simulations`, or `Metadata.Cache` |
| `type` | Application-defined type tag (e.g., `"ToolkitDataSource"`) |
| `resource` | Path to the actual data file or inline value |
| `dataFormat` | How to interpret the resource (e.g., `"parquet"`, `"JSON_dict"`) |
| `desc` | Free-form metadata dictionary |

### 3.3 DataSources

A DataSource is a registered, named, versioned reference to external data within a toolkit:

```python
# Register a datasource
toolkit.addDataSource("YAVNEEL", "/data/YAVNEEL.parquet", "parquet", version=[0, 0, 1])

# Retrieve data by name — no need to remember paths or formats
df = toolkit.getDataSourceData("YAVNEEL")

# List all datasources
toolkit.getDataSourceList()

# View as a DataFrame
toolkit.getDataSourceTable()
```

Datasources support semantic versioning `[major, minor, patch]`:

- **Major** `[X, 0, 0]` — Breaking changes, incompatible formats
- **Minor** `[0, X, 0]` — New fields, backward-compatible additions
- **Patch** `[0, 0, X]` — Bug fixes, minor corrections

Always set a default version for production use:
```python
toolkit.setDataSourceDefaultVersion("YAVNEEL", [0, 1, 0])
```

### 3.4 Repositories

A Repository is a JSON file that declares a collection of datasources, configurations, and documents organized by toolkit name. It acts as a blueprint for populating a project with data:

```json
{
    "MeteoLowFreq": {
        "Config": { "stationType": "IMS" },
        "Datasource": {
            "YAVNEEL": {
                "isRelativePath": "True",
                "item": {
                    "resource": "measurements/meteorology/YAVNEEL.parquet",
                    "dataFormat": "parquet"
                }
            }
        }
    }
}
```

When a repository is loaded into a project, all declared datasources are registered automatically.

### 3.5 Config and Counters

Projects store key-value configuration in a special Cache document:

```python
proj.setConfig(defaultSRTM="SRTMGL1", defaultCRS=4326)
config = proj.getConfig()
print(config["defaultSRTM"])  # "SRTMGL1"
```

Atomic counters are available for generating sequential IDs:

```python
proj.setCounter("experimentID", defaultValue=0)
current_id = proj.getCounterAndAdd("experimentID", addition=1)
```

---

## 4. Toolkit System

### 4.1 Built-in Toolkit Registry

| Constant | Toolkit Name | Category | Class |
|----------|-------------|----------|-------|
| `toolkitHome.GIS_RASTER_TOPOGRAPHY` | `GIS_Raster_Topography` | GIS | `TopographyToolkit` (raster) |
| `toolkitHome.GIS_VECTOR_TOPOGRAPHY` | `GIS_Vector_Topography` | GIS | `TopographyToolkit` (vector) |
| `toolkitHome.GIS_BUILDINGS` | `GIS_Buildings` | GIS | `BuildingsToolkit` |
| `toolkitHome.GIS_DEMOGRAPHY` | `GIS_Demography` | GIS | `DemographyToolkit` |
| `toolkitHome.GIS_LANDCOVER` | `GIS_LandCover` | GIS | `LandCoverToolkit` |
| `toolkitHome.GIS_TILES` | `GIS_Tiles` | GIS | `TilesToolkit` |
| `toolkitHome.METEOROLOGY_LOWFREQ` | `MeteoLowFreq` | Meteorology | `lowFreqToolKit` |
| `toolkitHome.METEOROLOGY_HIGHFREQ` | `MeteoHighFreq` | Meteorology | `HighFreqToolKit` |
| `toolkitHome.SIMULATIONS_OPENFOAM` | `OpenFOAM` | Simulations | `OFToolkit` |
| `toolkitHome.LSM` | `LSM` | Simulations | `LSMToolkit` |
| `toolkitHome.GAUSSIANDISPERSION` | `GaussianDispersion` | Simulations | `gaussianToolkit` |
| `toolkitHome.WINDPROFILE` | `WindProfile` | Simulations | `WindProfileToolkit` |
| `toolkitHome.RISKASSESSMENT` | `RiskAssessment` | Risk | `RiskToolkit` |
| `toolkitHome.EXPERIMENT` | `experiment` | Data | `experimentHome` |

Always use the constant rather than the string to avoid typos:

```python
# Recommended
topo = toolkitHome.getToolkit(toolkitHome.GIS_RASTER_TOPOGRAPHY, projectName="MY_PROJECT")
```

### 4.2 Dynamic (Custom) Toolkits

You can register custom toolkits at runtime:

```python
toolkitHome.registerToolkit(
    toolkit_name="myCustomToolkit",
    toolkit_path="/path/to/toolkit/directory",
    version=(1, 0, 0)
)
```

Or from the CLI:
```bash
hera-project addToolkit myCustomToolkit /path/to/toolkit/directory
```

The toolkit directory should contain a Python module following the convention `<toolkit_name>/<toolkit_name>.py` with a class that inherits from `abstractToolkit`.

---

## 5. GIS Toolkits

All GIS toolkits inherit from `abstractToolkit`. Vector-based toolkits share a common `VectorToolkit` base class.

### 5.1 GIS_Raster_Topography

Elevation data access and terrain analysis from SRTM (HGT) files.

- Load and query elevation grids
- Generate STL meshes for CFD simulation domains
- Terrain analysis (slope, aspect, viewshed)

### 5.2 GIS_Vector_Topography

Vector-based topography: contour lines and survey point data.

### 5.3 GIS_Buildings

Building footprint management and 3D geometry generation.

Key analysis methods:

| Method | Description |
|--------|-------------|
| `LambdaFromBuildingData(windDirection, resolution, ...)` | Compute block-averaged λp (plan area fraction), λf (frontal area density), hc (mean height) |
| `ConvexPolygons(regionNameOrData, buffer)` | Group nearby buildings into convex hulls |

Height resolution priority:
1. `BLDG_HT` column if present and > 0
2. `HI_PNT_Z - HT_LAND` if available
3. `building:levels * 3` as fallback
4. Buildings with `FTYPE` 14 or 16 are excluded from morphology calculations

### 5.4 GIS_Demography

Census population data management with spatial intersection capabilities.

Key methods:

| Method | Description |
|--------|-------------|
| `calculatePopulationInPolygon()` | Area-weighted spatial intersection |
| `createNewArea()` | Aggregate population within a custom polygon |

Presentation layer provides: density maps, population counts, age group distributions, polygon intersection plots, area annotation, and map overlays.

### 5.5 GIS_LandCover

MODIS MCD12Q1 land cover classification and surface roughness estimation.

### 5.6 GIS_Tiles

Tile-based raster data management. Renders map tile images with WGS84/ITM coordinate axes.

### 5.7 Coordinate Utilities

Available from `hera.measurements.GIS`:

```python
from hera.measurements.GIS import WSG84, ITM, convertCRS

# WSG84 = EPSG:4326 (WGS84)
# ITM   = EPSG:2039 (Israeli Transverse Mercator)

converted = convertCRS(points, inputCRS=WSG84, outputCRS=ITM)
```

---

## 6. Meteorology Toolkits

### 6.1 MeteoLowFreq

Manages hourly and daily meteorological station data (e.g., IMS stations). Supports Parquet-backed storage and rich statistical and visual analysis.

Analysis methods include:

- `addDatesColumns()` — Parse and attach date/time columns
- `calcHourlyDist()` — Compute hourly distributions
- `resampleSecondMoments()` — Statistical resampling

Presentation layer:

- `dailyPlots.plotScatter()` — Daily scatter plots
- `seasonalPlots.plotSeasonalHourly()` — Seasonal hourly heatmaps and wind roses

### 6.2 MeteoHighFreq

Manages high-frequency data from sonic anemometers and TRH sensors (typically at 10–20 Hz).

Supported file formats: Campbell Scientific binary, TOA5 CSV.

Analysis layer:

- `calculateMeanData()` — Compute block-averaged mean statistics
- Turbulence statistics: TKE, friction velocity, stability parameters

---

## 7. Simulation Toolkits

### 7.1 OpenFOAM (OFToolkit)

Manages the full lifecycle of OpenFOAM CFD simulations: templates, mesh generation, solver execution, and post-processing.

The toolkit uses **composition over inheritance** — solver-specific functionality (simpleFoam, buoyantReactingFoam, stochastic Lagrangian) is provided via extension objects rather than subclasses.

**Lifecycle phases:**

| Phase | Steps |
|-------|-------|
| Setup | Template management, case setup, mesh generation (blockMesh, snappyHexMesh) |
| Execution | Run simulation, parallel execution (MPI), monitor convergence |
| Post-Processing | Extract fields, generate plots, export VTK |

**Key usage:**

```python
of = toolkitHome.getToolkit(toolkitHome.SIMULATIONS_OPENFOAM, projectName="MY_PROJECT")

# Save a case template to the database
of.simpleFoam.saveTemplate("myTemplate", "/path/to/case")

# Load and run
of.simpleFoam.loadTemplate("myTemplate", toDirectory="/path/to/output", caseName="run_001")
of.runOFSimulation("wind_study_simpleFoam")
```

**Batch Slurm execution:**

```python
of.prepareSlurmWorkflowExecution(
    workflowName="simpleFoam_sweep",
    variations={"windSpeed": [3, 5, 8, 12]},
    jobName="wind_sweep"
)
```

**CLI:**

```bash
hera-openFoam simpleFoam templates list --projectName MY_PROJECT
hera-openFoam simpleFoam templates save myTemplate --projectName MY_PROJECT --directory /path
hera-openFoam simpleFoam templates load myTemplate --projectName MY_PROJECT --toDirectory /out --caseName run_001
```

### 7.2 LSM (Lagrangian Stochastic Model)

Manages atmospheric dispersion simulations using Lagrangian particle tracking. Handles the full lifecycle: templates, simulation runs (including batch Slurm execution), and results analysis (concentration fields, dosage calculations).

**Initialization:**

```python
lsm = toolkitHome.getToolkit(
    toolkitHome.LSM,
    projectName="MY_PROJECT",
    to_xarray=True,       # save results as xarray
    to_database=False,    # store runs in database
    forceKeep=False       # keep raw Lagrangian files
)
```

**Running simulations:**

```python
template.run(
    topography=topography_data,
    stations=weather_stations,
    simulationName="run_001",
    windSpeed=5.0,
    windDirection=270,
    releaseRate=1.0
)
```

**Querying results:**

```python
sim = lsm.getSimulations(windSpeed=5.0)[0]

concentration = sim.getConcentration(Q=1*kg)  # xarray Dataset with 'C' field
dosage = sim.getDosage(Q=1*kg)               # xarray Dataset with 'Dosage' field
```

**Integration with Risk Assessment:**

```python
concentration = sim.getConcentration(Q=1e6*mg)
risk = toolkitHome.getToolkit(toolkitHome.RISKASSESSMENT, projectName="MY_PROJECT")
agent = risk.getAgent("Chlorine")
toxic_loads = agent["RegularPopulation"].calculateToxicLoads(concentration, field="C")
```

**Slurm batch runs:**

```python
lsm.prepareSlurmLSMExecution(
    baseParameters={"windSpeed": 5.0, "releaseRate": 1.0},
    jsonVariations={"windDirection": [0, 90, 180, 270]},
    templateName="urban_dispersion",
    stations=stations_df,
    topography="/path/to/topography",
    jobName="dispersion_sweep"
)
```

### 7.3 GaussianDispersion

Gaussian puff and plume models for atmospheric pollutant transport.

- Cloud dispersion modeling
- Wind field and stability class integration
- Downwind concentration calculations

```python
gauss = toolkitHome.getToolkit(toolkitHome.GAUSSIANDISPERSION, projectName="MY_PROJECT")
```

### 7.4 WindProfile

Vertical wind profile modeling and analysis — used to characterize the atmospheric boundary layer structure for input to dispersion and CFD models.

```python
wp = toolkitHome.getToolkit(toolkitHome.WINDPROFILE, projectName="MY_PROJECT")
```

### 7.5 Hermes Workflow Toolkit

`hermesWorkflowToolkit` is the database-backed workflow orchestrator wrapping the Hermes library. It manages the full lifecycle of simulation workflows: creating them from JSON, storing them in MongoDB, building them into Luigi task DAGs, and executing them.

`OFToolkit` inherits from `hermesWorkflowToolkit`, gaining workflow support automatically.

**Class hierarchy:**

```
abstractToolkit
    └── hermesWorkflowToolkit
            ├── OFToolkit
            └── workflowToolkit (LSM variant)
```

---

## 8. Risk Assessment Toolkit

### 8.1 Overview

The `RiskAssessment` toolkit provides an agent-based risk modeling framework for hazardous material release scenarios. It integrates with GIS (population, buildings), meteorology (wind), and simulation (dispersion) toolkits.

**Inputs:**

- Population data (from `DemographyToolkit`)
- Hazard source (release scenario from LSM or Gaussian)
- Meteorology (from `MeteoLowFreq`)
- Building data (from `BuildingsToolkit`)

**Outputs:**

- Casualty estimates per severity level
- Risk maps
- Statistical analysis

### 8.2 Key Components

| Component | Description |
|-----------|-------------|
| **Agents** | Agent-based population model with spatial distribution |
| **Effects** | Injury level calculators (thermal, toxic, overpressure) |
| **Protection Policies** | Sheltering, evacuation, response models |
| **Analysis** | Casualty statistics, risk quantification |
| **Presentation** | Risk maps, casualty roses, bar plots |

### 8.3 Data Flow

```
InjuryLevel.calculateContours()
    → thresholdGeoDataFrame (contours at origin)
    → shiftLocationAndAngle(loc, angle)
    → thresholdGeoDataFrame (positioned on map)
    → project(demographic, loc, angle)
    → DataFrame (casualties per severity per time step)
    → groupby("severity").sum()
    → Total casualties per severity level
```

### 8.4 Analysis Methods

| Method | Description |
|--------|-------------|
| `getRiskAreas()` | Compute risk zones from dispersion data |
| `project(demographic, loc, angle)` | Project casualties onto population data at a given angle |

### 8.5 Presentation Methods

| Method | Description |
|--------|-------------|
| `plotCasualtiesRose()` | Radial casualty distribution across wind directions |
| Risk maps | Geographic risk distribution plots |

---

## 9. Experiment & Data Toolkits

### 9.1 experiment (experimentHome)

Manages experimental field campaign data with full lifecycle management — organizing raw data files into structured experiments with device metadata, trial sets, and transmission quality monitoring.

**Data engines (interchangeable):**

| Engine | Backend | Returns | Best for |
|--------|---------|---------|----------|
| `parquetDataEngineHera` | Hera data layer (Parquet) | Dask/pandas DataFrame | Local file-based experiments |
| `pandasDataEngineDB` | MongoDB direct | pandas DataFrame | Small-to-medium MongoDB datasets |
| `daskDataEngineDB` | MongoDB via Dask | Dask DataFrame | Large datasets requiring lazy evaluation |

**Analysis methods:**

| Method | Description |
|--------|-------------|
| `getDeviceLocations(...)` | Device location metadata as DataFrame |
| `getTurbulenceStatistics(...)` | Turbulence analysis for sonic data |
| `getDeviceTypeTransmissionFrequencyOfTrial(...)` | Data transmission frequency heatmap |
| `addMetadata(...)` | Merge device metadata into a dataset |
| `addTrialProperties(...)` | Add time-delta columns relative to trial start/release |

**Presentation:**

- `plotDevicesOnImage()` — Device locations on a map image
- `plotDeviceTypeFunctionality()` — Heatmap of normalized transmission frequency (red=none, green=good)
- `generateLatexTable()` — LaTeX/PDF report with device maps and metadata

### 9.2 dataToolkit

The special repository management toolkit used to load repositories into projects. See the Data Layer section for full details.

---

## 10. CLI Reference

Hera ships with several CLI entry points, all prefixed with `hera-`:

| Command | Description |
|---------|-------------|
| `hera-project` | Project management: create, list, delete projects; load repositories |
| `hera-toolkit` | List available toolkits, register dynamic toolkits |
| `hera-GIS` | GIS operations from the command line |
| `hera-LSM` | Load LSM templates, list simulations |
| `hera-openFoam` | Manage OpenFOAM templates and workflow groups |
| `hera-workflows` | Manage simulation workflow groups |
| `hera-rag-search` | Start/manage the RAG (AI search) API server |

**Examples:**

```bash
# Create a project
hera-project project create MY_PROJECT

# Load a repository
hera-project repository load MY_PROJECT /path/to/repository.json

# List toolkits
hera-toolkit list --project MY_PROJECT

# Register a custom toolkit
hera-project addToolkit myCustomToolkit /path/to/toolkit

# List LSM simulations
hera-LSM list MY_PROJECT

# OpenFOAM template management
hera-openFoam simpleFoam templates list --projectName MY_PROJECT
```

---

## 11. Best Practices

### Project Organization

- Use uppercase, descriptive project names with underscores: `WIND_ANALYSIS_2024`
- Create separate projects for different domains, time periods, or experiments
- Keep related analyses in one project when toolkits share datasources

### Datasource Naming

- Use clear, consistent names: `YAVNEEL`, `SRTMGL1`, `lamas_population`
- Match file names when possible (without extension)
- Always set a default version for production use

### Versioning

```python
# Initial version
toolkit.addDataSource("YAVNEEL", "v1.parquet", "parquet", version=[0, 0, 1])

# Set default after loading
toolkit.setDataSourceDefaultVersion("YAVNEEL", [0, 1, 0])
```

### Common Pitfalls

**Not setting default versions:**
```python
# Bad: ambiguous when multiple versions exist
data = toolkit.getDataSourceData("YAVNEEL")

# Good: always explicit
toolkit.setDataSourceDefaultVersion("YAVNEEL", [0, 1, 0])
data = toolkit.getDataSourceData("YAVNEEL")
```

**Materializing large datasets unnecessarily:**
```python
# Bad: loads everything into memory
df = toolkit.getDataSourceData("HUGE").compute()
result = df.head(100)

# Good: lazy evaluation
df = toolkit.getDataSourceData("HUGE")
result = df.head(100).compute()
```

**Mixing collection types:**
```python
# Bad
proj.addMeasurementsDocument(..., type="SimulationResult")

# Good
proj.addSimulationsDocument(..., type="SimulationResult")
```

**Using absolute paths in repositories:**
```json
// Bad: not portable
{ "resource": "/absolute/path/to/data.parquet" }

// Good: portable
{ "isRelativePath": "True", "item": { "resource": "data/data.parquet" } }
```

---

## 12. Project Structure

```
hera/
├── hera/                              # Main package
│   ├── __init__.py                    # Version, logging, toolkitHome singleton
│   ├── toolkit.py                     # ToolkitHome + abstractToolkit
│   ├── datalayer/                     # Project, Collection, DataHandler, datatypes
│   ├── measurements/                  # GIS, meteorology, experiment toolkits
│   │   ├── GIS/
│   │   │   ├── raster/                # TopographyToolkit, LandCoverToolkit, TilesToolkit
│   │   │   └── vector/                # VectorToolkit, BuildingsToolkit, DemographyToolkit
│   │   ├── meteorology/
│   │   │   ├── lowfreqdata/           # lowFreqToolKit
│   │   │   └── highfreqdata/          # HighFreqToolKit
│   │   └── experiment/                # experimentHome
│   ├── simulations/                   # OpenFOAM, LSM, Gaussian, WindProfile
│   │   ├── openFoam/
│   │   ├── LSM/
│   │   ├── gaussian/
│   │   ├── windProfile/
│   │   └── hermesWorkflowToolkit.py
│   ├── riskassessment/                # Risk assessment agents and toolkit
│   ├── utils/                         # Logging, unit handling, data utilities
│   └── bin/                           # CLI entry points
├── docs/                              # MkDocs documentation source
├── notebooks/                         # Jupyter exploration notebooks
├── pytest.ini
├── setup.py
├── requirements.txt
└── mkdocs.yml
```

---

## 13. Quick Start

### Installation

```bash
git clone https://github.com/KaplanOpenSource/hera
cd hera
source init_with_mongo.sh
```

### Basic Usage

```python
from hera import Project, toolkitHome

# Create or open a project
proj = Project(projectName="MY_PROJECT")

# Get a topography toolkit
topo = toolkitHome.getToolkit(toolkitHome.GIS_RASTER_TOPOGRAPHY, projectName="MY_PROJECT")

# Register a data source
topo.addDataSource("SRTMGL1", "/data/terrain.hgt", "HGT", version=[0, 0, 1])
topo.setDataSourceDefaultVersion("SRTMGL1", [0, 0, 1])

# Load data
elevation = topo.getDataSourceData("SRTMGL1")

# Get a meteorology toolkit
meteo = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_LOWFREQ, projectName="MY_PROJECT")

# Run a dispersion simulation
lsm = toolkitHome.getToolkit(toolkitHome.LSM, projectName="MY_PROJECT")
sim = lsm.getSimulations(windSpeed=5.0)[0]
concentration = sim.getConcentration(Q=1)  # kg

# Risk assessment
risk = toolkitHome.getToolkit(toolkitHome.RISKASSESSMENT, projectName="MY_PROJECT")
```

### Loading a Repository

```python
from hera import toolkitHome

# Load all data sources declared in a repository JSON
toolkitHome.getToolkit("heraData", projectName="MY_PROJECT").loadRepository(
    "/path/to/repository.json"
)
```

Or from the CLI:
```bash
hera-project repository load MY_PROJECT /path/to/repository.json
```

---

*Report generated from Hera v2.16.1 documentation.*
