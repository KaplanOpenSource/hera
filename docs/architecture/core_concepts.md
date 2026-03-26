# Core Concepts

This page provides a deep technical reference for the three foundational classes in Hera: **Project**, **ToolkitHome**, and **abstractToolkit**. Understanding these is essential for working with any part of the system.

---

## Class Hierarchy Overview

The diagram below shows the **core hierarchy** — `Project` at the base, `abstractToolkit` adding datasource management on top, and `ToolkitHome` acting as the registry/factory.

![Diagram](../images/diagrams/architecture_core_concepts_0_666f80ef.svg)

<!-- mermaid source (for editing, paste into mermaid.live):
```mermaid
classDiagram
    class Project {
        +str projectName
        +str FilesDirectory
        +Measurements_Collection measurements
        +Simulations_Collection simulations
        +Cache_Collection cache
        +addMeasurementsDocument(resource, dataFormat, type, desc)
        +getMeasurementsDocuments(**filters) list
        +deleteMeasurementsDocuments(**filters)
        +addSimulationsDocument(resource, dataFormat, type, desc)
        +getSimulationsDocuments(**filters) list
        +addCacheDocument(resource, dataFormat, type, desc)
        +getCacheDocuments(**filters) list
        +setConfig(**kwargs)
        +getConfig() dict
        +setCounter(counterName, defaultValue)
        +getCounterAndAdd(counterName, addition) int
        +saveData(name, data, desc, kind, type, dataFormat)
        +export(path)
        +load(proj, path, is_hard_import)
    }

    class abstractToolkit {
        +str _toolkitname
        +str _projectName
        +object _analysis
        +object _presentation
        +str toolkitName
        +str projectName
        +object analysis
        +object presentation
        +addDataSource(name, resource, dataFormat, version, overwrite)
        +getDataSourceDocument(name, version, **filters)
        +getDataSourceData(name, version, **filters)
        +getDataSourceList(**filters) list
        +getDataSourceMap(**filters) list
        +getDataSourceTable(**filters) DataFrame
        +deleteDataSource(name, version, **filters)
        +setDataSourceDefaultVersion(name, version)
    }

    class ToolkitHome {
        +dict _toolkits
        +getToolkit(toolkitName, projectName, filesDirectory) abstractToolkit
        +getToolkitTable(projectName) DataFrame
        +registerToolkit(toolkitclass, datasource_name, params, version, overwrite, projectName, repositoryName)
        +auto_register_and_get(toolkitName, projectName, repositoryJSON, repositoryName, params, version) abstractToolkit
        +getToolkitDocuments(name, projectName) list
        +setDefaultRepository(projectName, repositoryName)
        +getDefaultRepository(projectName) str
        +getDatasourceDocument(projectName, datasourceName, repositoryName, version)
        +import_toolkits_from_json(projectName, json_path, default_repository, overwrite) list
        +import_experiments_from_json(projectName, json_path) list
    }

    Project <|-- abstractToolkit : "inherits data layer"
    ToolkitHome ..> abstractToolkit : "instantiates via getToolkit()"
```
-->    Project <|-- abstractToolkit : "inherits data layer"
    ToolkitHome ..> abstractToolkit : "instantiates via getToolkit()"
```
-->
-->    Project <|-- abstractToolkit : "inherits data layer"
    ToolkitHome ..> abstractToolkit : "instantiates via getToolkit()"
```
-->
-->

| Class | Role | Key methods |
|-------|------|------------|
| `Project` | Central data access layer — CRUD for measurements, simulations, cache | `addMeasurementsDocument`, `getMeasurementsDocuments`, `setConfig`, `getCounterAndAdd`, `saveData`, `export` |
| `abstractToolkit` | Base for all domain toolkits — adds datasource management, analysis/presentation layers | `addDataSource`, `getDataSourceData`, `getDataSourceList`, `deleteDataSource`, `setDataSourceDefaultVersion` |
| `ToolkitHome` | Registry/factory — resolves toolkit names to classes | `getToolkit`, `getToolkitTable`, `registerToolkit`, `getToolkitDocuments` |

### Concrete Toolkit Hierarchy

All domain toolkits extend `abstractToolkit`. The diagram below shows the full inheritance tree including the special `dataToolkit` that manages repositories.

![Diagram](../images/diagrams/architecture_core_concepts_1_4e97226d.svg)

<!-- mermaid source (for editing, paste into mermaid.live):
```mermaid
classDiagram
    class abstractToolkit {
        +str toolkitName
        +object analysis
        +object presentation
        +addDataSource()
        +getDataSourceData()
    }

    class dataToolkit {
        +addRepository(repositoryName, repositoryPath)
        +getRepository(repositoryName) dict
        +getRepositoryTable() DataFrame
        +loadAllDatasourcesInRepositoryToProject(projectName, repositoryName, overwrite)
        +loadAllDatasourcesInRepositoryJSONToProject(projectName, repositoryJSON, basedir, overwrite, auto_register_missing)
        +resolveDataSourcePaths(repositoryJSON, basedir) dict
        +loadRepositoryFromPath(json_path) dict
    }

    class TopographyToolkit {
        +getPointElevation(lat, lon) float
        +getPointListElevation(points) list
        +getElevation(minx, miny, maxx, maxy, dxdy) xarray
        +getElevationOfXarray(ds) xarray
        +createElevationSTL(...)
        +convertPointsCRS(points, inputCRS, outputCRS) list
    }

    class lowFreqToolKit {
        +analysis LowFreqAnalysis
        +presentation LowFreqPresentation
        +str docType
    }

    class HighFreqToolKit {
        +analysis HighFreqAnalysis
        +str docType
    }

    class DemographyToolkit {
        +calculatePopulationInPolygon(polygon, datasource) float
        +createNewArea(polygon, datasource) GeoDataFrame
        +setDefaultDirectory(path)
    }

    class LandCoverToolkit {
        +getLandCoverAtPoint(lat, lon) str
        +getLandCover(minx, miny, maxx, maxy, dxdy) xarray
        +getRoughness(...) xarray
    }

    class BuildingsToolkit {
        +createBuildingsSTL(...)
        +getBuildingsGDF(...)
    }

    class OFToolkit {
        +runCase(casePath, nproc)
        +analysis OFAnalysis
    }

    class LSMToolkit {
        +analysis LSMAnalysis
        +presentation LSMPresentation
    }

    class RiskToolkit {
        +analysis RiskAnalysis
        +presentation RiskPresentation
    }

    abstractToolkit <|-- dataToolkit : "repository management"
    abstractToolkit <|-- TopographyToolkit : "GIS raster"
    abstractToolkit <|-- lowFreqToolKit : "meteorology"
    abstractToolkit <|-- HighFreqToolKit : "meteorology"
    abstractToolkit <|-- DemographyToolkit : "GIS vector"
    abstractToolkit <|-- LandCoverToolkit : "GIS raster"
    abstractToolkit <|-- BuildingsToolkit : "GIS vector"
    abstractToolkit <|-- OFToolkit : "simulations"
    abstractToolkit <|-- LSMToolkit : "simulations"
    abstractToolkit <|-- RiskToolkit : "risk assessment"
```
-->t : "simulations"
    abstractToolkit <|-- LSMToolkit : "simulations"
    abstractToolkit <|-- RiskToolkit : "risk assessment"
```
-->
-->t : "simulations"
    abstractToolkit <|-- LSMToolkit : "simulations"
    abstractToolkit <|-- RiskToolkit : "risk assessment"
```
-->
-->

| Toolkit class | Domain | Key capabilities |
|--------------|--------|-----------------|
| `dataToolkit` | Data | Repository management, loading data sources into projects |
| `TopographyToolkit` | GIS (raster) | Elevation data, STL generation, CRS conversion |
| `BuildingsToolkit` | GIS (vector) | Building footprints, 3D STL for CFD |
| `DemographyToolkit` | GIS (vector) | Population calculations within polygons |
| `LandCoverToolkit` | GIS (raster) | Land cover classification, roughness estimation |
| `lowFreqToolKit` | Meteorology | Hourly/daily station data with analysis + presentation |
| `HighFreqToolKit` | Meteorology | High-frequency sonic anemometer data |
| `OFToolkit` | Simulations | OpenFOAM lifecycle (templates, running, post-processing) |
| `LSMToolkit` | Simulations | Lagrangian Stochastic Model |
| `RiskToolkit` | Risk | Agent-based hazard and casualty modeling |

---

## Project

**Source:** `hera/datalayer/project.py`

The `Project` class is the central data access layer. Every interaction with stored data — measurements, simulations, cached results — goes through a `Project` instance.

### Responsibilities

1. **Document CRUD** — Add, query, and delete documents across three MongoDB collections: Measurements, Simulations, and Cache.
2. **Configuration** — Per-project key-value config stored as a special Cache document (type = `<projectName>__config__`).
3. **Counters** — Atomic counters for generating sequential IDs (stored inside the config document under `desc.counters`).
4. **File Management** — Each project has a `filesDirectory` where physical files (parquet, netcdf, etc.) are stored on disk.
5. **Import/Export** — Projects can be exported to a zip file and loaded into another MongoDB instance.

### Initialization Flow

![Diagram](../images/diagrams/architecture_core_concepts_2_003f2fd3.svg)

<!-- mermaid source (for editing, paste into mermaid.live):
```mermaid
flowchart TD
    Start["Project(projectName=...,\nconnectionName=...,\nfilesDirectory=...)"] --> CheckName{projectName\nis None?}

    CheckName -- "Yes" --> LoadCase["Load\ncaseConfiguration.json\nfrom current directory"]
    LoadCase --> FileExists{File\nexists?}
    FileExists -- "Yes" --> ExtractName["Extract projectName\nfrom JSON file"]
    FileExists -- "No" --> UseDefault["Use DEFAULTPROJECT\n(read-only shared project)"]

    CheckName -- "No" --> UseProvided["Use provided\nprojectName as-is"]

    ExtractName --> InitCollections
    UseDefault --> InitCollections
    UseProvided --> InitCollections

    subgraph InitCollections ["Initialize MongoDB Collections"]
        direction LR
        CreateMeas["Create\nMeasurements_Collection\n(measurements)"]
        CreateSim["Create\nSimulations_Collection\n(simulations)"]
        CreateCache["Create\nCache_Collection\n(cache)"]
        CreateAbstract["Create\nAbstractCollection\n(general queries)"]
    end

    InitCollections --> ResolveFiles["Resolve filesDirectory\npath on disk"]
    ResolveFiles --> SaveConfig["Save filesDirectory\nto project config\n(if not already saved)"]
    SaveConfig --> Ready["Project instance\nready for use"]
```
-->["Save filesDirectory\nto project config\n(if not already saved)"]
    SaveConfig --> Ready["Project instance\nready for use"]
```
-->
-->["Save filesDirectory\nto project config\n(if not already saved)"]
    SaveConfig --> Ready["Project instance\nready for use"]
```
-->
-->

!!! note "The Default Project"
    When `projectName=None` and no `caseConfiguration.json` exists, Hera uses a special `"defaultProject"` which is **read-only**. This is used by the `dataToolkit` to store repository metadata without polluting user projects.

### Key Methods

| Method | Description |
|--------|-------------|
| `addMeasurementsDocument(resource, dataFormat, type, desc)` | Insert a new measurement document into MongoDB |
| `getMeasurementsDocuments(**filters)` | Query measurement documents by any combination of fields |
| `deleteMeasurementsDocuments(**kwargs)` | Remove matching measurement documents |
| `setConfig(**kwargs)` | Update project configuration (atomic per-key updates) |
| `getConfig()` | Return the full config dict |
| `setCounter(name, default)` | Define a named counter |
| `getCounterAndAdd(name, addition)` | Atomically read and increment a counter |
| `saveData(name, data, desc, kind)` | Auto-detect format, save to disk, and create a document |
| `export(path)` | Export all project documents to a zip file |
| `Project.load(proj, path, is_hard_import)` | Import documents from an exported zip |

### Document Structure

Every document in Hera has the following top-level fields:

```python
{
    "projectName": "MY_PROJECT",
    "_cls": "Metadata.Measurements",   # or .Simulations / .Cache
    "type": "ToolkitDataSource",        # application-defined type tag
    "resource": "/path/to/data.parquet",
    "dataFormat": "parquet",
    "desc": {                           # free-form metadata dict
        "toolkit": "MeteoLowFreq",
        "datasourceName": "YAVNEEL",
        "version": [0, 0, 1],
        ...
    }
}
```

---

## ToolkitHome

**Source:** `hera/toolkit.py` (class `ToolkitHome`)

`ToolkitHome` is the **central registry** for all available toolkits. A singleton instance `toolkitHome` is created at import time in `hera/__init__.py`.

### Responsibilities

1. **Static Registry** — A dict mapping toolkit names (e.g., `"MeteoLowFreq"`) to their full Python class paths.
2. **Dynamic Registry** — Toolkits registered at runtime via `registerToolkit()` are stored as `ToolkitDataSource` documents in MongoDB.
3. **Instantiation** — `getToolkit(name, projectName)` resolves the class, imports it, and returns a new instance.
4. **Auto-Registration** — `auto_register_and_get()` can discover toolkit classes from repository JSON hints or DB documents.

### Toolkit Resolution Flow

![Diagram](../images/diagrams/architecture_core_concepts_3_4ef816a5.svg)

<!-- mermaid source (for editing, paste into mermaid.live):
```mermaid
flowchart TD
    Start["getToolkit(toolkitName,\nprojectName, filesDirectory)"] --> CheckStatic{Is toolkitName\nin static\nregistry?}

    CheckStatic -- "Yes" --> GetClasspath["Get classpath from\nstatic _toolkits dict"]
    GetClasspath --> ImportClass["pydoc.locate(classpath)\nDynamically import class"]
    ImportClass --> Instantiate["Create instance:\nToolkitClass(toolkitName,\nprojectName, ...)"]
    Instantiate --> ReturnInstance["Return toolkit\ninstance"]

    CheckStatic -- "No" --> CheckDB{Query MongoDB\nfor ToolkitDataSource\ndocument?}
    CheckDB -- "Found" --> LoadFromDB["Extract classpath\nfrom document desc"]
    LoadFromDB --> ImportClass

    CheckDB -- "Not found" --> AutoRegister{auto_register_missing\nenabled?}
    AutoRegister -- "Yes" --> TryJSON["Search repository JSON\nfor classpath hints"]
    TryJSON --> RegisterAndImport["registerToolkit()\nthen import"]
    RegisterAndImport --> Instantiate

    AutoRegister -- "No" --> ReturnNone["Return None\n(toolkit not available)"]
```
-->mport"]
    RegisterAndImport --> Instantiate

    AutoRegister -- "No" --> ReturnNone["Return None\n(toolkit not available)"]
```
-->
-->mport"]
    RegisterAndImport --> Instantiate

    AutoRegister -- "No" --> ReturnNone["Return None\n(toolkit not available)"]
```
-->
-->

!!! warning "Singleton Pattern"
    `ToolkitHome` is instantiated once at `hera/__init__.py` as `toolkitHome = ToolkitHome()`. All code should use this singleton rather than creating new instances. The static registry dict is populated in `__init__` and shared across the entire process.

### Static Toolkit Registry

The `_toolkits` dict is structured as:

```python
{
    "GIS_Raster_Topography": {
        "cls": "hera.measurements.GIS.raster.topography.TopographyToolkit",
        "desc": None,
        "type": "measurements"
    },
    "MeteoLowFreq": {
        "cls": "hera.measurements.meteorology.lowfreqdata.toolkit.lowFreqToolKit",
        "desc": None,
        "type": "measurements"
    },
    # ... 16 built-in toolkits total
}
```

### Save Modes

`ToolkitHome` defines standard save-mode constants used throughout the system:

| Constant | Value | Meaning |
|----------|-------|---------|
| `TOOLKIT_SAVEMODE_NOSAVE` | `None` | Do not persist |
| `TOOLKIT_SAVEMODE_ONLYFILE` | `"File"` | Save to disk only |
| `TOOLKIT_SAVEMODE_ONLYFILE_REPLACE` | `"File_overwrite"` | Save to disk, overwrite |
| `TOOLKIT_SAVEMODE_FILEANDDB` | `"DB"` | Save to disk + MongoDB document |
| `TOOLKIT_SAVEMODE_FILEANDDB_REPLACE` | `"DB_overwrite"` | Save to disk + MongoDB, overwrite |

---

## abstractToolkit

**Source:** `hera/toolkit.py` (class `abstractToolkit`)

`abstractToolkit` is the **base class for all domain toolkits**. It inherits from `Project`, meaning every toolkit instance **is** a project — it has full access to measurements, simulations, cache, config, and counters.

!!! info "Design Decision: Inheritance over Composition"
    `abstractToolkit` inherits from `Project` rather than wrapping it. This means every toolkit call like `toolkit.getMeasurementsDocuments(...)` goes directly to the `Project` implementation. The toolkit adds a **datasource management layer** on top, where each datasource is a measurement document tagged with `type="ToolkitDataSource"` and the toolkit's name.

### What It Adds Over Project

1. **Toolkit Identity** — `toolkitName` property, automatically injected into every document created by the toolkit.
2. **Datasource Versioning** — Each datasource has a name + version tuple `(major, minor, patch)`. The system supports multiple versions and a default-version mechanism.
3. **Analysis & Presentation Layers** — `self._analysis` and `self._presentation` hold domain-specific logic (e.g., statistical calculations, plotting).
4. **Automatic Tagging** — `addMeasurementsDocument`, `addSimulationsDocument`, and `addCacheDocument` are overridden to inject the toolkit name into every document's `desc`.

### Datasource Lifecycle

![Diagram](../images/diagrams/architecture_core_concepts_4_099a0520.svg)

<!-- mermaid source (for editing, paste into mermaid.live):
```mermaid
flowchart LR
    subgraph Register ["1. Register Datasource"]
        AddDS["addDataSource(\nname, resource,\ndataFormat, version,\noverwrite)"]
        StoreMongo["Store as\nToolkitDataSource\ndocument in MongoDB\n\ntype = ToolkitDataSource\ndesc.toolkit = toolkitName\ndesc.datasourceName = name\ndesc.version = version"]
    end

    subgraph Retrieve ["2. Retrieve Datasource"]
        GetDoc["getDataSourceDocument(\nname, version)"]
        QueryDB["Query MongoDB:\ntype=ToolkitDataSource\ntoolkit=toolkitName\ndatasourceName=name"]
        VersionResolve{"Version\nspecified?"}
        UseVersion["Filter by\nexact version"]
        UseDefault["Use default\nversion"]
    end

    subgraph Load ["3. Load Data"]
        GetData["getDataSourceData(\nname, version)"]
        CallGetData["doc.getData()\nvia DataHandler"]
        ReturnData["Returns:\nDataFrame, xarray,\ndict, string, etc.\nbased on dataFormat"]
    end

    AddDS --> StoreMongo
    StoreMongo -.-> GetDoc
    GetDoc --> QueryDB
    QueryDB --> VersionResolve
    VersionResolve -- "Yes" --> UseVersion
    VersionResolve -- "No" --> UseDefault
    UseVersion --> GetData
    UseDefault --> GetData
    GetData --> CallGetData
    CallGetData --> ReturnData
```
-->-> UseDefault
    UseVersion --> GetData
    UseDefault --> GetData
    GetData --> CallGetData
    CallGetData --> ReturnData
```
-->
-->-> UseDefault
    UseVersion --> GetData
    UseDefault --> GetData
    GetData --> CallGetData
    CallGetData --> ReturnData
```
-->
-->

### Key Methods

| Method | Description |
|--------|-------------|
| `addDataSource(name, resource, dataFormat, version)` | Register a new datasource (versioned) |
| `getDataSourceDocument(name, version)` | Get the MongoDB document for a datasource |
| `getDataSourceData(name, version)` | Load and return the actual data |
| `getDataSourceList(**filters)` | List all datasource names |
| `getDataSourceTable(**filters)` | DataFrame of all datasources with metadata |
| `deleteDataSource(name, version)` | Remove a datasource document |
| `setDataSourceDefaultVersion(name, version)` | Set which version is returned by default |

### Concrete Toolkit Example

A concrete toolkit like `lowFreqToolKit` extends `abstractToolkit` and provides:

- **`analysis`** — Methods like `addDatesColumns()`, `calcHourlyDist()`, `resampleSecondMoments()`
- **`presentation`** — Plot generators like `dailyPlots.plotScatter()`, `seasonalPlots.plotProbContourf_bySeason()`
- **`docType`** — A string identifying the data type this toolkit handles

```python
from hera import toolkitHome

# Instantiate via the registry
lf = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_LOWFREQ, projectName="MY_PROJECT")

# Access data
df = lf.getDataSourceData("YAVNEEL")

# Use analysis
enriched = lf.analysis.addDatesColumns(df, datecolumn="datetime")

# Use presentation
ax = lf.presentation.dailyPlots.plotScatter(enriched, plotField="RH")
```

---

## dataToolkit

**Source:** `hera/utils/data/toolkit.py`

`dataToolkit` extends `abstractToolkit` and serves as the **repository management layer**. It operates on the `defaultProject` (read-only for normal users) and manages the mapping between repository JSON files and project datasources.

### Responsibilities

1. **Repository Registration** — `addRepository(name, path)` stores a reference to a repository JSON file as a datasource in the default project.
2. **Repository Loading** — `loadAllDatasourcesInRepositoryJSONToProject()` iterates through a parsed repository JSON and, for each toolkit section, dispatches to specialized handlers.
3. **Static Helpers** — `resolveDataSourcePaths()` and `loadRepositoryFromPath()` allow working with repository JSONs without MongoDB.

### Repository JSON Loading Pipeline

![Diagram](../images/diagrams/architecture_core_concepts_5_98ea1dbf.svg)

<!-- mermaid source (for editing, paste into mermaid.live):
```mermaid
flowchart TD
    Start["loadAllDatasourcesIn\nRepositoryJSONToProject(\nprojectName, repositoryJSON,\nbasedir, overwrite,\nauto_register_missing)"] --> ParseJSON["Parse repositoryJSON\ndictionary"]

    ParseJSON --> IterateToolkits["Iterate over each\ntoolkit name in JSON"]

    IterateToolkits --> GetToolkit["toolkitHome.getToolkit(\ntoolkitName, projectName)"]

    GetToolkit --> ToolkitFound{Toolkit\nfound?}

    ToolkitFound -- "Yes" --> ProcessSections["Process sections\nfor this toolkit"]
    ToolkitFound -- "No" --> AutoReg{auto_register_missing\nis True?}
    AutoReg -- "Yes" --> TryAutoReg["auto_register_and_get(\ntoolkitName, projectName,\nrepositoryJSON)"]
    TryAutoReg --> ProcessSections
    AutoReg -- "No" --> SkipToolkit["Skip this toolkit\n(log warning)"]
    SkipToolkit --> IterateToolkits

    ProcessSections --> SectionType{Section\ntype?}

    SectionType -- "Config" --> HandleConfig["_handle_Config\ntoolkit.setConfig(...)"]
    SectionType -- "DataSource" --> HandleDS["_handle_DataSource\ntoolkit.addDataSource(...)\nResolves relative paths\nusing basedir"]
    SectionType -- "Measurements" --> HandleMeas["_DocumentHandler\ntoolkit.addMeasurements\nDocument(...)"]
    SectionType -- "Simulations" --> HandleSim["_DocumentHandler\ntoolkit.addSimulations\nDocument(...)"]
    SectionType -- "Cache" --> HandleCache["_DocumentHandler\ntoolkit.addCache\nDocument(...)"]
    SectionType -- "Function" --> HandleFunc["_handle_Function\nCall named function\nwith parameters"]

    HandleConfig --> NextSection["Next section"]
    HandleDS --> NextSection
    HandleMeas --> NextSection
    HandleSim --> NextSection
    HandleCache --> NextSection
    HandleFunc --> NextSection
    NextSection --> SectionType
```
-->on
    HandleSim --> NextSection
    HandleCache --> NextSection
    HandleFunc --> NextSection
    NextSection --> SectionType
```
-->
-->on
    HandleSim --> NextSection
    HandleCache --> NextSection
    HandleFunc --> NextSection
    NextSection --> SectionType
```
-->
-->

### Handler Dispatch Table

| JSON Section Key | Handler Method | Action |
|-----------------|----------------|--------|
| `Config` | `_handle_Config` | Calls `toolkit.setConfig(**docTypeDict)` |
| `Datasource` | `_handle_DataSource` | Adds each item as a versioned datasource |
| `Measurements` | `_DocumentHandler` | Adds raw measurement documents |
| `Simulations` | `_DocumentHandler` | Adds simulation documents |
| `Cache` | `_DocumentHandler` | Adds cache documents |
| `Function` | `_handle_Function` | Calls a named function with parameters |

---

## API Reference

::: hera.toolkit.ToolkitHome
    options:
      members:
        - getToolkit
        - getToolkitTable
        - registerToolkit
        - getToolkitDocuments

::: hera.toolkit.abstractToolkit
    options:
      members:
        - addDataSource
        - getDataSourceDocument
        - getDataSourceData
        - getDataSourceList
        - getDataSourceTable
        - deleteDataSource

::: hera.datalayer.project.Project
    options:
      members:
        - addMeasurementsDocument
        - getMeasurementsDocuments
        - deleteMeasurementsDocuments
        - setConfig
        - getConfig
        - setCounter
        - getCounterAndAdd
