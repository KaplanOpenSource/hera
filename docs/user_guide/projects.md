# Projects

A **Project** is Hera's unit of organization. Every document in the database — whether a measurement, simulation, or cached result — carries a `projectName` field that associates it with a project.

---

## Projects are defined by their documents

There is no master table of project names in Hera. A project **exists** as long as there are documents with its name. When the last document with `projectName="WindStudy"` is deleted, the project effectively ceases to exist. Creating a project is simply creating the first document with that name.

```python
from hera import Project

# This connects to (or implicitly creates) a project
proj = Project(projectName="WindStudy")

# List all projects that have at least one document
from hera.datalayer import getProjectList
getProjectList()
# ['WindStudy', 'CoastalSim', 'RiskAnalysis']
```

From the CLI:

```bash
# List all projects
hera-project project list

# Create a project directory with a caseConfiguration.json
hera-project project create WindStudy --directory /data/wind_study
```

---

## Directory-based projects

You don't always need to specify the project name explicitly. If you omit `projectName`, Hera looks for a `caseConfiguration.json` file in the current working directory:

```json
{
    "projectName": "WindStudy"
}
```

If the file is found, the project name is loaded automatically:

```python
# When run from a directory containing caseConfiguration.json:
proj = Project()  # projectName loaded from the file
```

This is the recommended workflow — create the project once with the CLI, then any script run from that directory connects to the right project:

```bash
# Create the project directory
hera-project project create WindStudy --directory /data/wind_study

# Work from that directory
cd /data/wind_study
python my_analysis.py   # Project() picks up "WindStudy" automatically
```

Toolkits follow the same convention — `toolkitHome.getToolkit("MeteoLowFreq")` without a `projectName` reads from `caseConfiguration.json` too.

If no `caseConfiguration.json` exists and no name is provided, Hera uses a read-only **default project** used internally for repository management.

---

## Three document collections

Each project organizes its data into three MongoDB collections:

| Collection | Role | Methods |
|-----------|------|---------|
| **Measurements** | Raw input data (station files, GIS data, sensor readings) | `addMeasurementsDocument`, `getMeasurementsDocuments`, `deleteMeasurementsDocuments` |
| **Simulations** | Computational model output (CFD results, dispersion runs) | `addSimulationsDocument`, `getSimulationsDocuments`, `deleteSimulationsDocuments` |
| **Cache** | Derived or intermediate results (statistics, function caches) | `addCacheDocument`, `getCacheDocuments`, `deleteCacheDocuments` |

All three collections share the same document structure and the same query interface. The separation is purely organizational — it helps you understand the provenance of each piece of data.

For detailed examples of adding, querying, and loading data, see [Working with Data](working_with_data.md).

---

## Project configuration

Each project has a **config** — a key-value store persisted in the database as a special cache document. Use it for project-level settings that should survive between sessions.

```python
proj = Project(projectName="WindStudy")

# Set configuration values
proj.setConfig(
    defaultStation="YAVNEEL",
    outputCRS=2039,
    domainSize={"width": 5000, "height": 5000}
)

# Read configuration
config = proj.getConfig()
print(config["defaultStation"])   # "YAVNEEL"
print(config["outputCRS"])        # 2039
print(config["domainSize"])       # {"width": 5000, "height": 5000}

# Update a single key (other keys are preserved)
proj.setConfig(outputCRS=4326)
```

**`initConfig`** sets values only if the keys don't already exist — useful for setting defaults without overwriting user choices:

```python
# These only take effect if the keys are not already set
proj.initConfig(
    defaultStation="BET_DAGAN",
    outputCRS=2039
)
```

Toolkits use the same config mechanism — each toolkit's settings are stored in the project config under toolkit-specific keys (e.g., `YAVNEEL_defaultVersion` for data source version defaults).

---

## Counters

Projects have built-in **atomic counters** — named integers stored in the database that increment safely even with concurrent access. A common use case is generating unique identifiers for output files.

```python
proj = Project(projectName="WindStudy")

# getCounterAndAdd returns the current value and increments atomically.
# On first call the counter is created starting at 0.
run_id = proj.getCounterAndAdd("simulation_run")  # 0
output = f"/data/results/run_{run_id}.nc"

run_id = proj.getCounterAndAdd("simulation_run")  # 1
output = f"/data/results/run_{run_id}.nc"
```

This is what `saveData` uses internally to generate unique file names.

| Method | Description |
|--------|-------------|
| `getCounterAndAdd(name, addition=1)` | Return current value and increment. Creates the counter at 0 if it doesn't exist. |
| `getCounter(name)` | Return current value without incrementing. Returns `None` if not defined. |
| `setCounter(name, defaultValue=0)` | Create or reset a counter. |

Counters are per-project and stored inside the project config document.

---

## Data manipulation methods

Projects provide methods for adding, querying, deleting, and saving data. Here is a summary — for full examples with code, see [Working with Data](working_with_data.md).

### Adding documents

| Method | Description |
|--------|-------------|
| `addMeasurementsDocument(resource, dataFormat, type, desc)` | Add a measurement document |
| `addSimulationsDocument(resource, dataFormat, type, desc)` | Add a simulation document |
| `addCacheDocument(resource, dataFormat, type, desc)` | Add a cache document |
| `saveData(name, data, desc, kind)` | Auto-detect format, save file to disk, create document |
| `saveMeasurementData(name, data, desc)` | Save as measurement (shorthand) |
| `saveSimulationData(name, data, desc)` | Save as simulation (shorthand) |
| `saveCacheData(name, data, desc)` | Save as cache (shorthand) |

See [Adding data](working_with_data.md#adding-data) for examples.

### Querying documents

| Method | Description |
|--------|-------------|
| `getMeasurementsDocuments(**filters)` | Query measurement documents |
| `getSimulationsDocuments(**filters)` | Query simulation documents |
| `getCacheDocuments(**filters)` | Query cache documents |
| `getAllDocuments(**filters)` | Query across all collections |
| `getDocumentByID(id)` | Get a single document by its MongoDB ID |
| `getMetadata()` | Return all document descriptions as a DataFrame |

See [Querying the database](working_with_data.md#querying) for examples of basic, nested, and structured queries.

### Deleting documents

| Method | Description |
|--------|-------------|
| `deleteMeasurementsDocuments(**filters)` | Delete matching measurement documents |
| `deleteSimulationsDocuments(**filters)` | Delete matching simulation documents |
| `deleteCacheDocuments(**filters)` | Delete matching cache documents |

### Export and import

| Method | Description |
|--------|-------------|
| `export(path)` | Export all project documents to a zip file |
| `Project.load(proj, path, is_hard_import)` | Import documents from an exported zip |

```bash
# CLI equivalents
hera-project project dump WindStudy --format json --fileName backup.json
hera-project project load WindStudy backup.json
```
