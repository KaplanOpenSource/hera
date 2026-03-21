# Repositories

A **Repository** is a JSON file that describes a collection of data sources and documents for one or more toolkits. Instead of adding data sources to a project one by one, you define them all in a repository file and load them in a single step.

---

## Why repositories?

Without repositories, setting up a project requires manually adding each data source:

```python
topo.addDataSource("Israel_DEM", resource="/data/srtm/dem.hgt", dataFormat="HDF")
lc.addDataSource("LandCover_2021", resource="/data/modis/lc.tif", dataFormat="geotiff")
meteo.addDataSource("YAVNEEL", resource="/data/meteo/yavneel.parquet", dataFormat="parquet")
# ... repeat for every data source
```

If you create a second project, you repeat the whole process. Repositories solve this — define once, load anywhere:

```bash
# Register once
hera-project repository add /path/to/my_repository.json

# Load into any project
hera-project project create WindStudy
# All data sources from all registered repositories are loaded automatically
```

---

## The user-level repository registry

Repositories are managed by a special toolkit called `dataToolkit`. It runs on the **default project** — a shared, user-level space that exists independently of your work projects.

When you register a repository with `hera-project repository add`, the path to the JSON file is stored as a data source in the default project. This means:

- Repositories are **per-user** — each user has their own registry of repositories
- Repositories persist across sessions — they survive until you remove them
- When you create a new project with the CLI, **all registered repositories are automatically loaded** into it

```bash
# Register a repository (stored in the default project)
hera-project repository add /data/repos/gis_data.json

# List all registered repositories
hera-project repository list

# Show contents of a repository
hera-project repository show gis_data

# Remove a repository from the registry
hera-project repository remove gis_data
```

### Loading repositories into existing projects

If you created a project before registering a repository, you can load it manually:

```bash
# Load a specific repository into a project
hera-project repository load gis_data WindStudy

# Load ALL registered repositories into a project
hera-project project updateRepositories --projectName WindStudy --overwrite
```

---

## Repository JSON format

A repository JSON maps **toolkit names** to sections of data. Each toolkit can have:

- **Config** — key-value settings applied to the toolkit's project config
- **DataSource** — named, versioned data sources (the most common section)
- **Measurements** / **Simulations** / **Cache** — raw documents added to the project
- **Function** — calls to named functions with parameters

### Basic example

```json
{
    "GIS_Raster_Topography": {
        "Config": {
            "defaultSRTM": "SRTMGL1"
        },
        "DataSource": {
            "SRTMGL1": {
                "isRelativePath": "True",
                "item": {
                    "resource": "data/srtm/srtmgl1.hgt",
                    "dataFormat": "HDF"
                }
            }
        }
    },
    "MeteoLowFreq": {
        "Config": {
            "defaultStation": "YAVNEEL"
        },
        "DataSource": {
            "YAVNEEL": {
                "isRelativePath": "True",
                "item": {
                    "resource": "data/meteo/yavneel.parquet",
                    "dataFormat": "parquet",
                    "desc": {
                        "station": "YAVNEEL",
                        "location": "Galilee"
                    }
                }
            }
        }
    }
}
```

### Section types

| Section | What it does | Example |
|---------|-------------|---------|
| **Config** | Sets toolkit configuration for the project | `{"defaultStation": "YAVNEEL"}` |
| **DataSource** | Registers named, versioned data sources | File paths with format and metadata |
| **Measurements** | Adds raw measurement documents | Sensor data, GIS files |
| **Simulations** | Adds simulation result documents | CFD output, model results |
| **Cache** | Adds cached/derived data documents | Pre-computed statistics |
| **Function** | Calls a named function with parameters | Custom loading logic |

### The `isRelativePath` flag

Each data item has an `isRelativePath` flag:

- `"True"` — the `resource` path is relative to the **repository JSON file's directory**. Hera resolves it to an absolute path when loading.
- `"False"` — the `resource` is an absolute path or URL, used as-is.

This makes repositories portable — you can move the repository JSON and its data together, and relative paths still work.

---

## What happens when a repository is loaded

When you load a repository into a project, Hera does the following for each toolkit section:

1. **Resolves the toolkit** — finds the toolkit class via `ToolkitHome.getToolkit()`
2. **Processes each section** in order:
   - **Config** → calls `toolkit.setConfig(**values)`
   - **DataSource** → calls `toolkit.addDataSource(name, resource, dataFormat)` for each item
   - **Measurements/Simulations/Cache** → calls the corresponding `addDocument` method
   - **Function** → calls the named method with the provided parameters
3. **Resolves paths** — converts relative paths to absolute using the repository file's directory

---

## How data sources are stored internally

A data source is **not a special object** — it's a regular measurement document with specific conventions. When a toolkit calls `addDataSource("YAVNEEL", ...)`, it creates a measurement document like this:

```json
{
    "_cls": "Metadata.Measurements",
    "projectName": "WindStudy",
    "type": "ToolkitDataSource",
    "resource": "/data/meteo/yavneel.parquet",
    "dataFormat": "parquet",
    "desc": {
        "toolkit": "lowFreqToolKit",
        "datasourceName": "YAVNEEL",
        "version": [1, 0, 0],
        "station": "YAVNEEL",
        "location": "Galilee"
    }
}
```

Key fields set by the toolkit:

| Field | Value | Purpose |
|-------|-------|---------|
| `type` | `"ToolkitDataSource"` | Marks this as a data source (not a regular document) |
| `desc.toolkit` | The toolkit name | Associates the data source with its toolkit |
| `desc.datasourceName` | The data source name | Used for lookup by name |
| `desc.version` | Version tuple as list | Versioning support |

This means you can also query data sources directly through the Project API if needed:

```python
# These are equivalent:
data = meteo.getDataSourceData("YAVNEEL")

# Direct query (what the toolkit does internally):
docs = proj.getMeasurementsDocuments(
    type="ToolkitDataSource",
    toolkit="lowFreqToolKit",
    datasourceName="YAVNEEL"
)
data = docs[0].getData()
```

The toolkit's `getDataSourceData` is a convenience wrapper that handles version resolution, default version lookup, and the query construction for you.

---

## Multiple repositories

You can register multiple repositories. They are all loaded when a project is created:

```bash
hera-project repository add /data/repos/gis_data.json
hera-project repository add /data/repos/meteo_data.json
hera-project repository add /data/repos/risk_agents.json

# All three are loaded into new projects
hera-project project create WindStudy
```

If two repositories define a data source with the same name for the same toolkit, the last one loaded wins (unless `--overwrite` is not set, in which case the existing one is kept).

---

## Next steps

- **[Working with Data Sources](working_with_data.md#data-sources)** — API details for versions, defaults, querying
- **[CLI Reference > Repository Management](../cli/reference.md#repository-management)** — full CLI command reference
- **[Projects > Project Lifecycle](projects.md#project-lifecycle)** — how repositories fit into the project setup flow
