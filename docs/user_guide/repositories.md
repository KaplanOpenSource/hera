# Repositories

A **Repository** is a JSON file that describes a collection of data sources and documents for one or more toolkits. Instead of adding data sources to a project one by one, you define them all in a repository file and load them in a single step.

---

## Why repositories?

Without repositories, setting up a project requires manually adding each data source:

```python
topo.addDataSource("SRTMGL1", resource="/data/srtm", dataFormat="string")   # directory of .hgt tiles
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

When you register a repository with `hera-project repository add`, the absolute path to the JSON file is stored as a data source (format `JSON_dict`) in the default project, under the file's base name (`gis_data.json` → `gis_data`). This means:

- Repositories are **per database connection** — by default the connection named after your OS user, so each user has their own registry of repositories
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
# Load a specific repository file into a project (takes the JSON file path, not the registered name)
hera-project repository load /data/repos/gis_data.json WindStudy

# Load ALL registered repositories into a project
hera-project project updateRepositories --projectName WindStudy --overwrite
```

---

## Repository JSON format

A repository JSON maps **toolkit names** to sections of data. Each toolkit can have:

- **Config** — key-value settings applied to the toolkit's project config
- **DataSource** — named, versioned data sources (the most common section)
- **Measurements** / **Simulations** / **Cache** — raw documents added to the project
- **Function** — calls to named `dataToolkit` methods with parameters

Section names are matched case-insensitively after title-casing (`DataSource` and `Datasource` both work). An unknown section name raises `ValueError`; a toolkit name that `ToolkitHome.getToolkit()` cannot resolve is skipped.

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
                    "resource": "data/srtm",
                    "dataFormat": "string"
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
| **Function** | Calls a named `dataToolkit` method; the value is a dict of keyword arguments or a list of such dicts, and `overwrite` is added to each call | Custom loading logic |

A **DataSource** item is passed to `toolkit.addDataSource(**item)`: `resource`, `dataFormat`, optional `version` (e.g. `[1, 0, 0]`, default `[0, 0, 1]`), and any other keys are stored in the document's `desc` (a `"desc"` key is stored as a nested `desc.desc` dictionary). Instead of `resource` you can give `resourceFilePath`: the path of a JSON file whose contents become the resource.

A **Measurements** / **Simulations** / **Cache** item is passed to `add<Type>Document(**item)`: `resource`, `dataFormat`, `type`, `desc`.

### The `isRelativePath` flag

Each **DataSource** entry has an `isRelativePath` flag next to `item`. It is required (`"True"`/`"False"` or a boolean; anything else fails the section):

- `"True"` — the `resource` path is relative to the **repository JSON file's directory**. Hera resolves it to an absolute path when loading.
- `"False"` — the `resource` is an absolute path, URL or inline value (e.g. a dict), used as-is.

For **Measurements** / **Simulations** / **Cache** entries the flag is not used: the `resource` is always joined to the repository directory, so a relative path is resolved and an absolute path stays absolute. Do not put `isRelativePath` inside `item` there; it would be passed to `add<Type>Document` and fail.

This makes repositories portable — you can move the repository JSON and its data together, and relative paths still work.

---

## What happens when a repository is loaded

When you load a repository into a project, Hera does the following for each toolkit section:

1. **Resolves the toolkit** — finds the toolkit via `ToolkitHome.getToolkit()` (static or dynamic); unknown toolkits are skipped
2. **Processes each section** in order, resolving relative paths against the repository file's directory:
   - **Config** → calls `toolkit.setConfig(**values)`
   - **DataSource** → calls `toolkit.addDataSource(dataSourceName=name, **item)` for each item. An existing data source with the same name is kept unless `overwrite` is set
   - **Measurements/Simulations/Cache** → queries for a matching document and calls the corresponding `add<Type>Document` if there is none (with `overwrite`, updates the existing document)
   - **Function** → calls the named `dataToolkit` method with the provided parameters

An error inside a section is logged and that section is skipped; loading continues with the next section.

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
        "toolkit": "lowFreqMeteorology",
        "datasourceName": "YAVNEEL",
        "version": [1, 0, 0],
        "desc": {
            "station": "YAVNEEL",
            "location": "Galilee"
        }
    }
}
```

Key fields set by the toolkit:

| Field | Value | Purpose |
|-------|-------|---------|
| `type` | `"ToolkitDataSource"` | Marks this as a data source (not a regular document) |
| `desc.toolkit` | The toolkit's internal `toolkitName` (e.g. `lowFreqMeteorology` for `MeteoLowFreq`) | Associates the data source with its toolkit |
| `desc.datasourceName` | The data source name | Used for lookup by name |
| `desc.version` | Version tuple as list | Versioning support |

This means you can also query data sources directly through the Project API if needed:

```python
# These are equivalent:
data = meteo.getDataSourceData("YAVNEEL")

# Direct query (what the toolkit does internally):
docs = proj.getMeasurementsDocuments(
    type="ToolkitDataSource",
    toolkit=meteo.toolkitName,     # "lowFreqMeteorology"
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

If two repositories define a data source with the same name for the same toolkit, the first one loaded is kept and the later one is skipped with an error in the log. With `--overwrite`, the last one loaded wins.

---

## Populating projects

New projects are automatically populated with repositories when created via `hera-project project create`. But if you:

- Added a new repository after creating projects
- Want to ensure all projects have the latest data sources
- Created projects before registering repositories

Use the **populate** command to bulk-load all repositories into all projects at once:

```bash
# Populate ALL existing projects with all registered repositories
make populate

# Populate a specific project
make populate-project PROJECT=WindStudy
```

Or via the CLI directly:

```bash
# All projects
hera-project project populate --overwrite

# Single project
hera-project project populate --projectName WindStudy --overwrite

# Without --overwrite: skips data sources that already exist
hera-project project populate
```

The `--overwrite` flag controls whether existing data sources are replaced. Without it, only new data sources are added.

### Typical workflow

```bash
# 1. Register your repositories (one time)
hera-project repository add /data/repos/gis_data.json
hera-project repository add /data/repos/meteo_data.json

# 2. Create some projects
hera-project project create WindStudy --directory /data/wind_study
hera-project project create CoastalSim --directory /data/coastal

# 3. Later, add a new repository
hera-project repository add /data/repos/risk_agents.json

# 4. Populate all projects with the new repository
make populate
# Both WindStudy and CoastalSim now have risk agent data sources
```

---

## Exporting documents to a repository

The reverse direction is also supported: `repository export` writes documents of a project into a repository JSON file under a given toolkit key, and registers the file (unless `--no-register`):

```bash
# Export ALL documents of the project under the "MeteoLowFreq" key
hera-project repository export meteo_data_v1 --toolkitName MeteoLowFreq --projectName WindStudy

# Export only selected documents (repeat --documentId), merging into an existing file
hera-project repository export /data/repos/meteo_data.json --toolkitName MeteoLowFreq \
    --projectName WindStudy --documentId 64f... --documentId 650...
```

`repositoryName` is a registered repository name or a path to a `.json` file (an unregistered bare name creates `<name>.json` in the current directory). Without `--projectName` the documents come from the default project. Options: `--idStrategy contentHash|objectId` (how duplicates are identified, default `contentHash`), `--mode add|override` (merge, or merge and then deduplicate the whole file; default `add`) and `--overwrite` (replace existing entries on identity match). The Python equivalent is `dataToolkit().exportDocumentsToRepository(toolkitName=..., repositoryName=..., projectName=..., documents=...)`.

---

## Next steps

- **[Working with Data Sources](working_with_data.md#data-sources)** — API details for versions, defaults, querying
- **[CLI Reference > Repository Management](../cli/reference.md#repository-management)** — full CLI command reference
- **[Projects > Project Lifecycle](projects.md#project-lifecycle)** — how repositories fit into the project setup flow
