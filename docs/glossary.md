# Glossary

Key terms and concepts used throughout the Hera system.

---

### Project

A **Project** is the central workspace in Hera. It represents a named container that groups all data (measurements, simulations, cached results) and configurations together. Every interaction with data goes through a `Project` instance, which manages three MongoDB collections and a local files directory.

```python
from hera import Project
proj = Project(projectName="MY_PROJECT")
```

See [Core Concepts: Project](architecture/core_concepts.md#project) for technical details.

---

### Toolkit

A **Toolkit** is a domain-specific module that extends the `Project` class with specialized functionality for a particular type of data or analysis. Examples include `TopographyToolkit` for elevation data, `lowFreqToolKit` for meteorological station data, and `OFToolkit` for OpenFOAM simulations.

Each toolkit provides:

- **Data access** — Managing datasources of its domain
- **Analysis** — Processing and computation methods
- **Presentation** — Visualization and plotting capabilities

Once data is loaded into a project, you connect a suitable toolkit to perform operations on it.

See [Core Concepts: abstractToolkit](architecture/core_concepts.md#abstracttoolkit) for the base class.

---

### DataSource

A **DataSource** is a registered data entry within a toolkit. It represents external data (a file, URL, or Python object) along with all metadata needed for the toolkit to read and understand it. Each datasource has:

- **Name** — A human-readable identifier (e.g., `"YAVNEEL"`, `"SRTMGL1"`)
- **Resource** — Path to the actual data file
- **Data Format** — How to read the data (e.g., `"parquet"`, `"geopandas"`)
- **Version** — A `(major, minor, patch)` tuple for version management

```python
# Register a datasource
toolkit.addDataSource("YAVNEEL", "/data/YAVNEEL.parquet", "parquet", version=[0, 0, 1])

# Retrieve data
df = toolkit.getDataSourceData("YAVNEEL")
```

---

### Repository

A **Repository** is a JSON file that declares a collection of datasources, configurations, and documents organized by toolkit name. It serves as a blueprint for populating a project with data.

When a repository is added to a project, its contents are automatically loaded, registering all declared datasources, setting configurations, and creating the necessary MongoDB documents.

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

See [Data Layer: Repository JSON](architecture/data_layer.md#repository-json-structure) for the full format.

---

### Document

A **Document** is a MongoDB record that represents a piece of data in Hera. Every document has:

| Field | Description |
|-------|-------------|
| `projectName` | The project it belongs to |
| `_cls` | Type discriminator: `Metadata.Measurements`, `Metadata.Simulations`, or `Metadata.Cache` |
| `type` | Application-defined type tag (e.g., `"ToolkitDataSource"`) |
| `resource` | Path to the actual data file or inline value |
| `dataFormat` | How to interpret the resource (e.g., `"parquet"`, `"JSON_dict"`) |
| `desc` | Free-form metadata dictionary |

Documents are the fundamental unit of data organization in Hera.

---

### Collection

A **Collection** is a MongoDB collection that stores documents of a particular type. Hera uses three collections:

| Collection | Class | Purpose |
|------------|-------|---------|
| **Measurements** | `Measurements_Collection` | Observational data, toolkit datasources |
| **Simulations** | `Simulations_Collection` | Simulation model outputs |
| **Cache** | `Cache_Collection` | Intermediate results, configurations |

Each collection provides `addDocument()`, `getDocuments()`, and `deleteDocuments()` methods.

---

### Config

**Config** is a per-project key-value store for settings and parameters. It is stored as a special Cache document with `type = "<projectName>__config__"`. Toolkits use config to store defaults (e.g., default datasource name, default CRS).

```python
# Set configuration
proj.setConfig(defaultSRTM="SRTMGL1", defaultCRS=4326)

# Get configuration
config = proj.getConfig()
print(config["defaultSRTM"])  # "SRTMGL1"
```

---

### Counter

A **Counter** is an atomic integer stored within the project config, used for generating sequential IDs. Counters are thread-safe and support atomic read-and-increment operations.

```python
# Define a counter
proj.setCounter("experimentID", defaultValue=0)

# Atomically get and increment
current_id = proj.getCounterAndAdd("experimentID", addition=1)
```

---

### Version

A **Version** is a three-element tuple `[major, minor, patch]` used to manage multiple versions of the same datasource. The system supports:

- **Explicit versioning** — Request a specific version: `getDataSourceData("YAVNEEL", version=[0, 0, 2])`
- **Default version** — Set a default for a datasource: `setDataSourceDefaultVersion("YAVNEEL", [0, 0, 2])`
- **Latest version** — If no version is specified and no default is set, the highest version is returned

---

### ToolkitHome

**ToolkitHome** is the singleton registry that manages all available toolkits. It maintains a static dictionary of built-in toolkits and supports dynamic registration of custom toolkits at runtime. Access it via:

```python
from hera import toolkitHome

# Get a toolkit instance
tk = toolkitHome.getToolkit("MeteoLowFreq", projectName="MY_PROJECT")
```

See [Core Concepts: ToolkitHome](architecture/core_concepts.md#toolkithome) for technical details.

---

### dataToolkit

The **dataToolkit** is a special toolkit responsible for repository management. It operates on the `defaultProject` and handles:

- Registering repository JSON files
- Loading all datasources from a repository into a project
- Resolving relative paths in repository JSONs

See [Data Layer: Repository Pipeline](architecture/data_layer.md#repository-json-structure) for details.

---

### datatypes

The **datatypes** class defines all supported data format constants and the dispatch logic for reading/writing data. Each format constant (e.g., `PARQUET`, `NETCDF_XARRAY`, `GEOPANDAS`) maps to a specific reader/writer implementation.

See [Data Layer: datatypes](architecture/data_layer.md#the-datatypes-system) for the complete format table.

---

### Analysis Layer

The **Analysis Layer** is a property of each toolkit (`toolkit.analysis`) that provides domain-specific data processing methods. For example:

- `lowFreqToolKit.analysis.addDatesColumns()` — Add temporal columns to meteorological data
- `TopographyToolkit.analysis.calculateStatistics()` — Compute elevation statistics

---

### Presentation Layer

The **Presentation Layer** is a property of each toolkit (`toolkit.presentation`) that provides visualization and plotting methods. For example:

- `lowFreqToolKit.presentation.dailyPlots.plotScatter()` — Scatter plot of daily data
- `lowFreqToolKit.presentation.seasonalPlots.plotProbContourf_bySeason()` — Seasonal probability contours

---

### Expected Output

An **Expected Output** is a serialized file containing the known-correct result of a test. These files are organized into **result sets** (named directories) and are used by the comparison helpers to validate test results.

See [Testing Flow: Expected Output Management](testing/flow.md#expected-output-management) for details.

---

### Result Set

A **Result Set** is a named collection of expected output files, stored in a directory under `expected/`. The default result set is called `BASELINE`. Alternative result sets can be created for regression testing.
