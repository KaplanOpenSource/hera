# Quick Reference / Cheat Sheet

A one-page quick reference for common Hera operations.

---

## CLI Commands

### Project Management

| Command | Description |
|--------|-------------|
| `hera-project project list` | List all projects |
| `hera-project project create MY_PROJECT` | Create a new project |
| `hera-project project dump MY_PROJECT` | Export project data |
| `hera-project project load MY_PROJECT file.json` | Import project data |
| `hera-project project measurements list --project MY_PROJECT` | List measurements |

### Repository Management

| Command | Description |
|--------|-------------|
| `hera-project repository list` | List registered repositories |
| `hera-project repository add myRepo /path/to/repo.json` | Register a repository |
| `hera-project repository load myRepo MY_PROJECT` | Load repository into project |
| `hera-project repository show myRepo` | Show repository contents |

### Toolkit Management

| Command | Description |
|--------|-------------|
| `hera-toolkit list --project MY_PROJECT` | List available toolkits |
| `hera-toolkit load --project MY_PROJECT --name MeteoLowFreq` | Load a toolkit |
| `hera-toolkit register --project MY_PROJECT --cls pkg.Class --name MyToolkit` | Register custom toolkit |

### Version Management

| Command | Description |
|--------|-------------|
| `hera-project project version display MY_PROJECT` | Show all datasource versions |
| `hera-project project version update MY_PROJECT YAVNEEL "0,0,2"` | Set default version |

---

## Python API Patterns

### Project Initialization

```python
from hera import Project

# Basic initialization
proj = Project(projectName="MY_PROJECT")

# With custom connection
proj = Project(projectName="MY_PROJECT", connectionName="production")

# With files directory
proj = Project(projectName="MY_PROJECT", filesDirectory="/data/myproject")
```

### Toolkit Access

```python
from hera import toolkitHome

# Get a toolkit instance
toolkit = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_LOWFREQ, projectName="MY_PROJECT")

# List all toolkits
table = toolkitHome.getToolkitTable("MY_PROJECT")
print(table)
```

### Datasource Operations

```python
# Get datasource data
data = toolkit.getDataSourceData("YAVNEEL")

# Get datasource document
doc = toolkit.getDataSourceDocument("YAVNEEL", version=[0, 0, 1])

# List all datasources
sources = toolkit.getDataSourceList()

# Add a datasource
toolkit.addDataSource(
    dataSourceName="NEW_DATA",
    resource="/path/to/data.parquet",
    dataFormat="parquet",
    version=[0, 0, 1]
)

# Set default version
toolkit.setDataSourceDefaultVersion("YAVNEEL", [0, 0, 2])
```

### Document Operations

```python
# Add measurement document
proj.addMeasurementsDocument(
    resource="/path/to/file.parquet",
    dataFormat="parquet",
    type="MyType",
    desc={"key": "value"}
)

# Query documents
docs = proj.getMeasurementsDocuments(type="ToolkitDataSource", toolkit="MeteoLowFreq")

# Delete documents
proj.deleteMeasurementsDocuments(type="MyType")
```

### Configuration

```python
# Set config
proj.setConfig(defaultSRTM="SRTMGL1", defaultCRS=4326)

# Get config
config = proj.getConfig()
print(config["defaultSRTM"])

# Counters
proj.setCounter("experimentID", defaultValue=0)
id = proj.getCounterAndAdd("experimentID", addition=1)
```

---

## Repository JSON Structure

```json
{
    "<ToolkitName>": {
        "Config": {
            "key": "value"
        },
        "DataSource": {
            "<name>": {
                "isRelativePath": "True",
                "item": {
                    "resource": "path/to/data.parquet",
                    "dataFormat": "parquet",
                    "version": [0, 0, 1],
                    "desc": {}
                }
            }
        },
        "Measurements": {
            "<name>": {
                "isRelativePath": "True",
                "item": {
                    "resource": "path/to/file.shp",
                    "dataFormat": "geopandas",
                    "type": "MyType",
                    "desc": {}
                }
            }
        }
    }
}
```

---

## Data Format Constants

| Constant | Value | Use Case |
|----------|-------|----------|
| `STRING` | `"string"` | Plain text paths |
| `PARQUET` | `"parquet"` | Tabular data (pandas/dask) |
| `NETCDF_XARRAY` | `"netcdf_xarray"` | NetCDF files (xarray) |
| `GEOPANDAS` | `"geopandas"` | Shapefiles, GeoPackages |
| `JSON_DICT` | `"JSON_dict"` | JSON configuration files |
| `CSV_PANDAS` | `"csv_pandas"` | CSV files (pandas) |
| `PICKLE` | `"pickle"` | Python pickle files |
| `CLASS` | `"Class"` | Dynamic Python classes |

```python
from hera.datalayer import datatypes

# Usage
toolkit.addDataSource(
    dataSourceName="DATA",
    resource="/path/to/file.parquet",
    dataFormat=datatypes.PARQUET,
    version=[0, 0, 1]
)
```

---

## Environment Variables

| Variable | Default | Purpose |
|----------|---------|---------|
| `TEST_HERA` | `~/hera_unittest_data` | Test data root |
| `RESULT_SET` | `BASELINE` | Expected output result set |
| `PREPARE_EXPECTED_OUTPUT` | unset | Generate test baselines |
| `MPLBACKEND` | system | Matplotlib backend |
| `GDF_TOL_AREA` | `1e-7` | Geometry comparison tolerance |

---

## Common Patterns

### Load Repository into Project

```python
from hera.utils.data.toolkit import dataToolkit

dt = dataToolkit()
dt.loadAllDatasourcesInRepositoryToProject(
    projectName="MY_PROJECT",
    repositoryName="myRepo",
    overwrite=True
)
```

### Load Repository JSON Directly

```python
import json
from hera.utils.data.toolkit import dataToolkit

with open("repository.json") as f:
    repo_json = json.load(f)

dt = dataToolkit()
dt.loadAllDatasourcesInRepositoryJSONToProject(
    projectName="MY_PROJECT",
    repositoryJSON=repo_json,
    basedir="/path/to/repo/dir",
    overwrite=True
)
```

### Toolkit Analysis Workflow

```python
# Get toolkit
lf = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_LOWFREQ, projectName="MY_PROJECT")

# Load data
df = lf.getDataSourceData("YAVNEEL").compute()

# Analysis
enriched = lf.analysis.addDatesColumns(df, datecolumn="datetime")
hourly = lf.analysis.calcHourlyDist(enriched, field="T", density=True)

# Presentation
ax = lf.presentation.dailyPlots.plotScatter(enriched, plotField="RH")
```

### Export/Import Project

```python
# Export
proj.export("/path/to/backup.zip")

# Import
Project.load(proj, "/path/to/backup.zip", is_hard_import=False)
```

---

## MongoDB Config File

Location: `~/.pyhera/config.json`

```json
{
    "connectionName": {
        "username": "user",
        "password": "pass",
        "dbIP": "localhost:27017",
        "dbName": "hera_db"
    }
}
```

---

## Quick Links

- [Installation Guide](../getting_started/installation.md)
- [CLI Reference](../cli/reference.md)
- [Repository Schema](repository_schema.md)
- [Best Practices](../guides/best_practices.md)
- [Troubleshooting](../troubleshooting.md)
