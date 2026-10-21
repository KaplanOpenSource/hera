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
| `hera-project repository add /path/to/myRepo.json` | Register a repository (named after the file: `myRepo`) |
| `hera-project repository export myRepo --toolkitName MeteoLowFreq` | Export project documents into a repository JSON |
| `hera-project repository load /path/to/myRepo.json MY_PROJECT` | Load a repository JSON file into a project (takes the file path, not the registered name) |
| `hera-project repository show myRepo` | Show repository contents |

### Toolkit Management

| Command | Description |
|--------|-------------|
| `hera-project addToolkit MyToolkit /path/to/toolkit` | Register a custom toolkit |
| `hera-project addToolkit MyToolkit /path/to/toolkit --version "0,0,2" --overwrite` | Re-register with a new version |

There is no `hera-toolkit` command; list and load toolkits from Python (`toolkitHome.getToolkitTable(...)`, `toolkitHome.getToolkit(...)`).

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
                    "version": [0, 0, 1]
                }
            }
        },
        "Measurements": {
            "<name>": {
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
| `GEOTIFF` | `"geotiff"` | GeoTIFF rasters |
| `ZARR_XARRAY` | `"zarr_xarray"` | Zarr stores (xarray) |
| `HDF` | `"HDF"` | HDF files |
| `JSON_PANDAS` / `JSON_GEOPANDAS` | `"JSON_pandas"` / `"JSON_geopandas"` | JSON-serialized (Geo)DataFrames |
| `NUMPY_ARRAY` / `NUMPY_DICT_ARRAY` | `"numpy_array"` / `"numpy_dict_array"` | NumPy arrays / dict of arrays |
| `TIME`, `DICT`, `IMAGE` | `"time"`, `"dict"`, `"image"` | Timestamps, in-document dicts, images |

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
| `MPLBACKEND` | system | Matplotlib backend |
| `GDF_TOL_AREA` | `1e-7` | Geometry comparison tolerance |
| `HERA_REPO_ROOT` | unset | Base dir for relative dynamic-toolkit / experiment paths |
| `RAG_*` | see env vars page | RAG search settings (e.g. `RAG_OLLAMA_BASE_URL`) |

Full list: [Environment Variables](../configuration/env_vars.md).

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
x_mid, y_mid, hist = lf.analysis.calcHourlyDist(enriched, Field="T", normalization="density")

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

- [Installation Guide](../user_guide/installation.md)
- [CLI Reference](../cli/reference.md)
- [Repository Schema](repository_schema.md)
- [Best Practices](../guides/best_practices.md)
- [Troubleshooting](../troubleshooting.md)
