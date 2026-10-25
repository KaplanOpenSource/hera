# Working with Data

This page covers the practical details of storing, querying, and loading data in Hera. For the high-level overview, see [Key Concepts](concepts.md).

> **Hands-on tutorials:** See the [DataSource](../tutorials/index.md) and [Repository](../tutorials/index.md) notebook tutorials for interactive walkthroughs with real output.

---

## Adding data

### Manual document creation

Use `addMeasurementsDocument`, `addSimulationsDocument`, or `addCacheDocument` to register a file with metadata:

```python
from hera import Project

proj = Project(projectName="WindStudy")

proj.addMeasurementsDocument(
    resource="/data/station_A.parquet",
    dataFormat="parquet",
    type="WeatherStation",
    desc={
        "station": "A",
        "location": "Haifa",
        "elevation": 120,
        "variables": ["temperature", "wind_speed"]
    }
)
```

### Using `Project.datatypes` for format names

You don't need to remember the exact format strings. Every `Project` instance (and every toolkit) exposes a `datatypes` object with all format constants:

```python
proj = Project(projectName="WindStudy")

# Use the datatypes constants instead of raw strings
proj.addMeasurementsDocument(
    resource="/data/station_A.parquet",
    dataFormat=proj.datatypes.PARQUET,       # instead of "parquet"
    type="WeatherStation",
    desc={"station": "A"}
)

proj.addSimulationsDocument(
    resource="/data/result.nc",
    dataFormat=proj.datatypes.NETCDF_XARRAY, # instead of "netcdf_xarray"
    type="DispersionRun",
    desc={"scenario": "baseline"}
)

proj.addCacheDocument(
    resource="/data/stats.json",
    dataFormat=proj.datatypes.JSON_DICT,     # instead of "JSON_dict"
    type="DailyStats",
    desc={"period": "2024"}
)
```

You can also import `datatypes` directly:

```python
from hera.datalayer import datatypes

# Same constants available standalone
datatypes.PARQUET          # "parquet"
datatypes.NETCDF_XARRAY    # "netcdf_xarray"
datatypes.GEOPANDAS        # "geopandas"
```

### Saving data with auto-detection

For common Python objects, `saveData` / `saveCacheData` / `saveMeasurementData` / `saveSimulationData` auto-detect the format, save the file under `<filesDirectory>/cache/<name>_<counter>.<ext>`, and create the document in one call (`type` defaults to `name`; pass `dataFormat=` to override detection):

```python
import pandas as pd

df = pd.DataFrame({"temp": [20, 21, 22], "wind": [5, 6, 7]})

# Hera detects DataFrame -> parquet, saves the file, creates the document
proj.saveCacheData(
    name="daily_summary",
    data=df,
    desc={"period": "2024-01", "station": "A"}
)
```

The auto-detection mapping:

| Python type | Data format | File extension |
|------------|-------------|---------------|
| `str` | `string` | `.txt` |
| `pandas.DataFrame` | `parquet` | `.parquet` |
| `pandas.Series` | `JSON_pandas` | `.json` |
| `dask.DataFrame` | `parquet` | `.parquet` |
| `geopandas.GeoDataFrame` | `geopandas` | `.gpkg` |
| `xarray.DataArray` | `netcdf_xarray` | `.zarr` (written as NetCDF despite the extension) |
| `numpy.ndarray` | `numpy_array` | `.npy` |
| `dict`, `list`, `bytes`, any other object | `pickle` | `.pckle` |

`xarray.Dataset` is not in the mapping and falls back to `pickle`; pass `dataFormat=proj.datatypes.NETCDF_XARRAY` to store it as NetCDF.

### Using counters for unique file names

Projects have built-in **atomic counters** — named integers stored in the database that increment safely even when multiple processes run in parallel. A common use case is generating unique file names for output data.

```python
proj = Project(projectName="WindStudy")

# getCounterAndAdd increments atomically and returns the new value.
# The first call only creates the counter at 0 and returns 0.
run_id = proj.getCounterAndAdd("simulation_run")  # 0
output_path = f"/data/results/dispersion_run_{run_id}.nc"

# Next call returns 1, then 2, etc.
run_id = proj.getCounterAndAdd("simulation_run")  # 1
output_path = f"/data/results/dispersion_run_{run_id}.nc"

# Save with the generated file name
proj.addSimulationsDocument(
    resource=output_path,
    dataFormat=proj.datatypes.NETCDF_XARRAY,
    type="DispersionRun",
    desc={"run_id": run_id, "scenario": "baseline"}
)
```

This is exactly what `saveData` does internally — it uses a counter named after the data to generate unique file paths automatically.

**Counter methods:**

| Method | Description |
|--------|-------------|
| `getCounterAndAdd(name, addition=1)` | Add `addition` atomically and return the new value. If the counter doesn't exist, create it at 0 and return 0. |
| `getCounter(name)` | Return the current value without incrementing. Returns `None` if the counter doesn't exist. |
| `setCounter(name, defaultValue=0)` | Create or reset a counter to a specific value. |

Counters are stored per-project in the project's config document, so each project has its own independent counters.

---

## The `type` field

The `type` field is a string label that you define to categorize documents within a collection. It has no fixed vocabulary — you choose type names that make sense for your domain. Hera uses `type` as the primary grouping mechanism for documents.

For example, within the Measurements collection you might have:

| `type` value | What it represents |
|-------------|-------------------|
| `"WeatherStation"` | Meteorological station data files |
| `"ToolkitDataSource"` | Versioned data sources managed by toolkits |
| `"Experiment_rawData"` | Raw experiment data files |
| `"BuildingFootprints"` | GIS building vector data |
| `"ElevationGrid"` | DEM / topography raster data |

Toolkits use `type` internally to organize their data — for instance, all toolkit data sources are stored with `type="ToolkitDataSource"`. When you query documents, `type` is typically the first filter you apply.

---

## Data formats (resource types) {#resource-types}

The `dataFormat` field tells Hera how to read the `resource`. When you call `doc.getData()`, Hera dispatches to the correct handler based on this field. The supported formats are:

| Format | `dataFormat` value | `datatypes` constant | Python type returned | File extension |
|--------|-------------------|---------------------|----------------------|---------------|
| **Parquet** | `"parquet"` | `PARQUET` | `dask.DataFrame` by default, `pandas.DataFrame` with `getData(usePandas=True)` (set automatically when the file was saved from pandas) | `.parquet` |
| **CSV** | `"csv_pandas"` | `CSV_PANDAS` | `pandas.DataFrame` | `.csv` |
| **HDF5** | `"HDF"` | `HDF` | `dask.DataFrame` (`pandas` with `usePandas=True`); resource is `{"path": ..., "key": ...}`, read-only | `.hdf` |
| **NetCDF** | `"netcdf_xarray"` | `NETCDF_XARRAY` | `xarray.Dataset` | `.nc` |
| **Zarr** | `"zarr_xarray"` | `ZARR_XARRAY` | `xarray.Dataset` (loading is currently broken in master: `getData` raises `AttributeError`) | `.zarr` |
| **GeoPackage** | `"geopandas"` | `GEOPANDAS` | `geopandas.GeoDataFrame` | `.gpkg` |
| **GeoJSON** | `"JSON_geopandas"` | `JSON_GEOPANDAS` | `geopandas.GeoDataFrame` | `.json` |
| **GeoTIFF** | `"geotiff"` | `GEOTIFF` | GDAL dataset (read-only; requires GDAL) | `.tif` |
| **JSON (dict)** | `"JSON_dict"` | `JSON_DICT` | `dict` | `.json` |
| **JSON (pandas)** | `"JSON_pandas"` | `JSON_PANDAS` | `pandas.DataFrame` | `.json` |
| **NumPy array** | `"numpy_array"` | `NUMPY_ARRAY` | `numpy.ndarray` | `.npy` |
| **NumPy dict** | `"numpy_dict_array"` | `NUMPY_DICT_ARRAY` | dict of `numpy.ndarray` | `.npz` |
| **Image** | `"image"` | `IMAGE` | `numpy.ndarray` (pixel data) | `.png` |
| **Pickle** | `"pickle"` | `PICKLE` | any Python object | `.pckle` |
| **String** | `"string"` | `STRING` | `str` | `.txt` |
| **Timestamp** | `"time"` | `TIME` | `pandas.Timestamp` | — |
| **Dict** | `"dict"` | `DICT` | `dict` (stored inline in resource) | — |
| **Class** | `"Class"` | `CLASS` | class instance or class object | — |

---

## Querying the database {#querying}

Hera uses [MongoEngine](http://mongoengine.org/) under the hood. When you call `getMeasurementsDocuments()`, `resource`, `dataFormat` and `type` match the top-level fields and every other keyword argument is a `desc` field. The `desc` arguments are passed through `dictToMongoQuery`, so they are flattened using double-underscore (`__`) notation — the same convention MongoEngine uses for nested document queries.

### Basic queries

Filter by top-level fields:

```python
# Find by type
docs = proj.getMeasurementsDocuments(type="WeatherStation")

# Find by type and data format
docs = proj.getMeasurementsDocuments(type="WeatherStation", dataFormat="parquet")
```

### Querying nested `desc` fields

Pass desc fields directly as keyword arguments:

```python
# Find stations in Haifa
docs = proj.getMeasurementsDocuments(type="WeatherStation", location="Haifa")

# Find stations at exactly 120 m elevation
docs = proj.getMeasurementsDocuments(type="WeatherStation", elevation=120)

# MongoEngine operators work too: stations above 100 m
docs = proj.getMeasurementsDocuments(type="WeatherStation", elevation__gt=100)
```

Behind the scenes, `location="Haifa"` becomes a MongoDB query on `desc.location`, using MongoEngine's `__` syntax: `desc__location="Haifa"`. A list value is expanded element by element (`variables=["a", "b"]` becomes `desc__variables__0="a", desc__variables__1="b"`).

### Structured (nested) metadata queries

When your `desc` has nested dictionaries, the double-underscore notation traverses the hierarchy:

```python
# Store a document with nested metadata
proj.addSimulationsDocument(
    resource="/data/sim_result.nc",
    dataFormat=proj.datatypes.NETCDF_XARRAY,
    type="DispersionRun",
    desc={
        "scenario": {
            "source": "factory_A",
            "release_rate": 10.0,
            "wind": {
                "speed": 5.0,
                "direction": 270
            }
        },
        "status": "completed"
    }
)

# Query by nested fields — Hera flattens them with __ notation
docs = proj.getSimulationsDocuments(
    type="DispersionRun",
    scenario__source="factory_A"            # desc.scenario.source
)

docs = proj.getSimulationsDocuments(
    type="DispersionRun",
    scenario__wind__speed=5.0               # desc.scenario.wind.speed
)

# Combine multiple nested filters
docs = proj.getSimulationsDocuments(
    type="DispersionRun",
    scenario__source="factory_A",
    scenario__wind__direction=270,
    status="completed"
)

# Equivalent: pass the nested dict directly — it is flattened the same way
docs = proj.getSimulationsDocuments(
    type="DispersionRun",
    scenario={"source": "factory_A", "wind": {"direction": 270}},
    status="completed"
)
```

### Loading data from query results

Once you have documents from a query, call `getData()` to load the actual data:

```python
docs = proj.getMeasurementsDocuments(type="WeatherStation", location="Haifa")

for doc in docs:
    print(f"Station: {doc['desc']['station']}, Format: {doc['dataFormat']}")
    df = doc.getData()  # returns pandas.DataFrame, xarray.Dataset, etc.
    print(df.head())
```

### Building queries from dictionaries: `dictToMongoQuery`

Sometimes you have a query as a Python dictionary — loaded from a JSON config file or built programmatically. The utility function `dictToMongoQuery` converts a nested dictionary into MongoEngine's `__` query format:

```python
from hera.utils import dictToMongoQuery

# A nested query dict (e.g., loaded from JSON)
query_dict = {
    "scenario": {
        "source": "factory_A",
        "wind": {
            "speed": 5.0,
            "direction": 270
        }
    },
    "status": "completed"
}

# Convert to MongoEngine query format
mongo_query = dictToMongoQuery(query_dict)
# Result:
# {
#     "scenario__source": "factory_A",
#     "scenario__wind__speed": 5.0,
#     "scenario__wind__direction": 270,
#     "status": "completed"
# }

# Use it directly in a query
docs = proj.getSimulationsDocuments(type="DispersionRun", **mongo_query)
```

The conversion rules:

| Input | Output |
|-------|--------|
| `{"field": "value"}` | `{"field": "value"}` |
| `{"field": {"sub": 1}}` | `{"field__sub": 1}` |
| `{"a": {"b": {"c": 3}}}` | `{"a__b__c": 3}` |
| `{"items": [10, 20]}` | `{"items__0": 10, "items__1": 20}` |

The query methods already apply this conversion to their `desc` keyword arguments, so `dictToMongoQuery` is mainly useful when you need the flat form yourself (for example to inspect it, or to query a collection directly). It is especially useful when your query parameters come from an external source (a JSON file, CLI arguments, or a repository configuration) rather than hardcoded keyword arguments.

Hera also provides `ConfigurationToJSON` (`hera.utils.jsonutils`) which handles the reverse direction — converting Python objects (including pint/unum physical quantities) into JSON-safe dictionaries suitable for storing in `desc` fields. Together, these utilities form the bridge between structured Python data and MongoDB queries.

---

## Working with data sources {#data-sources}

While the Project API gives you low-level control over documents, **data sources** are the recommended way to work with external data through toolkits. A data source is a versioned, named dataset that a toolkit manages for you.

### What is a data source?

Instead of remembering query filters to find your data (`type="ToolkitDataSource", toolkit="MeteoLowFreq", datasourceName="YAVNEEL"`), you simply ask the toolkit by name:

```python
from hera import toolkitHome

meteo = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_LOWFREQ, projectName="WindStudy")

# Get data by name — no need to know how it's stored
df = meteo.getDataSourceData("YAVNEEL")
```

Behind the scenes, data sources are stored as measurement documents with `type="ToolkitDataSource"`. The toolkit handles all the querying and loading for you.

### Listing available data sources

```python
# List data source names
meteo.getDataSourceList()
# ['YAVNEEL', 'HAIFA_PORT', 'BET_DAGAN']

# View as a table with full details (name, version, format, resource path)
meteo.getDataSourceTable()
```

From the CLI:

```bash
# List all data sources across all toolkits
hera-project project measurements list --project WindStudy --shortcut ds

# Filter by name
hera-project project measurements list --project WindStudy --shortcut ds --contains YAVNEEL
```

### Versions

Each data source can have multiple versions, identified by a 3-tuple `(major, minor, patch)`:

```python
# Add a data source with a specific version
meteo.addDataSource(
    dataSourceName="YAVNEEL",
    resource="/data/meteo/yavneel_v1.parquet",
    dataFormat=meteo.datatypes.PARQUET,
    version=(1, 0, 0)
)

# Add a newer version of the same data source
meteo.addDataSource(
    dataSourceName="YAVNEEL",
    resource="/data/meteo/yavneel_v2.parquet",
    dataFormat=meteo.datatypes.PARQUET,
    version=(2, 0, 0)
)
```

### Getting data from a specific version

```python
# Get a specific version
df_v1 = meteo.getDataSourceData("YAVNEEL", version=(1, 0, 0))
df_v2 = meteo.getDataSourceData("YAVNEEL", version=(2, 0, 0))

# Get without specifying version — returns the default or latest
df = meteo.getDataSourceData("YAVNEEL")
```

When no version is specified, Hera resolves the data source in this order:

1. **Default version** — if one is stored in the project config (`<name>_defaultVersion`, set via `setDataSourceDefaultVersion`)
2. **Latest version** — the highest version number among all versions. This fallback is **not** saved as the default, so adding a newer version changes what you get. Pin a version with `setDataSourceDefaultVersion` when you need stable results.

### Setting a default version

```python
# Pin a specific version as the default for this project
# (raises ValueError if that name/version does not exist)
meteo.setDataSourceDefaultVersion("YAVNEEL", version=(1, 0, 0))

# Now getDataSourceData("YAVNEEL") returns version (1, 0, 0)
# regardless of newer versions existing
```

This is useful when you want reproducible results — pin the version your analysis depends on, and new data uploads won't change your output.

### Updating and deleting data sources

```python
# Overwrite an existing version
meteo.addDataSource(
    dataSourceName="YAVNEEL",
    resource="/data/meteo/yavneel_v1_fixed.parquet",
    dataFormat=meteo.datatypes.PARQUET,
    version=(1, 0, 0),
    overwrite=True  # required when the version already exists
)

# Delete a data source version (without version: the default/latest one)
meteo.deleteDataSource("YAVNEEL", version=(1, 0, 0))
```

`addDataSource` raises `ValueError` if the name/version already exists and `overwrite` is `False`. The version defaults to `(0, 0, 1)`, and any extra keyword arguments are stored in the document's `desc`.

### Loading data sources from a repository

Instead of adding data sources one by one, you can define them in a repository JSON and load them all at once. See [Repositories](concepts.md#repositories) in Key Concepts for details.

### Data source methods summary

| Method | Description |
|--------|-------------|
| `getDataSourceData(datasourceName, version=None)` | Load and return the data (`None` if not found) |
| `getDataSourceList(**filters)` | List data source names |
| `getDataSourceMap(**filters)` | List of dicts (format, resource, desc fields) for every data source version |
| `getDataSourceTable(**filters)` | DataFrame of all data sources with metadata |
| `getDataSourceDocument(datasourceName, version=None)` | Get the raw MongoDB document (`None` if not found) |
| `getDataSourceDocuments(datasourceName, version=None)` | Same, as a list (empty if not found) |
| `addDataSource(dataSourceName, resource, dataFormat, version=(0, 0, 1), overwrite=False, **kwargs)` | Register a new data source |
| `deleteDataSource(datasourceName, version=None)` | Remove a data source version |
| `setDataSourceDefaultVersion(name, version)` | Pin a default version |
