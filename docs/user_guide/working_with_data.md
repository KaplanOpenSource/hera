# Working with Data

This page covers the practical details of storing, querying, and loading data in Hera. For the high-level overview, see [Key Concepts](concepts.md).

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

For common Python objects, `saveData` / `saveCacheData` / `saveMeasurementData` auto-detect the format, save the file to disk, and create the document in one call:

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
| `xarray.DataArray` | `zarr_xarray` | `.zarr` |
| `numpy.ndarray` | `numpy_array` | `.npy` |
| `dict`, `list`, `bytes` | `pickle` | `.pckle` |

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
| **Parquet** | `"parquet"` | `PARQUET` | `pandas.DataFrame` or `dask.DataFrame` | `.parquet` |
| **CSV** | `"csv_pandas"` | `CSV_PANDAS` | `pandas.DataFrame` | `.csv` |
| **HDF5** | `"HDF"` | `HDF` | `pandas.DataFrame` or `dask.DataFrame` | `.hdf` |
| **NetCDF** | `"netcdf_xarray"` | `NETCDF_XARRAY` | `xarray.Dataset` | `.nc` |
| **Zarr** | `"zarr_xarray"` | `ZARR_XARRAY` | `xarray.Dataset` | `.zarr` |
| **GeoPackage** | `"geopandas"` | `GEOPANDAS` | `geopandas.GeoDataFrame` | `.gpkg` |
| **GeoJSON** | `"JSON_geopandas"` | `JSON_GEOPANDAS` | `geopandas.GeoDataFrame` | `.json` |
| **GeoTIFF** | `"geotiff"` | `GEOTIFF` | GDAL dataset | `.tif` |
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

Hera uses [MongoEngine](http://mongoengine.org/) under the hood. When you call `getMeasurementsDocuments()`, keyword arguments are translated into MongoDB queries. The `desc` fields are flattened using double-underscore (`__`) notation — the same convention MongoEngine uses for nested document queries.

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

# Find stations above 100m elevation
docs = proj.getMeasurementsDocuments(type="WeatherStation", elevation=120)
```

Behind the scenes, `location="Haifa"` becomes a MongoDB query on `desc.location`, using MongoEngine's `__` syntax: `desc__location="Haifa"`.

### Structured (nested) metadata queries

When your `desc` has nested dictionaries, the double-underscore notation traverses the hierarchy:

```python
# Store a document with nested metadata
proj.addMeasurementsDocument(
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
docs = proj.getMeasurementsDocuments(type="DispersionRun", **mongo_query)
```

The conversion rules:

| Input | Output |
|-------|--------|
| `{"field": "value"}` | `{"field": "value"}` |
| `{"field": {"sub": 1}}` | `{"field__sub": 1}` |
| `{"a": {"b": {"c": 3}}}` | `{"a__b__c": 3}` |
| `{"items": [10, 20]}` | `{"items__0": 10, "items__1": 20}` |

This is especially useful when your query parameters come from an external source (a JSON file, CLI arguments, or a repository configuration) rather than hardcoded keyword arguments.

Hera also provides `ConfigurationToJSON` (`hera.utils.jsonutils`) which handles the reverse direction — converting Python objects (including physical units via `Unum`) into JSON-safe dictionaries suitable for storing in `desc` fields. Together, these utilities form the bridge between structured Python data and MongoDB queries.
