# Data Layer

This page covers the MongoDB document model, the `datatypes` system for format dispatch, and the repository JSON pipeline that loads data into projects.

---

## MongoDB Document Model

Hera stores all metadata in MongoDB using a single base model (`Metadata`) with three subtypes. Each document represents a pointer to data — the actual data lives on disk (or inline for small values).

![Diagram](../images/diagrams/architecture_data_layer_0_a442aadb.svg)

<!-- mermaid source (for editing, paste into mermaid.live):
```mermaid
classDiagram
    class Metadata {
        +str projectName
        +str type
        +str resource
        +str dataFormat
        +dict desc
        +getData() object
        +delete()
        +save()
        +asDict() dict
    }

    class Measurements {
        <<"Metadata.Measurements">>
        Observational and toolkit\ndatasource documents
    }

    class Simulations {
        <<"Metadata.Simulations">>
        Simulation model\noutput documents
    }

    class Cache {
        <<"Metadata.Cache">>
        Intermediate results\nand configuration documents
    }

    class AbstractCollection {
        +getDocuments(projectName, **filters) list
        +addDocument(projectName, resource, dataFormat, type, desc)
        +deleteDocuments(projectName, **filters)
    }

    class Measurements_Collection {
        Wraps Metadata.Measurements
    }

    class Simulations_Collection {
        Wraps Metadata.Simulations
    }

    class Cache_Collection {
        Wraps Metadata.Cache
    }

    Metadata <|-- Measurements
    Metadata <|-- Simulations
    Metadata <|-- Cache
    AbstractCollection <|-- Measurements_Collection
    AbstractCollection <|-- Simulations_Collection
    AbstractCollection <|-- Cache_Collection
    Measurements_Collection ..> Measurements : "manages"
    Simulations_Collection ..> Simulations : "manages"
    Cache_Collection ..> Cache : "manages"
```
--> Measurements : "manages"
    Simulations_Collection ..> Simulations : "manages"
    Cache_Collection ..> Cache : "manages"
```
-->
--> Measurements : "manages"
    Simulations_Collection ..> Simulations : "manages"
    Cache_Collection ..> Cache : "manages"
```
-->
-->

### Document Fields

| Field | Type | Description |
|-------|------|-------------|
| `projectName` | `str` | The project this document belongs to |
| `_cls` | `str` | Discriminator: `"Metadata.Measurements"`, `"Metadata.Simulations"`, or `"Metadata.Cache"` |
| `type` | `str` | Application-defined type tag (e.g., `"ToolkitDataSource"`, `"Experiment_rawData"`) |
| `resource` | `str` | Path to the data file on disk, or inline value for small data |
| `dataFormat` | `str` | One of the `datatypes` constants (see below) |
| `desc` | `dict` | Free-form metadata dictionary — toolkit name, version, parameters, etc. |

### Collection Architecture

Each collection type wraps a MongoEngine document class and provides the standard CRUD interface:

![Diagram](../images/diagrams/architecture_data_layer_1_126c6f18.svg)

<!-- mermaid source (for editing, paste into mermaid.live):
```mermaid
flowchart LR
    subgraph ProjectAPI ["Project API Layer"]
        direction TB
        AddMeas["addMeasurementsDocument(\nresource, dataFormat,\ntype, desc)"]
        GetMeas["getMeasurementsDocuments(\n**filters)"]
        DelMeas["deleteMeasurementsDocuments(\n**filters)"]
    end

    subgraph CollectionLayer ["Measurements_Collection"]
        direction TB
        AddDoc["addDocument(\nprojectName, ...)"]
        GetDoc["getDocuments(\nprojectName, **filters)"]
        DelDoc["deleteDocuments(\nprojectName, **filters)"]
    end

    subgraph MongoLayer ["MongoDB"]
        direction TB
        MetadataCol["metadata collection\n\nFilters by:\n_cls = Metadata.Measurements\nprojectName = ...\n+ user filters"]
    end

    AddMeas --> AddDoc
    GetMeas --> GetDoc
    DelMeas --> DelDoc
    AddDoc --> MetadataCol
    GetDoc --> MetadataCol
    DelDoc --> MetadataCol
```
--> GetMeas --> GetDoc
    DelMeas --> DelDoc
    AddDoc --> MetadataCol
    GetDoc --> MetadataCol
    DelDoc --> MetadataCol
```
-->
--> GetMeas --> GetDoc
    DelMeas --> DelDoc
    AddDoc --> MetadataCol
    GetDoc --> MetadataCol
    DelDoc --> MetadataCol
```
-->
-->

!!! note "Three Parallel APIs"
    The `Project` class exposes identical method sets for all three collection types: `addMeasurementsDocument` / `addSimulationsDocument` / `addCacheDocument`, and similarly for get and delete. Under the hood, each delegates to its own `Collection` instance which filters by the `_cls` discriminator.

---

## The datatypes System

**Source:** `hera/datalayer/datahandler.py` (class `datatypes`)

The `datatypes` class defines all supported data format constants and provides the dispatch logic to read and write data in each format.

### Supported Formats

![Diagram](../images/diagrams/architecture_data_layer_2_6b3cb736.svg)

<!-- mermaid source (for editing, paste into mermaid.live):
```mermaid
flowchart TD
    subgraph TabularFormats ["Tabular Data"]
        direction LR
        CSV["CSV_PANDAS\ncsv_pandas\npandas.read_csv"]
        PARQUET["PARQUET\nparquet\ndask.read_parquet"]
        HDF["HDF\nHDF\npandas.read_hdf"]
        JSON_PD["JSON_PANDAS\nJSON_pandas\npandas.read_json"]
    end

    subgraph SpatialFormats ["Spatial / Raster Data"]
        direction LR
        GEO["GEOPANDAS\ngeopandas\ngeopandas.read_file"]
        JSON_GEO["JSON_GEOPANDAS\nJSON_geopandas\ngeopandas.read_file"]
        GEOTIFF["GEOTIFF\ngeotiff\nrasterio.open"]
        NETCDF["NETCDF_XARRAY\nnetcdf_xarray\nxarray.open_dataset"]
        ZARR["ZARR_XARRAY\nzarr_xarray\nxarray.open_zarr"]
    end

    subgraph SerialFormats ["Serialization"]
        direction LR
        PICKLE["PICKLE\npickle\npickle.load"]
        NUMPY["NUMPY_ARRAY\nnumpy_array\nnumpy.load"]
        NUMPY_DICT["NUMPY_DICT_ARRAY\nnumpy_dict_array\nnumpy.load"]
    end

    subgraph PlainFormats ["Plain / Inline"]
        direction LR
        STRING["STRING\nstring\nReturn as-is"]
        DICT["DICT\ndict\nReturn directly"]
        JSON_D["JSON_DICT\nJSON_dict\njson.load"]
        IMAGE["IMAGE\nimage\nmatplotlib.imread"]
    end

    subgraph DynamicFormats ["Dynamic"]
        direction LR
        CLASS["CLASS\nClass\npydoc.locate + instantiate"]
    end
```
--> subgraph DynamicFormats ["Dynamic"]
        direction LR
        CLASS["CLASS\nClass\npydoc.locate + instantiate"]
    end
```
-->
--> subgraph DynamicFormats ["Dynamic"]
        direction LR
        CLASS["CLASS\nClass\npydoc.locate + instantiate"]
    end
```
-->
-->

| Constant | Value | Description |
|----------|-------|-------------|
| `STRING` | `"string"` | Plain text / path string |
| `CSV_PANDAS` | `"csv_pandas"` | CSV file read via pandas |
| `HDF` | `"HDF"` | HDF5 file |
| `NETCDF_XARRAY` | `"netcdf_xarray"` | NetCDF file read via xarray |
| `ZARR_XARRAY` | `"zarr_xarray"` | Zarr archive read via xarray |
| `JSON_DICT` | `"JSON_dict"` | JSON file parsed to dict |
| `JSON_PANDAS` | `"JSON_pandas"` | JSON file read via pandas |
| `JSON_GEOPANDAS` | `"JSON_geopandas"` | GeoJSON file read via geopandas |
| `GEOPANDAS` | `"geopandas"` | Shapefile / GeoPackage read via geopandas |
| `GEOTIFF` | `"geotiff"` | GeoTIFF raster read via rasterio |
| `PARQUET` | `"parquet"` | Parquet file read via dask/pandas |
| `IMAGE` | `"image"` | Image file read via matplotlib |
| `PICKLE` | `"pickle"` | Python pickle file |
| `DICT` | `"dict"` | Inline dictionary (stored in resource) |
| `NUMPY_ARRAY` | `"numpy_array"` | NumPy .npy/.npz file |
| `NUMPY_DICT_ARRAY` | `"numpy_dict_array"` | Dict of NumPy arrays |
| `CLASS` | `"Class"` | Dynamic Python class (imported at runtime) |

### Format Dispatch Flow

When `document.getData()` is called, the system resolves the handler based on `dataFormat`:

![Diagram](../images/diagrams/architecture_data_layer_3_a857df8c.svg)

<!-- mermaid source (for editing, paste into mermaid.live):
```mermaid
flowchart TD
    Start["doc.getData()"] --> ReadFormat["Read doc.dataFormat\nand doc.resource"]

    ReadFormat --> FormatSwitch{dataFormat\nvalue?}

    FormatSwitch -- "string" --> ReturnString["Return resource\npath string as-is"]
    FormatSwitch -- "parquet" --> ReadParquet["dask.dataframe\n.read_parquet(resource)\n.compute()"]
    FormatSwitch -- "netcdf_xarray" --> ReadNC["xarray\n.open_dataset(resource)"]
    FormatSwitch -- "geopandas" --> ReadGeo["geopandas\n.read_file(resource)"]
    FormatSwitch -- "JSON_dict" --> ReadJSON["json.load(\nopen(resource))"]
    FormatSwitch -- "csv_pandas" --> ReadCSV["pandas\n.read_csv(resource)"]
    FormatSwitch -- "Class" --> LoadClass["pydoc.locate(resource)\nInstantiate class"]
    FormatSwitch -- "pickle" --> ReadPickle["pickle.load(\nopen(resource, rb))"]
    FormatSwitch -- "geotiff" --> ReadTiff["rasterio.open(resource)\nConvert to xarray"]
    FormatSwitch -- "image" --> ReadImg["matplotlib.image\n.imread(resource)"]

    ReadParquet --> Return["Return loaded\ndata object"]
    ReadNC --> Return
    ReadGeo --> Return
    ReadJSON --> Return
    ReadCSV --> Return
    LoadClass --> Return
    ReadPickle --> Return
    ReadTiff --> Return
    ReadImg --> Return
    ReturnString --> Return
```
-->  LoadClass --> Return
    ReadPickle --> Return
    ReadTiff --> Return
    ReadImg --> Return
    ReturnString --> Return
```
-->
-->  LoadClass --> Return
    ReadPickle --> Return
    ReadTiff --> Return
    ReadImg --> Return
    ReturnString --> Return
```
-->
-->

!!! tip "Auto-Detection"
    The `datatypes.getDataFormatName(data)` static method can auto-detect the format from a Python object (DataFrame -> `"parquet"`, xarray.Dataset -> `"netcdf_xarray"`, dict -> `"JSON_dict"`, etc.). This is used by `Project.saveData()` to automatically choose the right format and file extension.

---

## Repository JSON Structure

A **repository JSON** is the standard way to declare and load data into a Hera project. It maps toolkit names to their configuration, datasources, and documents.

### Format

```json
{
    "<ToolkitName>": {
        "Config": {
            "key1": "value1",
            "key2": "value2"
        },
        "Datasource": {
            "<datasource_name>": {
                "isRelativePath": "True",
                "item": {
                    "resource": "relative/path/to/data.parquet",
                    "dataFormat": "parquet",
                    "version": [0, 0, 1],
                    "desc": { ... }
                }
            }
        },
        "Measurements": {
            "<measurement_name>": {
                "isRelativePath": "True",
                "item": {
                    "resource": "relative/path/to/file.shp",
                    "dataFormat": "geopandas",
                    "type": "SomeType",
                    "desc": { ... }
                }
            }
        }
    }
}
```

### Loading Pipeline

![Diagram](../images/diagrams/architecture_data_layer_4_1a682ce1.svg)

<!-- mermaid source (for editing, paste into mermaid.live):
```mermaid
sequenceDiagram
    participant User as User / CLI
    participant DT as dataToolkit
    participant TKHome as ToolkitHome
    participant Toolkit as Concrete Toolkit
    participant DB as MongoDB

    User->>DT: loadAllDatasourcesIn<br/>RepositoryJSONToProject(<br/>projectName, repoJSON,<br/>basedir, overwrite,<br/>auto_register_missing)

    note over DT: Parse JSON and iterate<br/>over each toolkit entry

    loop For each toolkit in JSON
        DT->>TKHome: getToolkit(toolkitName,<br/>projectName)

        alt Toolkit found in registry
            TKHome-->>DT: Return toolkit instance
        else Toolkit not found
            note over DT: If auto_register_missing:<br/>search JSON for classpath<br/>hints and auto-register
            DT->>TKHome: auto_register_and_get(<br/>toolkitName, projectName,<br/>repositoryJSON)
            TKHome-->>DT: Return newly registered<br/>toolkit instance
        end

        note over DT: Process each section<br/>type for this toolkit

        loop For each section in toolkit dict
            alt Config section
                DT->>Toolkit: setConfig(**configDict)
                Toolkit->>DB: Upsert config document
            else DataSource section
                note over DT: Resolve relative paths<br/>using basedir if<br/>isRelativePath is True
                DT->>Toolkit: addDataSource(name,<br/>resource, dataFormat,<br/>version)
                Toolkit->>DB: Insert ToolkitDataSource<br/>document
            else Measurements section
                DT->>Toolkit: addMeasurementsDocument(<br/>resource, dataFormat,<br/>type, desc)
                Toolkit->>DB: Insert Measurements doc
            else Simulations section
                DT->>Toolkit: addSimulationsDocument(<br/>resource, dataFormat,<br/>type, desc)
                Toolkit->>DB: Insert Simulations doc
            else Cache section
                DT->>Toolkit: addCacheDocument(<br/>resource, dataFormat,<br/>type, desc)
                Toolkit->>DB: Insert Cache doc
            else Function section
                DT->>Toolkit: Call named function<br/>with parameters
            end
        end
    end

    DT-->>User: Loading complete
```
-->DT->>Toolkit: Call named function<br/>with parameters
            end
        end
    end

    DT-->>User: Loading complete
```
-->
-->DT->>Toolkit: Call named function<br/>with parameters
            end
        end
    end

    DT-->>User: Loading complete
```
-->
-->

### Path Resolution

Each item in the repository JSON has an `isRelativePath` flag:

- **`"True"`** — The `resource` path is relative to the JSON file's directory. The loader prepends `basedir` to make it absolute.
- **`"False"`** — The `resource` is already an absolute path and is used as-is.

!!! warning "String Booleans"
    The `isRelativePath` field accepts both string `"True"/"False"` and Python booleans `true/false`. The loader checks for both forms. Always be explicit to avoid ambiguity.

### Static Loading (No MongoDB)

For testing or lightweight scripts, `dataToolkit` provides two static methods that work without MongoDB:

```python
from hera.utils.data.toolkit import dataToolkit

# Load and resolve all paths in one call
repo = dataToolkit.loadRepositoryFromPath("/path/to/repository.json")

# Or resolve paths on an already-parsed dict
resolved = dataToolkit.resolveDataSourcePaths(repo_dict, basedir="/data/root")
```

These methods perform a deep copy of the JSON and resolve all relative `resource` paths to absolute, but do **not** insert anything into MongoDB.

---

## ToolkitDataSource Documents

When a datasource is registered via `abstractToolkit.addDataSource()`, it creates a special document:

```python
{
    "projectName": "MY_PROJECT",
    "_cls": "Metadata.Measurements",
    "type": "ToolkitDataSource",
    "resource": "/data/meteorology/YAVNEEL.parquet",
    "dataFormat": "parquet",
    "desc": {
        "toolkit": "MeteoLowFreq",
        "datasourceName": "YAVNEEL",
        "version": [0, 0, 1]
    }
}
```

!!! info "Querying Datasources"
    The `abstractToolkit` methods always filter by `type="ToolkitDataSource"` and `toolkit=self.toolkitName`. This ensures that each toolkit only sees its own datasources, even though all documents share the same MongoDB collection.

### Version Resolution

When `getDataSourceDocument(name)` is called without a version:

![Diagram](../images/diagrams/architecture_data_layer_5_40c0c156.svg)

<!-- mermaid source (for editing, paste into mermaid.live):
```mermaid
flowchart TD
    Start["getDataSourceDocument(\nname, version=None)"] --> CheckDefault{Default version\nset in project\nconfig?}

    CheckDefault -- "Yes" --> QueryDefault["Query MongoDB with\ndefault version from config\nkey: name_defaultVersion"]
    CheckDefault -- "No" --> QueryAll["Query all documents\nmatching name\n(any version)"]

    QueryAll --> CountResults{How many\nresults?}

    CountResults -- "0 results" --> ReturnNone["Return None\n(datasource not found)"]
    CountResults -- "1 result" --> ReturnDoc["Return the\nsingle document"]
    CountResults -- "Multiple results" --> PickMax["Sort by version\ntuple and pick\nhighest version"]

    PickMax --> ReturnDoc
    QueryDefault --> ReturnDoc
```
-->s" --> PickMax["Sort by version\ntuple and pick\nhighest version"]

    PickMax --> ReturnDoc
    QueryDefault --> ReturnDoc
```
-->
-->s" --> PickMax["Sort by version\ntuple and pick\nhighest version"]

    PickMax --> ReturnDoc
    QueryDefault --> ReturnDoc
```
-->
-->

---

## API Reference

::: hera.datalayer.datahandler.datatypes
    options:
      show_root_heading: true
      members:
        - STRING
        - CSV_PANDAS
        - NETCDF_XARRAY
        - JSON_DICT
        - GEOPANDAS
        - PARQUET
        - CLASS

::: hera.utils.data.toolkit.dataToolkit
    options:
      show_root_heading: true
      members:
        - addRepository
        - getRepository
        - loadAllDatasourcesInRepositoryJSONToProject
        - resolveDataSourcePaths
        - loadRepositoryFromPath
