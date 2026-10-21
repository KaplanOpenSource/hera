# Repository JSON Examples

Complete, annotated examples of repository JSON files for different use cases.

---

## Minimal Example

A single toolkit with one datasource:

```json
{
    "MeteoLowFreq": {
        "DataSource": {
            "YAVNEEL": {
                "isRelativePath": "True",
                "item": {
                    "resource": "measurements/meteorology/lowfreqdata/YAVNEEL.parquet",
                    "dataFormat": "parquet",
                    "version": [0, 0, 1],
                    "desc": {
                        "stationName": "YAVNEEL",
                        "type": "lowfreq"
                    }
                }
            }
        }
    }
}
```

**Explanation:**
- `MeteoLowFreq` — Toolkit name (must be resolvable by `toolkitHome.getToolkit()`; unknown names are silently skipped)
- `DataSource` — Section type (creates a ToolkitDataSource document)
- `YAVNEEL` — Datasource name (used in `getDataSourceData("YAVNEEL")`)
- `isRelativePath: "True"` — Path is relative to repository JSON file location
- `resource` — Path to the data file
- `dataFormat: "parquet"` — How to read the file
- `version: [0, 0, 1]` — Version tuple (major, minor, patch)
- `desc` — Free-form metadata. Every key of `item` other than `resource`, `dataFormat`, `version` (and `overwrite`) is stored in the document's `desc`, so this dictionary ends up as `doc.desc["desc"]`. Put keys directly in `item` instead if you want them at the top level of `desc`.

---

## Multi-Toolkit Example

Multiple toolkits with different section types:

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
                    "resource": "measurements/GIS/raster/",
                    "dataFormat": "string",
                    "version": [0, 0, 1],
                    "desc": {
                        "defaultSRTM": "SRTMGL1"
                    }
                }
            }
        }
    },
    "MeteoLowFreq": {
        "DataSource": {
            "YAVNEEL": {
                "isRelativePath": "True",
                "item": {
                    "resource": "measurements/meteorology/lowfreqdata/YAVNEEL.parquet",
                    "dataFormat": "parquet",
                    "version": [0, 0, 1],
                    "desc": {
                        "stationName": "YAVNEEL",
                        "type": "lowfreq"
                    }
                }
            }
        }
    },
    "GIS_Demography": {
        "DataSource": {
            "lamas_population": {
                "isRelativePath": "True",
                "item": {
                    "resource": "measurements/GIS/vector/population_lamas.shp",
                    "dataFormat": "geopandas",
                    "version": [0, 0, 1],
                    "desc": {
                        "source": "LAMAS",
                        "year": 2020
                    }
                }
            }
        }
    }
}
```

**Key Points:**
- Each toolkit section is independent
- `Config` section sets toolkit configuration (via `toolkit.setConfig()`)
- Different data formats: `string`, `parquet`, `geopandas`
- All paths are relative (`isRelativePath: "True"`)

---

## Versioned Datasources

A repository JSON can hold **only one version per datasource name**: JSON objects cannot contain the same key twice (most parsers silently keep the last one), and the loader looks datasources up by name only. When a datasource with that name already exists, it is skipped unless `overwrite` is set, in which case the existing documents of that name are deleted before the new item is added.

Declare the version you want the repository to provide:

```json
{
    "MeteoLowFreq": {
        "DataSource": {
            "YAVNEEL": {
                "isRelativePath": "True",
                "item": {
                    "resource": "measurements/meteorology/lowfreqdata/YAVNEEL_v2.parquet",
                    "dataFormat": "parquet",
                    "version": [0, 0, 2],
                    "processedDate": "2024-03-20",
                    "qualityControl": "passed"
                }
            }
        }
    }
}
```

To keep several versions side by side, add them through the API (or use distinct datasource names):

```python
lf = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_LOWFREQ, projectName="MY_PROJECT")
lf.addDataSource("YAVNEEL", "/data/YAVNEEL_v1.parquet", "parquet", version=(0, 0, 1))
lf.addDataSource("YAVNEEL", "/data/YAVNEEL_v2.parquet", "parquet", version=(0, 0, 2))
lf.setDataSourceDefaultVersion("YAVNEEL", (0, 0, 2))
```

---

## All Section Types

Complete example showing all supported section types:

```json
{
    "MeteoLowFreq": {
        "Config": {
            "defaultStation": "YAVNEEL",
            "timezone": "Asia/Jerusalem"
        },
        "DataSource": {
            "YAVNEEL": {
                "isRelativePath": "True",
                "item": {
                    "resource": "measurements/meteorology/lowfreqdata/YAVNEEL.parquet",
                    "dataFormat": "parquet",
                    "version": [0, 0, 1],
                    "desc": {
                        "stationName": "YAVNEEL"
                    }
                }
            }
        },
        "Measurements": {
            "raw_export_2024": {
                "item": {
                    "resource": "measurements/meteorology/exports/raw_2024.parquet",
                    "dataFormat": "parquet",
                    "type": "Experiment_rawData",
                    "desc": {
                        "exportDate": "2024-12-01",
                        "source": "IMS"
                    }
                }
            }
        },
        "Simulations": {
            "wind_simulation_001": {
                "item": {
                    "resource": "simulations/wind/wind_001.nc",
                    "dataFormat": "netcdf_xarray",
                    "type": "WindProfile",
                    "desc": {
                        "simulationDate": "2024-11-15",
                        "solver": "simpleFoam"
                    }
                }
            }
        },
        "Cache": {
            "processed_stats": {
                "item": {
                    "resource": "cache/statistics.json",
                    "dataFormat": "JSON_dict",
                    "type": "ProcessedStatistics",
                    "desc": {
                        "computedDate": "2024-11-20",
                        "method": "hourly_distribution"
                    }
                }
            }
        },
        "Function": {
            "addRepository": {
                "repositoryName": "shared_gis",
                "repositoryPath": "/data/repositories/shared_gis.json"
            }
        }
    }
}
```

**Section Types Explained:**

| Section | Handler | Action |
|---------|---------|--------|
| `Config` | `_handle_Config` | Calls `toolkit.setConfig(**values)` |
| `DataSource` | `_handle_DataSource` | Calls `toolkit.addDataSource(...)` |
| `Measurements` | `_DocumentHandler` | Calls `toolkit.addMeasurementsDocument(...)` |
| `Simulations` | `_DocumentHandler` | Calls `toolkit.addSimulationsDocument(...)` |
| `Cache` | `_DocumentHandler` | Calls `toolkit.addCacheDocument(...)` |
| `Function` | `_handle_Function` | Calls `dataToolkit.<name>(**params, overwrite=...)` (a method of the loader, not of the toolkit) |

Measurements/Simulations/Cache `item`s may contain only `resource`, `dataFormat`, `type` and `desc`; their `resource` is always joined to the repository directory (absolute paths are kept), so they need no `isRelativePath`.

---

## Relative vs Absolute Paths

Example showing both path resolution patterns:

```json
{
    "GIS_Raster_Topography": {
        "DataSource": {
            "SRTMGL1_relative": {
                "isRelativePath": "True",
                "item": {
                    "resource": "measurements/GIS/raster/",
                    "dataFormat": "string",
                    "version": [0, 0, 1],
                    "desc": {}
                }
            },
            "SRTMGL1_absolute": {
                "isRelativePath": "False",
                "item": {
                    "resource": "/data/shared/GIS/SRTM/",
                    "dataFormat": "string",
                    "version": [0, 0, 1],
                    "desc": {}
                }
            }
        }
    }
}
```

**Path Resolution:**

If the repository JSON is at `/home/user/repos/my_repo.json`:

- **Relative path** (`isRelativePath: "True"`):
  - `resource: "measurements/GIS/raster/"`
  - Resolved to: `/home/user/repos/measurements/GIS/raster/`

- **Absolute path** (`isRelativePath: "False"`):
  - `resource: "/data/shared/GIS/SRTM/"`
  - Used as-is: `/data/shared/GIS/SRTM/`

!!! tip "Best Practice"
    Use relative paths (`isRelativePath: "True"`) when the repository JSON and data files are in the same directory tree. This makes the repository portable and easier to share.

---

## Real-World Example: Complete Project Setup

A realistic repository for a meteorological analysis project:

```json
{
    "GIS_Raster_Topography": {
        "Config": {
            "defaultSRTM": "SRTMGL1",
            "defaultCRS": 4326
        },
        "DataSource": {
            "SRTMGL1": {
                "isRelativePath": "True",
                "item": {
                    "resource": "measurements/GIS/raster/SRTM/",
                    "dataFormat": "string",
                    "version": [0, 0, 1],
                    "desc": {
                        "resolution": "30m",
                        "source": "NASA"
                    }
                }
            }
        }
    },
    "GIS_LandCover": {
        "Config": {
            "defaultLandCover": "lc_mcd12q1"
        },
        "DataSource": {
            "lc_mcd12q1": {
                "isRelativePath": "True",
                "item": {
                    "resource": "measurements/GIS/raster/lc_mcd12q1.tif",
                    "dataFormat": "string",
                    "version": [0, 0, 1],
                    "desc": {
                        "year": 2020,
                        "source": "MODIS"
                    }
                }
            }
        }
    },
    "MeteoLowFreq": {
        "DataSource": {
            "YAVNEEL": {
                "isRelativePath": "True",
                "item": {
                    "resource": "measurements/meteorology/lowfreqdata/YAVNEEL.parquet",
                    "dataFormat": "parquet",
                    "version": [0, 0, 1],
                    "desc": {
                        "stationName": "YAVNEEL",
                        "type": "lowfreq",
                        "latitude": 31.7683,
                        "longitude": 35.2137
                    }
                }
            },
            "TEL_AVIV": {
                "isRelativePath": "True",
                "item": {
                    "resource": "measurements/meteorology/lowfreqdata/TEL_AVIV.parquet",
                    "dataFormat": "parquet",
                    "version": [0, 0, 1],
                    "desc": {
                        "stationName": "TEL_AVIV",
                        "type": "lowfreq",
                        "latitude": 32.0853,
                        "longitude": 34.7818
                    }
                }
            }
        }
    },
    "GIS_Demography": {
        "DataSource": {
            "lamas_population": {
                "isRelativePath": "True",
                "item": {
                    "resource": "measurements/GIS/vector/population_lamas.shp",
                    "dataFormat": "geopandas",
                    "version": [0, 0, 1],
                    "desc": {
                        "source": "LAMAS",
                        "year": 2020,
                        "crs": 2039
                    }
                }
            }
        }
    }
}
```

**Usage:**

```bash
# Register the repository (registered under the file's base name)
hera-project repository add /path/to/meteo_project.json

# Load into a project
hera-project repository load /path/to/meteo_project.json MY_PROJECT --overwrite
```

Or via Python:

```python
from hera.utils.data.toolkit import dataToolkit
import json

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

---

## Common Patterns

In these patterns `MyToolkit` stands for any toolkit name that `toolkitHome.getToolkit()` can resolve (a built-in name or a toolkit registered with `hera-project addToolkit`); otherwise the whole key is skipped.

### Pattern 1: Multiple Data Formats

```json
{
    "MyToolkit": {
        "DataSource": {
            "csv_data": {
                "isRelativePath": "True",
                "item": {
                    "resource": "data/input.csv",
                    "dataFormat": "csv_pandas",
                    "version": [0, 0, 1],
                    "desc": {}
                }
            },
            "netcdf_data": {
                "isRelativePath": "True",
                "item": {
                    "resource": "data/output.nc",
                    "dataFormat": "netcdf_xarray",
                    "version": [0, 0, 1],
                    "desc": {}
                }
            },
            "geojson_data": {
                "isRelativePath": "True",
                "item": {
                    "resource": "data/boundaries.geojson",
                    "dataFormat": "JSON_geopandas",
                    "version": [0, 0, 1],
                    "desc": {}
                }
            }
        }
    }
}
```

### Pattern 2: Inline Configuration

```json
{
    "MyToolkit": {
        "Config": {
            "setting1": "value1",
            "setting2": 42,
            "setting3": {
                "nested": "config"
            }
        }
    }
}
```

The `Config` section is passed directly to `toolkit.setConfig(**configDict)`.

### Pattern 3: String Resources (Directory Paths)

Some toolkits accept directory paths as strings:

```json
{
    "GIS_Raster_Topography": {
        "DataSource": {
            "SRTMGL1": {
                "isRelativePath": "True",
                "item": {
                    "resource": "measurements/GIS/raster/",
                    "dataFormat": "string",
                    "version": [0, 0, 1],
                    "desc": {}
                }
            }
        }
    }
}
```

The toolkit will search this directory for files (e.g., `.hgt` files for SRTM).

---

## Validation and Error Handling

The repository loader validates:

- **Toolkit existence** — Keys that `toolkitHome.getToolkit()` cannot resolve are skipped (INFO log only); no auto-registration is attempted
- **Section names** — After `key.title()`, must be one of: `Config`, `Datasource`, `Measurements`, `Simulations`, `Cache`, `Function` (an unknown section raises `ValueError`)
- **Required fields** — `isRelativePath` on each DataSource entry, and `resource` (or `resourceFilePath`) and `dataFormat` in its `item`
- **Version format** — A list of 3 integers: `[major, minor, patch]`
- Paths are **not** checked at load time; a wrong path fails only when the data is read

**Common Errors:**

- `Unkonw Handler X. The handler must be ...` — Section name doesn't match a handler
- `Skipping key 'X' in repository JSON` (INFO) — Toolkit cannot be resolved
- `Source X already exists in <project>. Use --overwrite to force update` — Datasource exists and `overwrite=False` (logged, item skipped)
- `isRelativePath must be defined as 'True' or 'False'.` — Missing/invalid `isRelativePath` on a DataSource entry

---

## See Also

- [Repository JSON Schema Reference](../reference/repository_schema.md) — Complete schema documentation
- [Best Practices: Repository Structure](../guides/best_practices.md#repository-structure) — Organization guidelines
- [Data Layer: Repository Pipeline](../architecture/data_layer.md#repository-json-structure) — How repositories are loaded
