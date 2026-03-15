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
- `MeteoLowFreq` — Toolkit name (must match registered toolkit)
- `DataSource` — Section type (creates a ToolkitDataSource document)
- `YAVNEEL` — Datasource name (used in `getDataSourceData("YAVNEEL")`)
- `isRelativePath: "True"` — Path is relative to repository JSON file location
- `resource` — Path to the data file
- `dataFormat: "parquet"` — How to read the file
- `version: [0, 0, 1]` — Version tuple (major, minor, patch)
- `desc` — Free-form metadata dictionary

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

Multiple versions of the same datasource:

```json
{
    "MeteoLowFreq": {
        "DataSource": {
            "YAVNEEL": {
                "isRelativePath": "True",
                "item": {
                    "resource": "measurements/meteorology/lowfreqdata/YAVNEEL_v1.parquet",
                    "dataFormat": "parquet",
                    "version": [0, 0, 1],
                    "desc": {
                        "stationName": "YAVNEEL",
                        "processedDate": "2024-01-15"
                    }
                }
            },
            "YAVNEEL": {
                "isRelativePath": "True",
                "item": {
                    "resource": "measurements/meteorology/lowfreqdata/YAVNEEL_v2.parquet",
                    "dataFormat": "parquet",
                    "version": [0, 0, 2],
                    "desc": {
                        "stationName": "YAVNEEL",
                        "processedDate": "2024-03-20",
                        "qualityControl": "passed"
                    }
                }
            }
        }
    }
}
```

**Note:** Both entries have the same datasource name (`"YAVNEEL"`) but different versions. The loader will create two separate documents. Use `setDataSourceDefaultVersion()` to choose which version is returned by default.

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
                "isRelativePath": "True",
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
                "isRelativePath": "True",
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
                "isRelativePath": "True",
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
            "initializeToolkit": {
                "params": {
                    "autoLoadDefaults": true
                }
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
| `Function` | `_handle_Function` | Calls a named function with parameters |

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
# Register the repository
hera-project repository add meteo_project /path/to/repository.json

# Load into a project
hera-project repository load meteo_project MY_PROJECT --overwrite
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

- **Toolkit existence** — Toolkit must be registered or auto-registrable
- **Section names** — Must be one of: `Config`, `DataSource`, `Measurements`, `Simulations`, `Cache`, `Function`
- **Required fields** — `resource`, `dataFormat` are required for DataSource items
- **Path resolution** — Relative paths must be resolvable (directory must exist)
- **Version format** — Must be a list of 3 integers: `[major, minor, patch]`

**Common Errors:**

- `Unknown Handler X` — Section name doesn't match expected handlers
- `Toolkit X not found` — Toolkit not registered and auto-registration failed
- `Source X already exists` — Datasource exists and `overwrite=False`

---

## See Also

- [Repository JSON Schema Reference](../reference/repository_schema.md) — Complete schema documentation
- [Best Practices: Repository Structure](../guides/best_practices.md#repository-structure) — Organization guidelines
- [Data Layer: Repository Pipeline](../architecture/data_layer.md#repository-json-structure) — How repositories are loaded
