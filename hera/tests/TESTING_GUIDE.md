# Hera Test Suite — Documentation & Usage Guide

## Overview

The Hera test infrastructure was migrated from a legacy system (custom JSON runner + unittest) to **native Pytest**.  
All tests now live under `hera/tests/` and can be executed with a single `pytest` command.

---

## Architecture

```
hera/
├── pytest.ini                          # Pytest configuration
├── hera/
│   ├── utils/data/toolkit.py           # dataToolkit (enhanced with direct-load methods)
│   └── tests/
│       ├── conftest.py                 # Shared fixtures, comparison helpers, CLI options
│       ├── test_datalayer.py           # Project CRUD tests
│       ├── test_repository.py          # Repository add/get/load, path resolution
│       ├── test_topography.py          # TopographyToolkit tests
│       ├── test_landcover.py           # LandCoverToolkit tests
│       ├── test_lowfreq.py             # lowFreqToolKit + analysis + presentation
│       ├── test_highfreq.py            # HighFreqToolKit + calculators + turbulence
│       ├── test_demography.py          # DemographyToolkit tests
│       ├── UNIT_TEST_DYNAMIC_TOOLKITS/ # Dynamic toolkit tests (kept as-is)
│       ├── repository/testCases/       # Test JSON data for repository tests
│       └── datalayer/testCases/        # Test JSON data for datalayer tests
└── ~/hera_unittest_data/               # External test data repository
    ├── data_config.json                # Data configuration metadata
    ├── measurements/                   # Raw test data files
    │   ├── GIS/raster/                 # HGT, TIF, CSV, NC files
    │   ├── GIS/vector/                 # SHP, GeoJSON files
    │   └── meteorology/               # Parquet files (low/high freq)
    └── expected/                       # Expected output result sets
        ├── BASELINE/
        ├── REGRESSION_20251113_1556/
        └── demo/
```

---

## Prerequisites

### 1. Python Environment

```bash
cd /home/ilay/hera
source heraenv/bin/activate
pip install pytest   # if not already installed
```

### 2. Test Data Repository

The tests rely on external data files stored in `~/hera_unittest_data/`.  
This directory must contain:

- `data_config.json` — metadata about paths, assets, and result sets
- `measurements/` — raw data files (HGT, TIF, SHP, Parquet, etc.)
- `expected/` — expected output files organized by result set

### 3. Environment Variables

| Variable | Required | Default | Description |
|---|---|---|---|
| `TEST_HERA` | No | `~/hera_unittest_data` | Path to the test data repository root |
| `RESULT_SET` | No | `BASELINE` | Name of the expected-output result set |
| `PREPARE_EXPECTED_OUTPUT` | No | (unset) | Set to `"1"` to generate expected outputs instead of comparing |

---

## How to Run Tests

### Run All Tests

```bash
cd /home/ilay/hera
export TEST_HERA=~/hera_unittest_data
pytest hera/tests/ -v
```

### Run a Specific Test Module

```bash
pytest hera/tests/test_datalayer.py -v
pytest hera/tests/test_repository.py -v
pytest hera/tests/test_topography.py -v
pytest hera/tests/test_landcover.py -v
pytest hera/tests/test_lowfreq.py -v
pytest hera/tests/test_highfreq.py -v
pytest hera/tests/test_demography.py -v
```

### Run a Specific Test Class or Function

```bash
# Run all tests in a class
pytest hera/tests/test_topography.py::TestGetPointElevation -v

# Run a single test
pytest hera/tests/test_topography.py::TestGetPointElevation::test_basic -v
```

### Choose a Result Set

```bash
# Via CLI option
pytest hera/tests/ --result-set BASELINE -v

# Via environment variable
export RESULT_SET=REGRESSION_20251113_1556
pytest hera/tests/ -v
```

### Run with Short Traceback

```bash
pytest hera/tests/ -v --tb=short
```

### Run Only Fast Tests (skip slow)

```bash
pytest hera/tests/ -v -m "not slow"
```

### Run with Parallel Workers (requires pytest-xdist)

```bash
pip install pytest-xdist
pytest hera/tests/ -v -n auto
```

---

## Test Modules — Detailed Description

### `test_datalayer.py`

Tests for `hera.datalayer.project.Project` CRUD operations.

| Test | Description |
|---|---|
| `test_project_init` | Verify Project creation and basic properties |
| `test_add_measurements_document` | Add a document, verify it persists |
| `test_get_measurements_documents` | Query documents by resource/format/type |
| `test_delete_measurements_documents` | Delete all documents, verify removal |
| `test_add_and_read_counters` | Read/write Counter documents via setConfig/getConfig |

**Requires:** MongoDB connection

---

### `test_repository.py`

Tests for `hera.utils.data.toolkit.dataToolkit` (repository management).

| Test | Description |
|---|---|
| `test_add_repository` | Register a repository JSON via `addRepository` |
| `test_get_repository` | Retrieve and verify loaded JSON content |
| `test_load_datasources_to_project` | Full round-trip: load repository JSON, assert correct document count |
| `test_resolve_relative_paths` | Verify `isRelativePath` handling produces absolute paths |
| `test_absolute_paths_unchanged` | Verify absolute paths are not modified |
| `test_load_repository_from_path` | Test the direct-load method (no MongoDB) |
| `test_load_repository_nonexistent` | Verify FileNotFoundError for missing files |

**Requires:** MongoDB connection (for add/get/load tests), test JSON in `repository/testCases/`

---

### `test_topography.py`

Tests for `hera.measurements.GIS.raster.topography.TopographyToolkit`.

| Test | Description |
|---|---|
| `test_basic` (getPointElevation) | Single point elevation lookup |
| `test_second_file` | Elevation from a different HGT tile |
| `test_matches_hgt_file` | Verify toolkit result matches raw HGT binary read |
| `test_basic` (getPointListElevation) | Elevation for multiple points |
| `test_matches_hgt_files` | Multi-point comparison against raw HGT data |
| `test_basic` (getElevationOfXarray) | Elevation grid via xarray Dataset |
| `test_matches_hgt_file` (xarray) | Xarray grid comparison against raw HGT data |
| `test_basic` (getElevation) | Area elevation via bounding box |
| `test_matches_hgt_file` (area) | Area elevation comparison against raw HGT data |
| `test_basic` (convertPointsCRS) | CRS conversion (WGS84 -> ITM) |
| `test_basic` (createElevationSTL) | STL string generation |
| `test_basic` (getElevationSTL) | STL from existing Dataset |
| `test_basic` (calculateStatistics) | Mean, min, max statistics |

**Requires:** HGT files in `~/hera_unittest_data/measurements/GIS/raster/`

---

### `test_landcover.py`

Tests for `hera.measurements.GIS.raster.landcover.LandCoverToolkit`.

| Test | Description |
|---|---|
| `test_basic` (getLandCoverAtPoint) | Land cover value at a single point |
| `test_against_raster` | Compare toolkit result with raw rasterio read |
| `test_basic` (getLandCover) | Land cover map for a bounding box |
| `test_map_vs_raster` | Sampled map values vs. raster file |
| `test_at_point` (roughness) | Roughness at a point |
| `test_area` (roughness) | Roughness map for a bounding box |
| `test_values_in_range` | Verify roughness values are within expected range |
| `test_roughnesslength2sandgrainroughness` | Conversion function |
| `test_known_landcover` | Known land cover value -> expected roughness |
| `test_out_of_bounds` | IndexError for out-of-bounds coordinates |
| `test_get_coding_map` | Coding map structure and values |

**Requires:** Landcover TIF in `~/hera_unittest_data/measurements/GIS/raster/`

---

### `test_lowfreq.py`

Tests for `hera.measurements.meteorology.lowfreqdata.toolkit.lowFreqToolKit`, analysis, and presentation layers.

| Category | Tests |
|---|---|
| Toolkit Init | `test_has_analysis`, `test_has_presentation`, `test_has_docType`, `test_docType_value` |
| Analysis | `test_basic` (addDatesColumns), `test_max_normalized`, `test_density`, `test_y_normalized_behaviour`, `test_basic` (resampleSecondMoments) |
| Presentation | `test_plotScatter`, `test_dateLinePlot`, `test_plotProbContourf`, `test_plotProbContourf_bySeason` |
| Data Matching | `test_dateLinePlot_matches_data`, `test_plotScatter_matches_data` |
| Edge Cases | `test_scatter_empty_dataframe`, `test_scatter_nan_and_outliers`, `test_scatter_WS_field` |
| Distribution | `test_contourf_distribution_ranges` |
| Save | `test_scatter_creates_non_empty_image` |

**Requires:** `YAVNEEL.parquet` in `~/hera_unittest_data/measurements/meteorology/lowfreqdata/`

---

### `test_highfreq.py`

Tests for `hera.measurements.meteorology.highfreqdata` toolkit, analysis calculators, and turbulence statistics.

| Category | Tests |
|---|---|
| Toolkit | `test_docType_property` |
| Data Reading | `test_read_sonic_data`, `test_read_trh_data`, `test_read_nonexistent_file` |
| Time Range | `test_sonic_time_range`, `test_trh_time_range` |
| Specific Points | `test_sonic_first_row`, `test_trh_first_row` |
| Error Paths | `test_campbelToParquet_nonexistent`, `test_asciiToParquet_nonexistent` |
| AbstractCalculator | `test_init_basic`, `test_sampling_window`, `test_compute_methods_exist`, `test_set_save_properties` |
| MeanDataCalculator | `test_calculate_mean`, `test_hour_and_timeWithinDay`, `test_horizontalSpeed`, `test_sigma_sigmaH`, `test_Ustar_and_uStarOverWindSpeed`, `test_compute_returns_dataframe` |
| Advanced MeanData | `test_TKE`, `test_MOLength` |
| RawdataAnalysis | `test_singlePointTurbulenceStatistics_returns_instance`, `test_raises_on_invalid`, `test_AveragingCalculator`, `test_AveragingCalculator_raises_on_invalid` |
| Turbulence Stats | `test_instantiation`, `test_invalid_input_type`, `test_fluctuations`, `test_secondMoments`, `test_sigma`, `test_horizontalSpeed`, `test_Ustar`, `test_TKE`, `test_MOLength_Sonic` |

**Requires:** Sonic/TRH parquet files in `~/hera_unittest_data/measurements/meteorology/highfreqdata/`

---

### `test_demography.py`

Tests for `hera.measurements.GIS.vector.demography.DemographyToolkit`.

| Test | Description |
|---|---|
| `test_basic` (calculatePopulationInPolygon) | Basic polygon intersection |
| `test_partial_intersection` | Partial polygon overlap |
| `test_outside_bounds` | Polygon completely outside data extent |
| `test_invalid_datasource` | ValueError for non-existing data source |
| `test_with_known_values` | Synthetic data with known population values |
| `test_simple` (createNewArea) | Create new area and verify total population |
| `test_creates_and_sets_path` (setDefaultDirectory) | Directory creation and path assignment |

**Requires:** `population_lamas.shp` in `~/hera_unittest_data/measurements/GIS/vector/`

---

## Shared Fixtures (conftest.py)

| Fixture | Scope | Description |
|---|---|---|
| `test_hera_root` | session | Validated path to `~/hera_unittest_data` |
| `data_config` | session | Parsed `data_config.json` dict |
| `result_set` | session | Active result-set name |
| `expected_dir` | session | Path to `expected/<result_set>/` |
| `gis_raster_path` | session | Path to GIS raster data |
| `gis_vector_path` | session | Path to GIS vector data |
| `lowfreq_path` | session | Path to low-frequency meteorology data |
| `highfreq_path` | session | Path to high-frequency meteorology data |
| `project_fixture` | function | Temporary Project with cleanup |
| `data_toolkit_fixture` | session | dataToolkit instance |

---

## Comparison Helpers

Available in `conftest.py` for use in tests:

```python
from hera.tests.conftest import compare_dataframes, compare_dataarrays, compare_outputs

# DataFrame comparison with numeric tolerance
assert compare_dataframes(result_df, expected_df, rtol=1e-6, atol=1e-6)

# DataArray comparison
assert compare_dataarrays(result_da, expected_da)

# Type-based comparison (supports: dataframe, geodataframe, xarray, float, dict, etc.)
assert compare_outputs(result, expected, "dataframe")
```

---

## New dataToolkit Methods

Two methods were added to `hera.utils.data.toolkit.dataToolkit` to support direct loading without MongoDB:

### `loadRepositoryFromPath(json_path)` (static)

```python
from hera.utils.data.toolkit import dataToolkit

repo = dataToolkit.loadRepositoryFromPath("/path/to/repository.json")
# Returns dict with all relative resource paths resolved to absolute
```

### `resolveDataSourcePaths(repositoryJSON, basedir)` (static)

```python
resolved = dataToolkit.resolveDataSourcePaths(repo_dict, basedir="/data/root")
# Deep-copies the dict and resolves all relative resource paths
```

---

## Troubleshooting

### Tests are skipped

- **"TEST_HERA directory not found"** — Set `TEST_HERA` env var or create `~/hera_unittest_data/`
- **"No .hgt files found"** — Ensure HGT files exist under `measurements/GIS/raster/`
- **"YAVNEEL.parquet not found"** — Ensure parquet files exist under `measurements/meteorology/lowfreqdata/`

### MongoDB connection errors

Tests in `test_datalayer.py` and `test_repository.py` require an active MongoDB instance.  
Ensure MongoDB is running and accessible before running those tests.

### Matplotlib backend issues

Presentation tests (plots) may require a non-interactive backend:

```bash
export MPLBACKEND=Agg
pytest hera/tests/test_lowfreq.py -v
```

---

## Migration Summary

### What was removed

| File/Directory | Replaced by |
|---|---|
| `hera/tests/run_all_definitions.py` | Native pytest modules |
| `hera/tests/run_all_json_tests.sh` | `pytest` CLI |
| `hera/tests/env.template` | `conftest.py` fixtures |
| `hera/tests/json_definitions/` (8 files) | Individual `test_*.py` modules |
| `hera/tests/datalayer/project.py` | `test_datalayer.py` |
| `hera/tests/repository/repository.py` | `test_repository.py` |
| `hera/measurements/**/test_unit_*.py` (5 files) | `test_topography.py`, `test_landcover.py`, `test_lowfreq.py`, `test_highfreq.py`, `test_demography.py` |

### What was kept

- `hera/tests/UNIT_TEST_DYNAMIC_TOOLKITS/` — already native pytest, kept as-is
- `hera/tests/repository/testCases/` — test data JSONs used by `test_repository.py`
- `hera/tests/datalayer/testCases/` — test data JSONs used by `test_datalayer.py`
- `hera/tests/DEMO/` — demo data and repositories
