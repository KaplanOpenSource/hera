# Testing Flow

This page provides a complete technical walkthrough of the Hera test infrastructure — how tests are organized, how data flows from JSON files through MongoDB into toolkit instances, and how results are compared against expected outputs.

---

## Overview

The Hera test suite uses **native Pytest** with a **project-based data access** pattern. The core principle:

!!! tip "The Golden Rule"
    **Tests never access files directly by path.** They interact only with the `Project` and `Toolkit` APIs, exactly as production code does. Data is loaded into MongoDB once per session, and toolkits read it back through the standard datasource mechanism.

### Test Directory Structure

```
hera/tests/
├── conftest.py                      # Session fixtures, comparison helpers
├── test_datalayer.py                # Project CRUD tests
├── test_repository.py               # Repository add/get/load tests
├── test_topography.py               # TopographyToolkit tests
├── test_landcover.py                # LandCoverToolkit tests
├── test_lowfreq.py                  # lowFreqToolKit tests
├── test_highfreq.py                 # HighFreqToolKit tests
├── test_demography.py               # DemographyToolkit tests
├── repository/testCases/            # Test JSON data for repository tests
├── datalayer/testCases/             # Test JSON data for datalayer tests
├── expected/
│   ├── BASELINE/                    # Default expected outputs
│   └── REGRESSION_2025_11_11/       # Alternative result set
└── TESTING_GUIDE.md                 # Human-readable test guide
```

External test data lives in a separate directory:

```
~/hera_unittest_data/                # Configured via TEST_HERA env var
├── data_config.json                 # Data configuration metadata
├── test_repository.json             # Hera-format repository mapping
├── measurements/                    # Raw test data files
│   ├── GIS/raster/                  # HGT, TIF files
│   ├── GIS/vector/                  # SHP files
│   └── meteorology/                 # Parquet files
└── expected/                        # Expected output result sets
    ├── BASELINE/
    └── REGRESSION_20251113_1556/
```

---

## Session Lifecycle — The Complete Flow

```mermaid
sequenceDiagram
    participant Pytest as Pytest Runner
    participant Conftest as conftest.py
    participant RepoJSON as test_repository.json
    participant DT as dataToolkit
    participant TKHome as ToolkitHome
    participant MongoDB as MongoDB
    participant TKFixture as Toolkit Fixtures
    participant TestModule as Test Modules
    participant Expected as expected/ directory

    rect rgb(240, 248, 255)
    note over Pytest, MongoDB: PHASE 1 -- Session Setup (runs once)

    Pytest ->> Conftest: Collect session-scoped fixtures
    Conftest ->> Conftest: test_hera_root fixture:<br/>Read TEST_HERA env var<br/>Default: ~/hera_unittest_data
    Conftest ->> Conftest: data_config fixture:<br/>Parse data_config.json<br/>from test_hera_root
    Conftest ->> Conftest: result_set fixture:<br/>Priority: CLI flag ><br/>env var > config > BASELINE
    Conftest ->> Conftest: expected_dir fixture:<br/>Build path to<br/>expected/RESULT_SET/

    note over Conftest, RepoJSON: Load repository JSON

    Conftest ->> RepoJSON: hera_test_project fixture:<br/>Open test_repository.json
    RepoJSON -->> Conftest: Parsed JSON dictionary<br/>with all toolkit entries

    note over Conftest, MongoDB: Create project and populate data

    Conftest ->> Conftest: Create Project instance<br/>projectName = PYTEST_HERA_PROJECT
    Conftest ->> DT: loadAllDatasourcesIn<br/>RepositoryJSONToProject(<br/>projectName, repoJSON,<br/>basedir, overwrite=True)

    loop For each toolkit in repository JSON
        DT ->> TKHome: getToolkit(toolkitName,<br/>projectName)
        TKHome -->> DT: Return toolkit instance

        note over DT, MongoDB: Process Config, DataSource,<br/>Measurements sections

        DT ->> MongoDB: toolkit.setConfig(**values)
        DT ->> MongoDB: toolkit.addDataSource(<br/>name, resource, format,<br/>version)
    end

    Conftest -->> Pytest: Yield populated<br/>Project instance
    end

    rect rgb(240, 255, 240)
    note over Conftest, TKFixture: PHASE 2 -- Fixture Resolution (runs once per session)

    Conftest ->> TKFixture: topo_toolkit =<br/>TopographyToolkit(<br/>PYTEST_HERA_PROJECT)
    Conftest ->> TKFixture: lc_toolkit =<br/>LandCoverToolkit(<br/>PYTEST_HERA_PROJECT)
    Conftest ->> TKFixture: demo_toolkit =<br/>DemographyToolkit(<br/>PYTEST_HERA_PROJECT)
    Conftest ->> TKFixture: lf_toolkit =<br/>lowFreqToolKit(<br/>PYTEST_HERA_PROJECT)
    Conftest ->> TKFixture: hf_toolkit =<br/>HighFreqToolKit(<br/>PYTEST_HERA_PROJECT)
    end

    rect rgb(255, 248, 235)
    note over TKFixture, Expected: PHASE 3 -- Test Execution (per test function)

    Pytest ->> TestModule: Run test functions in<br/>test_topography.py,<br/>test_lowfreq.py, etc.

    TestModule ->> TKFixture: toolkit.getDataSourceData(<br/>"YAVNEEL")
    TKFixture ->> MongoDB: Query ToolkitDataSource<br/>document by name
    MongoDB -->> TKFixture: Return document with<br/>resource path + dataFormat
    TKFixture -->> TestModule: Return loaded data:<br/>DataFrame / xarray / etc.

    note over TestModule: Call analysis and<br/>presentation methods

    TestModule ->> TestModule: result = toolkit.analysis<br/>.someMethod(data, ...)

    note over TestModule, Expected: Compare with expected output

    TestModule ->> Expected: load_expected_output(<br/>filename, output_type,<br/>expected_dir)
    Expected -->> TestModule: expected data object
    TestModule ->> TestModule: compare_outputs(<br/>result, expected,<br/>output_type)
    TestModule -->> Pytest: PASS or FAIL
    end

    rect rgb(255, 240, 240)
    note over Pytest, MongoDB: PHASE 4 -- Teardown (runs once)

    Pytest ->> Conftest: Session teardown triggered
    Conftest ->> MongoDB: Delete ALL documents<br/>where projectName =<br/>PYTEST_HERA_PROJECT
    MongoDB -->> Conftest: Deletion confirmed
    Conftest -->> Pytest: Cleanup complete
    end
```

---

## Phase 1: Session Setup (conftest.py)

### The hera_test_project Fixture

This is the **single most important fixture** in the test suite. It runs once per session and populates a fresh Hera project with all test data.

```python
# Simplified from hera/tests/conftest.py

@pytest.fixture(scope="session")
def hera_test_project(test_hera_root):
    from hera.datalayer.project import Project
    from hera.utils.data.toolkit import dataToolkit

    # 1. Read the repository JSON
    repo_json_path = test_hera_root / "test_repository.json"
    with open(repo_json_path) as fh:
        repo_json = json.load(fh)

    # 2. Create the project
    proj = Project(projectName="PYTEST_HERA_PROJECT")

    # 3. Load ALL datasources into the project
    dt = dataToolkit()
    dt.loadAllDatasourcesInRepositoryJSONToProject(
        projectName="PYTEST_HERA_PROJECT",
        repositoryJSON=repo_json,
        basedir=str(test_hera_root),
        overwrite=True,
    )

    yield proj

    # 4. Cleanup: remove all documents
    for doc in proj.getMeasurementsDocuments():
        doc.delete()
```

### What loadAllDatasourcesInRepositoryJSONToProject Does

```mermaid
flowchart TD
    Start["Start:\nrepositoryJSON dict\n+ basedir path"] --> IterToolkits["Iterate over each\ntoolkitName key\nin the JSON"]

    IterToolkits --> GetToolkit["toolkitHome.getToolkit(\ntoolkitName, projectName)"]
    GetToolkit --> ToolkitOK{Toolkit\ninstance\nobtained?}

    ToolkitOK -- "Yes" --> IterSections["Iterate over each\nsection key in\ntoolkit sub-dict"]
    ToolkitOK -- "No" --> TryAuto{auto_register_missing\nenabled?}
    TryAuto -- "Yes" --> AutoReg["auto_register_and_get(\ntoolkitName, projectName,\nrepositoryJSON)"]
    AutoReg --> IterSections
    TryAuto -- "No" --> SkipWarn["Log warning and\nskip this toolkit"]
    SkipWarn --> IterToolkits

    IterSections --> SectionType{Section\ntype?}

    SectionType -- "Config" --> HandleConfig["_handle_Config:\ntoolkit.setConfig(\n**configDict)"]
    SectionType -- "DataSource" --> IterDS["Iterate over each\ndatasource item"]
    SectionType -- "Measurements" --> HandleMeas["_DocumentHandler:\ntoolkit.addMeasurements\nDocument(resource,\ndataFormat, type, desc)"]
    SectionType -- "Simulations" --> HandleSim["_DocumentHandler:\ntoolkit.addSimulations\nDocument(resource,\ndataFormat, type, desc)"]
    SectionType -- "Cache" --> HandleCache["_DocumentHandler:\ntoolkit.addCache\nDocument(resource,\ndataFormat, type, desc)"]
    SectionType -- "Function" --> HandleFunc["_handle_Function:\nCall named function\nwith provided parameters"]

    IterDS --> CheckRelPath{isRelativePath\n== True?}
    CheckRelPath -- "Yes" --> ResolvePath["Prepend basedir\nto resource path:\nresource = basedir /\nresource"]
    CheckRelPath -- "No" --> UseAbsPath["Use resource\npath as-is"]
    ResolvePath --> AddDS["toolkit.addDataSource(\nname, resource,\ndataFormat, version,\noverwrite)"]
    UseAbsPath --> AddDS

    HandleConfig --> NextSection["Next section\nor next toolkit"]
    AddDS --> NextSection
    HandleMeas --> NextSection
    HandleSim --> NextSection
    HandleCache --> NextSection
    HandleFunc --> NextSection
    NextSection --> SectionType
```

!!! info "Overwrite Mode"
    The `overwrite=True` parameter ensures that running the test suite multiple times does not accumulate stale documents. Existing documents with the same datasource name are deleted before the new ones are inserted.

---

## Phase 2: Fixture Resolution

### Session-Scoped Toolkit Fixtures

Each toolkit test module depends on a session-scoped fixture that instantiates the real toolkit class, connected to the test project:

| Fixture | Toolkit Class | Depends On | Data Sources |
|---------|---------------|------------|--------------|
| `topo_toolkit` | `TopographyToolkit` | `hera_test_project` | `SRTMGL1` (HGT directory) |
| `lc_toolkit` | `LandCoverToolkit` | `hera_test_project` | `lc_mcd12q1` (TIF path) |
| `demo_toolkit` | `DemographyToolkit` | `hera_test_project` | `lamas_population` (SHP -> GeoDataFrame) |
| `lf_toolkit` | `lowFreqToolKit` | `hera_test_project` | `YAVNEEL` (parquet -> dask/pandas) |
| `hf_toolkit` | `HighFreqToolKit` | `hera_test_project` | `slicedYamim_sonic`, `slicedYamim_TRH` (parquet) |

### Function-Scoped Fixtures

| Fixture | Scope | Description |
|---------|-------|-------------|
| `project_fixture` | function | Temporary `Project` with cleanup (for `test_datalayer.py`) |
| `data_toolkit_fixture` | session | `dataToolkit` instance (for `test_repository.py`) |

### Fixture Dependency Graph

```mermaid
flowchart TD
    subgraph EnvLayer ["Environment Layer"]
        TestHera["test_hera_root\n\nSource: TEST_HERA env var\nDefault: ~/hera_unittest_data\nScope: session"]
    end

    subgraph ConfigLayer ["Configuration Layer"]
        DataCfg["data_config\n\nSource: data_config.json\nContains: default_result_set,\nmetadata\nScope: session"]
        ResultSet["result_set\n\nPriority:\n1. --result-set CLI flag\n2. RESULT_SET env var\n3. data_config default\n4. 'BASELINE' fallback\nScope: session"]
        ExpDir["expected_dir\n\nPath: test_hera_root /\nexpected / result_set\nScope: session"]
    end

    subgraph ProjectLayer ["Project Layer"]
        ProjName["hera_project_name\n\nConstant string:\nPYTEST_HERA_PROJECT\nScope: session"]
        HeraPrj["hera_test_project\n\nFull Project instance\npopulated with all\ntest data from\ntest_repository.json\nScope: session"]
    end

    subgraph ToolkitLayer ["Toolkit Fixtures"]
        direction LR
        TopoTK["topo_toolkit\n\nTopography\nToolkit"]
        LcTK["lc_toolkit\n\nLandCover\nToolkit"]
        DemoTK["demo_toolkit\n\nDemography\nToolkit"]
        LfTK["lf_toolkit\n\nlowFreq\nToolKit"]
        HfTK["hf_toolkit\n\nHighFreq\nToolKit"]
    end

    TestHera --> DataCfg
    TestHera --> HeraPrj
    TestHera --> ExpDir
    DataCfg --> ResultSet
    ResultSet --> ExpDir

    HeraPrj --> TopoTK
    HeraPrj --> LcTK
    HeraPrj --> DemoTK
    HeraPrj --> LfTK
    HeraPrj --> HfTK
```

---

## Phase 3: Test Execution

### Test Modules Overview

| Module | Toolkit | Tests | Key Patterns |
|--------|---------|-------|--------------|
| `test_datalayer.py` | (Project directly) | 5 | CRUD operations, counters, config |
| `test_repository.py` | `dataToolkit` | 7 | Repository add/get/load, path resolution |
| `test_topography.py` | `TopographyToolkit` | 13 | Point/list/grid elevation, STL, CRS conversion |
| `test_landcover.py` | `LandCoverToolkit` | 11 | Land cover at point/area, roughness, coding map |
| `test_lowfreq.py` | `lowFreqToolKit` | 18 | Analysis, presentation, data matching, edge cases |
| `test_highfreq.py` | `HighFreqToolKit` | 24 | Sonic/TRH data, calculators, turbulence statistics |
| `test_demography.py` | `DemographyToolkit` | 7 | Population calculations, area creation, defaults |

### Anatomy of a Toolkit Test

Here is the typical pattern, using `test_lowfreq.py` as an example:

```mermaid
flowchart TD
    Start["Test function receives\nlf_toolkit fixture\n(session-scoped)"] --> LoadData["Module fixture: lowfreq_df\nlf_toolkit.getDataSourceData(\n'YAVNEEL')"]

    LoadData --> Compute["Returns dask DataFrame\nCall .compute() to\nmaterialize to pandas"]

    Compute --> FixDates["Fix datetime column:\nParse dates, set index\nif needed"]

    FixDates --> CallMethod["Call analysis or\npresentation method:\nlf_toolkit.analysis\n.addDatesColumns(df, ...)"]

    CallMethod --> BasicAssert["Basic assertions:\nassert result is not None\nassert len(result) > 0\nassert expected columns exist"]

    BasicAssert --> NeedCompare{Compare with\nexpected output?}

    NeedCompare -- "Yes" --> CheckPrepare{PREPARE_EXPECTED_OUTPUT\nenv var set?}

    CheckPrepare -- "Yes (generation mode)" --> SaveOutput["save_expected_output(\nfilename, result,\noutput_type, expected_dir)"]
    SaveOutput --> PassTest["Test PASSES\n(output saved for\nfuture comparisons)"]

    CheckPrepare -- "No (comparison mode)" --> LoadExpected["load_expected_output(\nfilename, output_type,\nexpected_dir)"]
    LoadExpected --> Compare["compare_outputs(\nresult, expected,\noutput_type)\n\nDispatches to type-specific\ncomparison function"]
    Compare --> AssertMatch["assert comparison\nreturns True"]
    AssertMatch --> PassTest2["Test PASSES"]

    NeedCompare -- "No" --> DirectAssert["Direct assertions:\nassert value == expected\nassert shape == (n, m)"]
    DirectAssert --> PassTest3["Test PASSES"]
```

!!! note "Dask to Pandas"
    The `parquet` data handler returns a **dask DataFrame** for lazy loading. Test fixtures call `.compute()` to materialize it into a pandas DataFrame before running assertions.

### Test Data Mapping (test_repository.json)

| Toolkit Key | Config Entries | Datasources |
|-------------|----------------|-------------|
| `GIS_Raster_Topography` | `defaultSRTM: SRTMGL1` | `SRTMGL1` -> `measurements/GIS/raster` (string) |
| `GIS_LandCover` | `defaultLandCover: lc_mcd12q1` | `lc_mcd12q1` -> `measurements/GIS/raster/lc_mcd12q1.tif` (string) |
| `GIS_Demography` | — | `lamas_population` -> `measurements/GIS/vector/population_lamas.shp` (geopandas) |
| `MeteoLowFreq` | — | `YAVNEEL` -> `measurements/meteorology/lowfreqdata/YAVNEEL.parquet` (parquet) |
| `MeteoHighFreq` | — | `slicedYamim_sonic` + `slicedYamim_TRH` -> `measurements/meteorology/highfreqdata/` (parquet) |

---

## Comparison Helpers

The `conftest.py` module provides a rich set of comparison functions for validating test outputs against expected baselines.

### compare_outputs Dispatcher

```mermaid
flowchart TD
    Start["compare_outputs(\nresult, expected,\noutput_type)"] --> TypeSwitch{output_type\nvalue?}

    TypeSwitch -- "dataframe" --> CompareDF["compare_dataframes()\n\n1. Sort columns alphabetically\n2. Reset index\n3. Float cols: np.allclose\n   rtol=1e-6, atol=1e-6\n4. Datetime: strip timezone\n5. Other: direct equality"]

    TypeSwitch -- "geodataframe" --> CompareGDF["compare_dataframes()\n+ geometry comparison\n\n1. Sort by numeric columns\n2. Compare geometries via\n   symmetric_difference().area\n3. Area tolerance: GDF_TOL_AREA"]

    TypeSwitch -- "ndarray" --> CompareNP["np.allclose(\nresult, expected,\nrtol=1e-6, atol=1e-6)"]

    TypeSwitch -- "xarray" --> CompareXR["result.equals(expected)\n\nFull dataset comparison\nincluding coordinates"]

    TypeSwitch -- "dataarray" --> CompareDA["compare_dataarrays()\n\nnp.allclose on .values\n+ coordinate comparison"]

    TypeSwitch -- "float / int" --> CompareNum["math.isclose(\nresult, expected,\nrel_tol=1e-6)"]

    TypeSwitch -- "dict / list / tuple" --> CompareDeep["deep_compare_with_tolerance()\n\nRecursive comparison:\n- Floats: math.isclose\n- DataFrames: compare_dataframes\n- Arrays: np.allclose\n- Nested: recurse"]

    TypeSwitch -- "str / bytes" --> CompareStr["Direct equality\nresult == expected"]

    TypeSwitch -- "npz" --> CompareNPZ["For each array key:\nnp.allclose(\nresult, expected)"]

    CompareDF --> ReturnBool["Return True / False"]
    CompareGDF --> ReturnBool
    CompareNP --> ReturnBool
    CompareXR --> ReturnBool
    CompareDA --> ReturnBool
    CompareNum --> ReturnBool
    CompareDeep --> ReturnBool
    CompareStr --> ReturnBool
    CompareNPZ --> ReturnBool
```

### compare_dataframes — Deep Comparison

The `compare_dataframes` function handles several complex scenarios:

1. **Column Alignment** — Sorts columns alphabetically and resets index
2. **Numeric Tolerance** — Uses `np.allclose(rtol=1e-6, atol=1e-6)` for float columns
3. **Datetime Handling** — Strips timezone info before comparison
4. **Geometry Handling** — For GeoDataFrames, compares geometries via `symmetric_difference().area`
5. **Sort Stability** — For GeoDataFrames, sorts by preferred numeric columns to ensure deterministic comparison

### deep_compare_with_tolerance — Recursive Comparison

For nested structures (dicts, lists, tuples), the `deep_compare_with_tolerance` function recursively compares:

- **Floats** — `math.isclose(rel_tol=1e-6, abs_tol=1e-6)`
- **DataFrames** — Delegates to `compare_dataframes`
- **NumPy arrays** — `np.allclose`
- **Lists/Tuples** — Element-wise recursive comparison
- **Dicts** — Key-set comparison + recursive value comparison
- **Everything else** — Direct equality

---

## Expected Output Management

### Result Sets

Expected outputs are organized into **result sets** — named directories under `expected/`:

```
expected/
├── BASELINE/                        # Default result set
│   ├── getPointElevation.json
│   ├── expected_lowfreq_addDatesColumns.parquet
│   ├── expected_lowfreq_calcHourlyDist_density.npz
│   ├── create_xarray.nc
│   └── ...
└── REGRESSION_2025_11_11/           # Alternative result set
    └── ...
```

### Choosing a Result Set

The active result set is determined by (in priority order):

1. **CLI option:** `pytest --result-set REGRESSION_2025_11_11`
2. **Environment variable:** `RESULT_SET=REGRESSION_2025_11_11`
3. **Config default:** `data_config.json -> default_result_set`
4. **Hardcoded fallback:** `"BASELINE"`

### Save / Load Helpers

| Function | Purpose |
|----------|---------|
| `save_expected_output(filename, data, output_type, expected_dir)` | Serialize test output to the expected directory |
| `load_expected_output(filename, output_type, expected_dir)` | Deserialize expected output for comparison |

!!! warning "PREPARE_EXPECTED_OUTPUT Mode"
    Setting the `PREPARE_EXPECTED_OUTPUT=1` environment variable switches tests into **generation mode**: instead of comparing results against expected outputs, they **write** the current results as the new expected outputs. This is used when establishing a new baseline after intentional changes.

### Supported Output Formats

| output_type | Save Format | Load Method |
|-------------|-------------|-------------|
| `dataframe` | `.parquet` or `.json` | `pd.read_parquet()` / `pd.read_json()` |
| `geodataframe` | `.geojson` | `gpd.read_file()` |
| `xarray` / `dataarray` | `.nc` (NetCDF) | `xr.open_dataset()` / `xr.open_dataarray()` |
| `dict` / `list` | `.json` | `json.load()` |
| `float` / `int` | `.json` | `json.load()` + cast |
| `ndarray` / `npz` | `.npz` | `np.load()` |
| `str` / `string` | plain text | `open().read()` |

---

## Environment Variables

| Variable | Required | Default | Description |
|----------|----------|---------|-------------|
| `TEST_HERA` | No | `~/hera_unittest_data` | Path to the test data repository root |
| `RESULT_SET` | No | `BASELINE` | Name of the expected-output result set |
| `PREPARE_EXPECTED_OUTPUT` | No | (unset) | Set to `"1"` to generate expected outputs |
| `MPLBACKEND` | No | (system default) | Set to `Agg` for headless matplotlib |
| `GDF_TOL_AREA` | No | `1e-7` | Tolerance for geometry comparison area |

---

## Running Tests

### Quick Reference

```bash
# Activate the environment
cd /home/ilay/hera
source heraenv/bin/activate

# Run all tests
pytest hera/tests/ -v

# Run a specific module
pytest hera/tests/test_lowfreq.py -v

# Run a specific test class
pytest hera/tests/test_topography.py::TestGetPointElevation -v

# Run a single test function
pytest hera/tests/test_lowfreq.py::TestLowFreqToolkitInit::test_has_analysis -v

# Choose a result set
pytest hera/tests/ --result-set BASELINE -v

# Skip slow tests
pytest hera/tests/ -v -m "not slow"

# Generate expected outputs
PREPARE_EXPECTED_OUTPUT=1 pytest hera/tests/ -v
```

### Pytest Configuration (pytest.ini)

```ini
[pytest]
testpaths = hera/tests
python_files = test_*.py
python_classes = Test*
python_functions = test_*
addopts = -v --tb=short
markers =
    slow: marks tests as slow (deselect with '-m "not slow"')
    integration: marks tests that require MongoDB
```

---

## Adding New Tests

### Step-by-Step Guide

1. **Add test data** to `~/hera_unittest_data/measurements/<subdir>/`
2. **Update `test_repository.json`** with a new entry under the appropriate toolkit key
3. **Add a fixture** in `conftest.py` (session-scoped, depends on `hera_test_project`)
4. **Create a test module** `hera/tests/test_<name>.py`
5. **Use the fixture** to get a real toolkit instance — no file paths in tests
6. **Compare outputs** using `compare_outputs()` and expected files under `expected/BASELINE/`

### Example: Adding a New Toolkit Test

```python
# In conftest.py — add a session-scoped fixture
@pytest.fixture(scope="session")
def my_toolkit(hera_test_project):
    from hera.my_module import MyToolkit
    return MyToolkit(projectName=PYTEST_PROJECT_NAME)

# In test_my_toolkit.py
class TestMyToolkit:
    def test_basic(self, my_toolkit):
        data = my_toolkit.getDataSourceData("my_datasource")
        assert data is not None
        # ... assertions ...
```
