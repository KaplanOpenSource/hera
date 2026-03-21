# Best Practices Guide

Recommended approaches and patterns for working with Hera effectively.

---

## Project Organization

### Naming Conventions

**Project Names:**
- Use descriptive, uppercase names: `MY_PROJECT`, `WIND_ANALYSIS_2024`
- Avoid special characters and spaces
- Use underscores for separation: `RISK_ASSESSMENT_TEL_AVIV`

**Datasource Names:**
- Use clear, consistent names: `YAVNEEL`, `SRTMGL1`, `lamas_population`
- Match file names when possible (without extension)
- Use lowercase with underscores for multi-word names: `tel_aviv_station`

**Repository Names:**
- Descriptive and versioned: `meteo_data_v1`, `gis_base_repository`
- Include purpose: `test_data_repository`, `production_gis_data`

---

### When to Create Separate Projects

Create separate projects for:

- **Different domains** — Keep GIS, meteorology, and simulations in separate projects if they're independent
- **Different time periods** — `PROJECT_2023`, `PROJECT_2024` for temporal separation
- **Different experiments** — Each experiment gets its own project
- **Different users/teams** — Isolation and access control

Keep in one project when:

- **Related analyses** — All data feeds into the same analysis
- **Shared datasources** — Multiple toolkits use the same base data
- **Workflow dependencies** — Toolkits depend on each other's outputs

---

### Project Structure

Organize your project directory:

```
my_project/
├── caseConfiguration.json          # Auto-detected project name
├── data/                           # Local data files
│   ├── measurements/
│   ├── simulations/
│   └── cache/
├── scripts/                        # Analysis scripts
│   ├── analysis.py
│   └── visualization.py
└── results/                        # Output directory
    ├── plots/
    └── reports/
```

---

## Datasource Management

### Naming Conventions

- **Be descriptive** — `YAVNEEL` is better than `station1`
- **Include metadata in desc** — Don't rely on names alone
- **Use consistent patterns** — `STATION_NAME`, `DATASET_YEAR`, `SOURCE_VERSION`

### Versioning Strategy

**When to version:**

- **Data updates** — New processing, corrections, or additions
- **Format changes** — Schema changes in the data
- **Source changes** — Different data sources for the same concept

**Version numbering:**

- **Major** `[X, 0, 0]` — Breaking changes, incompatible formats
- **Minor** `[0, X, 0]` — New fields, backward-compatible additions
- **Patch** `[0, 0, X]` — Bug fixes, minor corrections

**Example:**

```python
# Initial version
toolkit.addDataSource("YAVNEEL", "v1.parquet", "parquet", version=[0, 0, 1])

# Added quality control flags (backward compatible)
toolkit.addDataSource("YAVNEEL", "v2.parquet", "parquet", version=[0, 1, 0])

# Major format change (incompatible)
toolkit.addDataSource("YAVNEEL", "v3.parquet", "parquet", version=[1, 0, 0])
```

### Setting Default Versions

Always set a default version for production use:

```python
# After loading repository
toolkit.setDataSourceDefaultVersion("YAVNEEL", [0, 1, 0])
```

This ensures consistent behavior across scripts and users.

---

## Repository Structure

### Organizing Repository JSON Files

**Single repository per domain:**

```
repositories/
├── gis_base.json              # GIS toolkits and datasources
├── meteo_stations.json        # Meteorological data
├── simulation_templates.json  # Simulation configurations
└── test_data.json            # Test/development data
```

**Or one repository per project:**

```
repositories/
├── project_alpha.json
├── project_beta.json
└── project_gamma.json
```

### Relative vs Absolute Paths

**Use relative paths when:**
- Repository and data are in the same directory tree
- Repository will be shared or moved
- Data is project-specific

**Use absolute paths when:**
- Data is on a shared network drive
- Multiple repositories reference the same data
- Data location is fixed and won't change

**Best practice:** Prefer relative paths for portability.

### Repository Naming

- Include purpose: `meteo_analysis_repo`, `gis_base_data`
- Include version if repositories evolve: `meteo_repo_v2`
- Use descriptive names: `test_data_repo` not `repo1`

---

## Performance

### Working with Large Datasets

**Use dask for lazy loading:**

```python
# Parquet files return dask DataFrames
df = toolkit.getDataSourceData("LARGE_DATASET")  # Not materialized yet

# Process in chunks
for partition in df.repartition(partition_size="100MB"):
    result = process_partition(partition.compute())
    save_result(result)
```

**Materialize only when needed:**

```python
# Bad: Materializes entire dataset
df = toolkit.getDataSourceData("LARGE").compute()  # Loads everything
result = df.head(1000)

# Good: Only load what you need
df = toolkit.getDataSourceData("LARGE")
result = df.head(1000).compute()  # Only loads first 1000 rows
```

### MongoDB Query Optimization

**Use specific filters:**

```python
# Bad: Loads all documents
docs = proj.getMeasurementsDocuments()

# Good: Filter by type and toolkit
docs = proj.getMeasurementsDocuments(
    type="ToolkitDataSource",
    toolkit="MeteoLowFreq"
)
```

**Leverage toolkit methods:**

```python
# Bad: Query all and filter in Python
all_docs = proj.getMeasurementsDocuments()
my_docs = [d for d in all_docs if d.desc.get("toolkit") == "MeteoLowFreq"]

# Good: Use toolkit's filtered method
toolkit = toolkitHome.getToolkit("MeteoLowFreq", projectName="MY_PROJECT")
docs = toolkit.getDataSourceDocuments()  # Already filtered
```

### Caching Strategies

**Cache computed results:**

```python
# Check cache first
cached = proj.getCacheDocuments(type="ProcessedStatistics", name="hourly_dist")
if cached:
    stats = cached[0].getData()
else:
    # Compute
    stats = compute_statistics(data)
    # Save to cache
    proj.addCacheDocument(
        resource=stats,
        dataFormat="JSON_dict",
        type="ProcessedStatistics",
        desc={"name": "hourly_dist"}
    )
```

**Use project config for small values:**

```python
# Store small configuration in project config
proj.setConfig(lastProcessedDate="2024-11-20", processingVersion="2.1")

# Retrieve
config = proj.getConfig()
last_date = config.get("lastProcessedDate")
```

---

## Data Formats

### Choosing the Right Format

| Format | Best For | Pros | Cons |
|--------|----------|------|------|
| `parquet` | Tabular data (pandas/dask) | Fast, compressed, columnar | Not human-readable |
| `netcdf_xarray` | Multi-dimensional arrays, time series | Standard for scientific data | Larger file size |
| `geopandas` | Vector spatial data | Preserves geometry, CRS | Requires geopandas |
| `JSON_dict` | Small configs, metadata | Human-readable, portable | Not efficient for large data |
| `CSV_PANDAS` | Simple tabular data | Human-readable, universal | Slower, larger files |
| `string` | Directory paths, simple values | Lightweight | Limited to strings |

**Recommendations:**

- **Tabular data** → `parquet` (best performance)
- **Spatial vector data** → `geopandas`
- **Time series / grids** → `netcdf_xarray`
- **Config files** → `JSON_dict`
- **Large datasets** → `parquet` with dask

### Format Conversion

**Convert when loading:**

```python
# Load as one format, convert to another
data = toolkit.getDataSourceData("DATA")  # Returns parquet DataFrame

# Convert to xarray if needed
import xarray as xr
ds = xr.Dataset.from_dataframe(data)
```

**Save in optimal format:**

```python
# Auto-detect and save
proj.saveData(
    name="processed_results",
    data=result_dataframe,
    desc={"processed": "2024-11-20"},
    kind="Measurements",
    type="ProcessedData"
)
# Automatically chooses parquet for DataFrame
```

---

## Toolkit Development

### Extending abstractToolkit

**Basic structure:**

```python
from hera import toolkit

class MyCustomToolkit(toolkit.abstractToolkit):
    def __init__(self, projectName, filesDirectory=None, connectionName=None):
        super().__init__(
            projectName=projectName,
            toolkitName="MyCustomToolkit",
            filesDirectory=filesDirectory,
            connectionName=connectionName
        )
        # Initialize analysis layer
        self._analysis = MyAnalysis(self)
        # Initialize presentation layer
        self._presentation = MyPresentation(self, self.analysis)
```

### Analysis Layer Pattern

```python
class MyAnalysis:
    def __init__(self, datalayer):
        self._datalayer = datalayer
    
    @property
    def datalayer(self):
        return self._datalayer
    
    def processData(self, data, **kwargs):
        # Processing logic
        return processed_data
```

### Presentation Layer Pattern

```python
class MyPresentation:
    def __init__(self, datalayer, analysis):
        self._datalayer = datalayer
        self._analysis = analysis
    
    def plotResults(self, data, **kwargs):
        import matplotlib.pyplot as plt
        # Plotting logic
        return ax
```

### Registering Custom Toolkits

**Option 1: Static registration** (built-in toolkits)

Edit `hera/toolkit.py`:

```python
self._toolkits["MyCustomToolkit"] = {
    "cls": "my_package.MyCustomToolkit",
    "desc": None,
    "type": "measurements"
}
```

**Option 2: Dynamic registration** (runtime)

```python
from hera import toolkitHome

toolkitHome.registerToolkit(
    toolkitclass="my_package.MyCustomToolkit",
    datasource_name="MyCustomToolkit",
    projectName="MY_PROJECT",
    repositoryName="defaultRepo",
    version=[0, 0, 1]
)
```

---

## Testing

### Writing Tests for Custom Toolkits

**Follow the existing test pattern:**

```python
# In conftest.py
@pytest.fixture(scope="session")
def my_toolkit(hera_test_project):
    from my_package import MyCustomToolkit
    return MyCustomToolkit(projectName=PYTEST_PROJECT_NAME)

# In test_my_toolkit.py
class TestMyCustomToolkit:
    def test_basic(self, my_toolkit):
        data = my_toolkit.getDataSourceData("my_datasource")
        assert data is not None
        assert len(data) > 0
    
    def test_analysis(self, my_toolkit):
        data = my_toolkit.getDataSourceData("my_datasource")
        result = my_toolkit.analysis.processData(data)
        assert result is not None
```

### Using the Test Infrastructure

- **Use session-scoped fixtures** — Share project across tests
- **Load test data via repository** — Use `test_repository.json` pattern
- **Compare with expected outputs** — Use `compare_outputs()` helper
- **Clean up in teardown** — Project cleanup is automatic

See [Testing Flow](../testing/flow.md) for complete details.

---

## Common Pitfalls

### Pitfall 1: Not Setting Default Versions

**Problem:**
```python
# Multiple versions exist
data = toolkit.getDataSourceData("YAVNEEL")  # Which version?
```

**Solution:**
```python
toolkit.setDataSourceDefaultVersion("YAVNEEL", [0, 1, 0])
data = toolkit.getDataSourceData("YAVNEEL")  # Always uses [0, 1, 0]
```

### Pitfall 2: Materializing Large Datasets Unnecessarily

**Problem:**
```python
df = toolkit.getDataSourceData("HUGE").compute()  # Loads everything
result = df.head(100)  # Only uses 100 rows
```

**Solution:**
```python
df = toolkit.getDataSourceData("HUGE")
result = df.head(100).compute()  # Only loads 100 rows
```

### Pitfall 3: Not Using Relative Paths in Repositories

**Problem:**
```json
{
    "resource": "/absolute/path/to/data.parquet"  // Not portable
}
```

**Solution:**
```json
{
    "isRelativePath": "True",
    "item": {
        "resource": "data/data.parquet"  // Portable
    }
}
```

### Pitfall 4: Mixing Collection Types

**Problem:**
```python
# Storing simulation results in Measurements
proj.addMeasurementsDocument(..., type="SimulationResult")
```

**Solution:**
```python
# Use the correct collection
proj.addSimulationsDocument(..., type="SimulationResult")
```

---

## See Also

- [Repository Examples](../examples/repository_examples.md) — Complete repository JSON examples
- [Workflow Examples](../examples/workflows.md) — End-to-end workflow patterns
- [Testing Flow](../testing/flow.md) — Testing best practices
- [Troubleshooting](../troubleshooting.md) — Common issues and solutions
