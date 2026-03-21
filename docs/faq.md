# Frequently Asked Questions (FAQ)

Common questions and answers about using Hera.

---

## Project Management

### How do I share a project with someone else?

Use the export/import functionality:

```python
# On the source machine
from hera import Project
proj = Project(projectName="MY_PROJECT")
proj.export("/path/to/project_backup.zip")

# Transfer the zip file to the other machine

# On the destination machine
from hera import Project
proj = Project(projectName="MY_PROJECT")
Project.load(proj, "/path/to/project_backup.zip", is_hard_import=False)
```

Or use the CLI:

```bash
# Export
hera-project project dump MY_PROJECT --format json --fileName backup.json

# Import
hera-project project load OTHER_PROJECT backup.json
```

!!! note "MongoDB Connection"
    The destination machine must have MongoDB running and configured in `~/.pyhera/config.json`. The data files referenced by documents must also be accessible (same paths or update paths after import).

---

### How do I move a project to another computer?

1. **Export the project** (see above)
2. **Copy the data files** — All files referenced in `resource` fields must be copied to the new machine
3. **Update paths if needed** — If file paths changed, update the `resource` fields in MongoDB documents
4. **Import the project** on the new machine
5. **Verify MongoDB connection** — Ensure `~/.pyhera/config.json` is configured

---

### What's the difference between creating a project via CLI vs Python?

Both methods create the same project structure:

- **CLI**: `hera-project project create MY_PROJECT` — Also loads repositories if `--noRepositories` is not used
- **Python**: `Project(projectName="MY_PROJECT")` — Creates project on first access, no repository loading

The CLI method is more convenient for initial setup with repositories.

---

## Data Organization

### What's the difference between Measurements and Simulations collections?

| Collection | Purpose | Typical Content |
|------------|---------|-----------------|
| **Measurements** | Observational/input data | Station data, GIS files, toolkit datasources, raw experimental data |
| **Simulations** | Model outputs | OpenFOAM results, LSM concentration fields, simulation outputs |
| **Cache** | Intermediate/computed data | Processed results, configurations, temporary calculations |

The distinction is semantic — you can store any document type in any collection, but following conventions makes data easier to find and manage.

---

### When should I use Cache vs Measurements?

Use **Cache** for:
- Project configuration (`type = "<projectName>__config__"`)
- Intermediate processing results
- Temporary computed values
- Data that's derived from other documents

Use **Measurements** for:
- Raw observational data
- Toolkit datasources (`type = "ToolkitDataSource"`)
- Input data files
- Primary data sources

---

### How do I query documents across all collections?

Query each collection separately:

```python
proj = Project(projectName="MY_PROJECT")

# Query measurements
meas = proj.getMeasurementsDocuments(type="ToolkitDataSource")

# Query simulations
sims = proj.getSimulationsDocuments(type="MySimulationType")

# Query cache
cache = proj.getCacheDocuments(type="MyCacheType")
```

There's no cross-collection query — the `_cls` discriminator enforces separation.

---

## Toolkits

### How do I add a custom toolkit?

**Option 1: Register dynamically (runtime)**

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

**Option 2: Add to static registry**

Edit `hera/toolkit.py` and add to the `_toolkits` dict in `ToolkitHome.__init__()`:

```python
self._toolkits["MyCustomToolkit"] = {
    "cls": "my_package.MyCustomToolkit",
    "desc": None,
    "type": "measurements"
}
```

**Option 3: Via repository JSON**

Include a `Registry` section with classpath hints (see [Repository Schema](reference/repository_schema.md)).

---

### Why does `getToolkit()` return None?

Possible causes:

1. **Toolkit name misspelled** — Names are case-sensitive
2. **Toolkit not registered** — Check `toolkitHome.getToolkitTable(projectName)`
3. **Custom toolkit not loaded** — Use `registerToolkit()` first
4. **Import error** — The toolkit class can't be imported (check Python path)

**Debug steps:**

```python
from hera import toolkitHome

# List all available toolkits
table = toolkitHome.getToolkitTable("MY_PROJECT")
print(table)

# Check if toolkit exists
if "MeteoLowFreq" in table.index:
    print("Toolkit found")
else:
    print("Toolkit not registered")
```

---

### Can I use multiple toolkits in the same project?

Yes. Each toolkit operates independently within the same project:

```python
from hera import toolkitHome

# Get multiple toolkits
topo = toolkitHome.getToolkit(toolkitHome.GIS_RASTER_TOPOGRAPHY, projectName="MY_PROJECT")
meteo = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_LOWFREQ, projectName="MY_PROJECT")
risk = toolkitHome.getToolkit(toolkitHome.RISKASSESSMENT, projectName="MY_PROJECT")

# All share the same project's data layer
# Each toolkit only sees its own datasources (filtered by toolkit name)
```

---

## Repository

### How do I create a repository JSON?

1. **Start with a minimal structure:**

```json
{
    "MeteoLowFreq": {
        "DataSource": {
            "YAVNEEL": {
                "isRelativePath": "True",
                "item": {
                    "resource": "data/YAVNEEL.parquet",
                    "dataFormat": "parquet",
                    "version": [0, 0, 1],
                    "desc": {
                        "stationName": "YAVNEEL"
                    }
                }
            }
        }
    }
}
```

2. **Save as `repository.json`**
3. **Register it:**

```bash
hera-project repository add myRepo /path/to/repository.json
```

4. **Load into a project:**

```bash
hera-project repository load myRepo MY_PROJECT
```

See [Repository Examples](examples/repository_examples.md) for complete examples.

---

### What's the difference between DataSource and Measurements in repository JSON?

| Section | Purpose | Document Type | Typical Use |
|---------|---------|---------------|-------------|
| **DataSource** | Toolkit-managed datasources | `Measurements` with `type="ToolkitDataSource"` | Data that toolkits access via `getDataSourceData()` |
| **Measurements** | Raw measurement documents | `Measurements` with custom `type` | General measurement data, not toolkit-specific |
| **Simulations** | Simulation outputs | `Simulations` | Model results |
| **Cache** | Cached/computed data | `Cache` | Intermediate results |

**DataSource** entries are versioned and managed by the toolkit. **Measurements** entries are general-purpose documents.

---

### How does path resolution work in repository JSON?

The `isRelativePath` flag controls path resolution:

- **`"True"` or `true`** — Path is relative to the repository JSON file's directory
- **`"False"` or `false`** — Path is absolute (used as-is)

```json
{
    "MeteoLowFreq": {
        "DataSource": {
            "YAVNEEL": {
                "isRelativePath": "True",  // Relative to repo.json location
                "item": {
                    "resource": "data/YAVNEEL.parquet"  // Resolved to: /path/to/repo/data/YAVNEEL.parquet
                }
            },
            "OTHER": {
                "isRelativePath": "False",  // Absolute path
                "item": {
                    "resource": "/absolute/path/to/data.parquet"  // Used as-is
                }
            }
        }
    }
}
```

---

## Versioning

### How do I manage multiple versions of the same datasource?

1. **Add multiple versions:**

```python
toolkit.addDataSource("YAVNEEL", "/data/v1.parquet", "parquet", version=[0, 0, 1])
toolkit.addDataSource("YAVNEEL", "/data/v2.parquet", "parquet", version=[0, 0, 2])
toolkit.addDataSource("YAVNEEL", "/data/v3.parquet", "parquet", version=[0, 0, 3])
```

2. **Set default version:**

```python
toolkit.setDataSourceDefaultVersion("YAVNEEL", [0, 0, 2])
```

3. **Access specific version:**

```python
# Uses default version
data = toolkit.getDataSourceData("YAVNEEL")

# Explicit version
data = toolkit.getDataSourceData("YAVNEEL", version=[0, 0, 1])
```

4. **List all versions:**

```python
docs = toolkit.getDataSourceDocuments("YAVNEEL")
for doc in docs:
    print(f"Version: {doc.desc.get('version')}")
```

---

### What happens if I don't specify a version?

The system uses this priority:

1. **Default version** (if set via `setDataSourceDefaultVersion()`)
2. **Highest version** (if multiple versions exist, picks the one with the highest tuple)
3. **Single version** (if only one version exists, uses it)
4. **None** (if no datasource found)

---

## Performance

### How do I work with large datasets?

Use **dask** for lazy loading:

```python
# Parquet files return dask DataFrames
df = toolkit.getDataSourceData("LARGE_DATASET")  # Returns dask DataFrame

# Materialize only when needed
df_small = df.head(1000).compute()  # Only loads first 1000 rows

# Process in chunks
for chunk in df.repartition(partition_size="100MB"):
    result = process_chunk(chunk.compute())
```

See [Best Practices: Performance](guides/best_practices.md#performance) for more tips.

---

### How do I optimize MongoDB queries?

- **Use specific filters** — Narrow down queries with `type`, `toolkit`, etc.
- **Index frequently queried fields** — MongoDB automatically indexes `projectName` and `_cls`
- **Limit result sets** — Use `getDataSourceList()` instead of loading all documents
- **Cache config** — Project config is cached after first access

---

## Troubleshooting

### Common Issues

For detailed troubleshooting, see the [Troubleshooting Guide](troubleshooting.md).

**Quick fixes:**

- **`IOError: config file doesn't exist`** → Create `~/.pyhera/config.json`
- **`ConnectionRefusedError`** → Start MongoDB
- **`getToolkit() returns None`** → Check toolkit name spelling, verify registration
- **`getDataSourceData() returns None`** → Verify datasource name, check if repository was loaded

---

## Getting More Help

- **Documentation**: Browse the [full documentation index](index.md)
- **Examples**: See [Repository Examples](examples/repository_examples.md) and [Workflows](examples/workflows.md)
- **Best Practices**: Read the [Best Practices Guide](guides/best_practices.md)
- **Troubleshooting**: Check the [Troubleshooting Guide](troubleshooting.md)
