# Troubleshooting

Common errors and solutions when working with the Hera system.

---

## Installation & Configuration

### `IOError: The config file doesn't exist in the default directory`

**Cause:** The MongoDB configuration file `~/.pyhera/config.json` does not exist or is not filled in.

**Solution:**

1. Hera auto-creates a template at `~/.pyhera/config.json` on first run
2. Edit the file with valid MongoDB connection details:

```json
{
    "connectionName": {
        "username": "your_username",
        "password": "your_password",
        "dbIP": "localhost:27017",
        "dbName": "hera_db"
    }
}
```

3. Replace `connectionName` with a meaningful name (e.g., your OS username)

---

### `ConnectionRefusedError: [Errno 111] Connection refused`

**Cause:** MongoDB is not running or not accessible at the configured IP/port.

**Solution:**

```bash
# Check if MongoDB is running
sudo systemctl status mongod

# Start MongoDB
sudo systemctl start mongod

# Or start manually
mongod --dbpath /var/lib/mongodb --port 27017
```

Verify the `dbIP` in `~/.pyhera/config.json` matches your MongoDB server address.

---

### `ModuleNotFoundError: No module named 'hera'`

**Cause:** The Hera package is not installed in the current Python environment.

**Solution:**

```bash
# Activate the virtual environment
source heraenv/bin/activate

# Install in development mode
cd /path/to/hera
pip install -e .

# Verify
python -c "import hera; print(hera.__version__)"
```

---

### `command not found: hera-project`

**Cause:** The `hera/bin/` directory is not in your system `PATH`.

**Solution:**

```bash
# Add to PATH (temporary)
export PATH="/path/to/hera/hera/bin:$PATH"

# Add to PATH (permanent - add to ~/.bashrc or ~/.zshrc)
echo 'export PATH="/path/to/hera/hera/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

---

## Toolkit Issues

### `getToolkit() returns None`

**Cause:** The toolkit is not registered in the static registry or the database.

**Solution:**

1. Check if the toolkit name is spelled correctly (case-sensitive):

```python
from hera import toolkitHome
# List all available toolkits
print(toolkitHome.getToolkitTable("MY_PROJECT"))
```

2. If the toolkit is custom/dynamic, register it:

```python
toolkitHome.registerToolkit(
    toolkitclass="my_package.MyToolkit",
    datasource_name="MyCustomToolkit",
    projectName="MY_PROJECT"
)
```

3. Or via CLI:

```bash
hera-toolkit register --project MY_PROJECT --cls my_package.MyToolkit --name MyCustomToolkit
```

---

### `getDataSourceData() returns None`

**Cause:** The datasource is not loaded into the project, or the name/version doesn't match.

**Solution:**

1. List available datasources:

```python
toolkit = toolkitHome.getToolkit("MeteoLowFreq", projectName="MY_PROJECT")
print(toolkit.getDataSourceTable())
```

2. Check if the repository was loaded:

```bash
hera-project repository list
hera-project repository load myRepo MY_PROJECT --overwrite
```

3. Verify the datasource name matches exactly (case-sensitive).

---

### Version Conflicts in Datasources

**Cause:** Multiple versions of the same datasource exist and the default version is not set.

**Solution:**

```python
# List all versions
docs = toolkit.getDataSourceDocuments("YAVNEEL")
for doc in docs:
    print(f"Version: {doc.desc.get('version')}")

# Set the default version
toolkit.setDataSourceDefaultVersion("YAVNEEL", [0, 0, 2])
```

Or via CLI:

```bash
hera-project project version display MY_PROJECT --datasource YAVNEEL
hera-project project version update MY_PROJECT YAVNEEL "0,0,2"
```

---

## Testing Issues

### Tests Skip or Fail with Missing Data

**Cause:** The `TEST_HERA` environment variable is not set or points to a non-existent directory.

**Solution:**

```bash
# Check current setting
echo $TEST_HERA

# Set it
export TEST_HERA="$HOME/hera_unittest_data"

# Verify the directory contains required files
ls $TEST_HERA/test_repository.json
ls $TEST_HERA/data_config.json
ls $TEST_HERA/measurements/
```

---

### `PREPARE_EXPECTED_OUTPUT` Not Generating Files

**Cause:** The environment variable is not set correctly, or the expected directory doesn't exist.

**Solution:**

```bash
# Must be set to the string "1"
PREPARE_EXPECTED_OUTPUT=1 pytest hera/tests/ -v

# Verify expected directory exists
ls $TEST_HERA/expected/BASELINE/
```

---

### Geometry Comparison Fails with Small Differences

**Cause:** Floating-point precision differences in geometry operations.

**Solution:**

Increase the tolerance:

```bash
GDF_TOL_AREA=1e-5 pytest hera/tests/test_demography.py -v
```

Or set it in your environment permanently.

---

## Documentation Build Issues

### `griffe AliasResolutionError` During `mkdocs build`

**Cause:** The `griffe` static analyzer (used by `mkdocstrings`) fails to resolve wildcard imports like `from unum.units import *`.

**Solution:**

This has been fixed by replacing the wildcard import in `hera/utils/unitHandler.py` with an explicit module-level `globals().update()` pattern. If you encounter this error:

1. Ensure you have the latest version of `unitHandler.py`
2. Alternatively, remove the `:::` API reference blocks from the affected documentation pages

---

### `mkdocs serve` Fails with Plugin Errors

**Cause:** Missing documentation build dependencies.

**Solution:**

```bash
pip install -r requirements.txt
```

---

## Data Layer Issues

### `mongoengine.connection.ConnectionFailure`

**Cause:** Multiple incompatible MongoDB connections in the same process.

**Solution:**

Always use the same `connectionName` within a single script/session. If you need multiple connections, use distinct connection aliases.

---

### `TypeError: Object of type ObjectId is not JSON serializable`

**Cause:** Trying to serialize a MongoDB document directly to JSON.

**Solution:**

Use the `asDict()` method on documents:

```python
doc = proj.getMeasurementsDocuments(type="ToolkitDataSource")[0]
import json
json.dumps(doc.asDict(), default=str)
```

---

## Getting Help

If your issue is not covered here:

1. Check the [Testing Flow](testing/flow.md) page for test-specific documentation
2. Check the [Environment Variables](configuration/env_vars.md) page for configuration options
3. Review the source code — Hera uses Python `logging` extensively; enable DEBUG level for detailed output:

```python
import logging
logging.basicConfig(level=logging.DEBUG)
```
