# Troubleshooting

Common errors and solutions when working with the Hera system.

---

## Installation & Configuration

### `IOError: The config file doesn't exist in the default directory`

**Cause:** The MongoDB configuration file `~/.pyhera/config.json` does not exist or is not filled in.

**Solution:**

1. Hera auto-creates a template at `~/.pyhera/config.json` on first run (keyed by your OS user name, with `{username}`-style placeholders)
2. Edit the file with valid MongoDB connection details (or run `hera-project db create <name> --username ... --password ... --IP ... --databaseName ...`):

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

3. Replace `connectionName` with your OS user name — that connection is used when no `connectionName` is given. Other names are used only when passed explicitly (`Project(..., connectionName=...)`).

---

### `ConnectionRefusedError: [Errno 111] Connection refused` / `ServerSelectionTimeoutError`

**Cause:** MongoDB is not running or not accessible at the configured IP/port.

**Solution:**

If you use the Docker MongoDB set up by `init_with_mongo.sh` / `make install`:

```bash
make mongo-status   # is the hera-mongo container running?
make mongo-up       # start it
make mongo-logs     # inspect errors
```

For a system MongoDB service:

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
# Create the environment (venv at ~/.pyhera/environment, installs requirements + hera in editable mode)
cd /path/to/hera
make install-env

# Activate it (also sets HERA_REPO_ROOT and adds hera/bin to PATH)
source activate_hera.sh

# Verify
python -c "import hera; print(hera.__version__)"
```

---

### `command not found: hera-project`

**Cause:** The Hera environment is not activated, or `hera/bin/` is not in your system `PATH`.

**Solution:**

`source activate_hera.sh` adds `hera/bin` to `PATH`. Otherwise:

```bash
# Add to PATH (temporary)
export PATH="/path/to/hera/hera/bin:$PATH"

# Add to PATH (permanent - add to ~/.bashrc or ~/.zshrc)
echo 'export PATH="/path/to/hera/hera/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

---

## Toolkit Issues

### `ValueError: Toolkit '...' not found in static registry, ToolkitDataSource, or experiment toolkit`

**Cause:** The toolkit is not registered in the static registry or the database. (`getToolkit()` raises; it does not return `None`.)

**Solution:**

1. Check if the toolkit name is spelled correctly (case-sensitive):

```python
from hera import toolkitHome
# List all available toolkits
print(toolkitHome.getToolkitTable("MY_PROJECT"))
```

2. If the toolkit is custom/dynamic, register it:

```python
# The directory must contain MyCustomToolkit.py defining class MyCustomToolkit
toolkitHome.registerToolkit(
    toolkit_name="MyCustomToolkit",
    toolkit_path="/path/to/toolkits/MyCustomToolkit",
)
```

3. Or via CLI:

```bash
hera-project addToolkit MyCustomToolkit /path/to/toolkits/MyCustomToolkit
```

4. `ImportError: Cannot find toolkit class X.X` means the directory does not contain a module `X` with a class `X`; `ValueError: ... would shadow the already-importable module` means the toolkit name collides with an installed module — rename it.

---

### `getDataSourceData() returns None`

**Cause:** The datasource is not loaded into the project, or the name/version doesn't match.

**Solution:**

1. List available datasources:

```python
toolkit = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_LOWFREQ, projectName="MY_PROJECT")
print(toolkit.getDataSourceTable())
```

2. Check if the repository was loaded:

```bash
hera-project repository list
hera-project repository load /path/to/myRepo.json MY_PROJECT --overwrite
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

**Cause:** `PREPARE_EXPECTED_OUTPUT` is mentioned in the `hera/tests/conftest.py` docstring, but no test code reads it in master, so it never generates baseline files.

**Solution:**

There is no automatic baseline-generation mode. Write expected files explicitly with the `save_expected_output()` helper in `hera/tests/conftest.py`, and make sure the result-set directory exists:

```bash
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
pip install -r docs/requirements-docs.txt
# or: make docs-deps   (also installs hera with --no-deps)
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
