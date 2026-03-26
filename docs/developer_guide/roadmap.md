# Roadmap

Planned architectural improvements for future Hera releases.

---

## Contract-First Design with Typed Interfaces

**Status:** Planned

### Problem

Currently, Hera's internal interfaces rely on duck typing and dictionaries. Function signatures accept `str`, `dict`, or `geopandas.GeoDataFrame` interchangeably, and the `desc` metadata field is an untyped `dict`. This makes it hard to:

- Know what parameters a function actually expects without reading the source
- Validate inputs before they reach MongoDB
- Auto-generate API documentation with accurate type information
- Catch errors early (wrong column name, missing field, wrong CRS)

### Proposed changes

1. **Pydantic models for document metadata** — Define typed schemas for `desc` fields instead of free-form dicts:
   ```python
   class ToolkitDataSourceDesc(BaseModel):
       toolkit: str
       datasourceName: str
       version: tuple[int, int, int]
   ```

2. **Typed method signatures** — Replace `**kwargs` and `**desc` patterns with explicit typed parameters:
   ```python
   # Before
   def addDataSource(self, dataSourceName, resource, dataFormat, version=(0,0,1), **kwargs): ...

   # After
   def addDataSource(self, name: str, resource: Path, format: DataFormat, version: Version = (0,0,1), metadata: DataSourceMetadata | None = None): ...
   ```

3. **Enum for data formats** — Replace string constants with a proper enum:
   ```python
   class DataFormat(str, Enum):
       PARQUET = "parquet"
       NETCDF = "netcdf_xarray"
       GEOPANDAS = "geopandas"
       # ...
   ```

4. **Protocol classes for toolkit layers** — Define what `analysis` and `presentation` layers must implement:
   ```python
   class AnalysisProtocol(Protocol):
       @property
       def datalayer(self) -> abstractToolkit: ...

   class PresentationProtocol(Protocol):
       @property
       def datalayer(self) -> abstractToolkit: ...
   ```

### Migration path

- Phase 1: Add type annotations to all public methods (backward compatible)
- Phase 2: Introduce Pydantic models alongside existing dict interfaces
- Phase 3: Deprecate untyped interfaces with warnings
- Phase 4: Remove untyped interfaces in a major version

---

## Unified Toolkit Registry

**Status:** Planned

### Problem

Currently, toolkits are discovered from two different sources:

1. **Internal (built-in):** Hardcoded Python dict in `ToolkitHome.__init__` — requires source code edit to add/remove
2. **Dynamic (external):** Registered in MongoDB via `registerToolkit` CLI — no source code change needed

This creates:
- Two code paths for toolkit resolution (dict lookup vs DB query)
- Built-in toolkits can't be overridden without modifying source
- No single source of truth for "what toolkits are available"
- Adding a new built-in toolkit requires editing Python code

### Proposed changes

1. **Single registration mechanism** — All toolkits (built-in and external) are registered in the database using the same `ToolkitDataSource` document type

2. **Built-in registry JSON** — Ship a `toolkits_registry.json` with Hera:
   ```json
   {
       "GIS_Raster_Topography": {
           "classpath": "hera.measurements.GIS.raster.topography.TopographyToolkit",
           "type": "measurements"
       },
       "MeteoLowFreq": {
           "classpath": "hera.measurements.meteorology.lowfreqdata.toolkit.lowFreqToolKit",
           "type": "measurements"
       }
   }
   ```

3. **Registration command** — `make install` runs `hera-project registerBuiltins` which reads the JSON and registers all built-in toolkits in the DB

4. **Single resolution path** — `getToolkit(name)` always queries the DB. No hardcoded fallback dict.

5. **Populate command** — `make populate` loads all repositories into all projects (already implemented)

### Flow after unification

```
make install
    → mongo-up
    → hera-project registerBuiltins     ← reads toolkits_registry.json → DB
    → make populate                     ← loads repositories into all projects

Adding a new built-in toolkit:
    1. Create the class in hera/
    2. Add one line to toolkits_registry.json
    3. Run: make install

Adding an external toolkit:
    1. Create the class anywhere on disk
    2. Run: hera-project addToolkit myToolkit /path/to/toolkit
    3. Run: make populate (optional, to propagate to all projects)
```

### Migration path

- Phase 1: Create `toolkits_registry.json` and `registerBuiltins` command (keep hardcoded dict as fallback)
- Phase 2: `make install` runs registration automatically
- Phase 3: Remove hardcoded `_toolkits` dict — DB is the single source of truth
- Phase 4: Constants like `toolkitHome.GIS_RASTER_TOPOGRAPHY` remain as string aliases (no behavior change for users)

### Backward compatibility

- `toolkitHome.getToolkit("MeteoLowFreq")` continues to work — resolution just comes from DB instead of dict
- `toolkitHome.GIS_RASTER_TOPOGRAPHY` constant still works — it's just a string `"GIS_Raster_Topography"`
- Existing registered dynamic toolkits continue to work unchanged
- Users who never run `registerBuiltins` get the hardcoded fallback (during transition)

---

## Other Planned Improvements

### Environment variable configuration for MongoDB

Support `HERA_DB_*` environment variables as an override for `~/.pyhera/config.json`, enabling container and CI deployments without config files.

### Async support

Add `async` versions of database operations for use in web servers and Jupyter notebooks with event loops.

### Plugin system for data handlers

Allow third-party packages to register custom `DataHandler_*` classes without modifying `datahandler.py`.
