# Utilities

The `hera.utils` package provides a collection of modules for unit handling, data transformation, logging, and more. The public functions of the units, JSON, angle, logging and contour modules are re-exported from the top-level `hera.utils` namespace (loaded lazily on first access), so you can import them directly:

```python
from hera.utils import ureg, ConfigurationToJSON, toMeteorologicalAngle
```

The remaining modules (Slurm, statistics, DataFrame comparison, filtering, sensitivity analysis) are imported from their own submodule, as shown in the table.

## Available Modules

| Page | Module | Description |
|------|--------|-------------|
| [Unit Handling](units.md) | `hera.utils.unitHandler` (re-exported) | Shared pint registry `ureg` (with extra units such as `dunam` and `mmH2O`) and the `tonumber` / `tounit` conversion helpers. |
| [JSON Utilities](json.md) | `hera.utils.jsonutils` (re-exported) | Convert configurations with units to and from JSON, strip units, generate parameter sweeps (`JSONVariations`), compare and diff JSONs, set values by path. |
| [Angle Conversions](angles.md) | `hera.utils.angle` (re-exported) | Convert between mathematical, meteorological, and azimuth angle conventions. |
| [Logging Configuration](logging.md) | `hera.utils.logging` (re-exported) | Configure Hera's logging system, set per-logger levels, add file handlers and formatters, get loggers. |
| [Slurm Batch Jobs](slurm.md) | `hera.utils.slurm` | Generate a Slurm array-job submission script that runs a script for each directory in a list file. |
| [2D Distributions](statistics.md) | `hera.utils.statistics` | `calcDist2d`: normalized 2D histograms for contour-style visualizations. |
| [DataFrame Comparison](dataframes.md) | `hera.utils.dataframeutils` | Find the parameters whose values differ across named datasets, in long or wide format. |
| [Data Filtering](filtering.md) | `hera.utils.filter_immediate` | Chainable threshold and interval filters for pandas DataFrames. |
| [Sensitivity Analysis](sensitivity.md) | `hera.utils.SALibUtils` | Build SALib problem definitions with typed parameters and transform samples. Requires `pip install SALib` (not in `requirements.txt`). |
| [Contour to GIS](contours.md) | `hera.utils.matplotlibCountour` (`toGeopandas` re-exported) | Convert matplotlib contour output to GeoDataFrames for GIS export. Requires `unum`. |
