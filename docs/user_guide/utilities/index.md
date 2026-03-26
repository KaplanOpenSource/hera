# Utilities

The `hera.utils` package provides a collection of modules for unit handling, data transformation, logging, and more. Most public functions are re-exported from the top-level `hera.utils` namespace, so you can import directly:

```python
from hera.utils import ureg, ConfigurationToJSON, toMeteorologicalAngle
```

## Available Modules

| Module | Description |
|--------|-------------|
| [Unit Handling](units.md) | Physical-unit arithmetic and conversion using pint (`ureg`, `tonumber`, `tounit`). |
| [JSON Utilities](json.md) | Serialize configurations with units, generate parameter sweeps, compare JSON files. |
| [Angle Conversions](angles.md) | Convert between mathematical, meteorological, and azimuth angle conventions. |
| [Logging](logging.md) | Configure Hera's logging system, add file handlers, and set per-module levels. |
| [Slurm Batch Jobs](slurm.md) | Generate Slurm array-job submission scripts from a list of job directories. |
| [2D Statistics](statistics.md) | Compute normalized 2D histograms for contour-style visualizations. |
| [DataFrame Comparison](dataframes.md) | Compare parameter sets across datasets in long or wide format. |
| [Data Filtering](filtering.md) | Chainable threshold and interval filters for pandas DataFrames. |
| [Sensitivity Analysis](sensitivity.md) | Build SALib problem definitions with typed parameters and transform samples. |
| [Contour to GIS](contours.md) | Convert matplotlib contour output to GeoDataFrames for GIS export. |
