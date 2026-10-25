# MeteoHighFreq

**Toolkit constant:** `toolkitHome.METEOROLOGY_HIGHFREQ` (`"MeteoHighFreq"`)
**Toolkit name in the database:** `highFreqMeteorology`
**Implementation:** `hera.measurements.meteorology.highfreqdata.toolkit.HighFreqToolKit`

High-frequency (10–20 Hz) sonic anemometer and TRH sensor data for turbulence analysis. Provides a calculator-based analysis pipeline: parse raw binary/ASCII data, compute turbulence statistics, derive mean-field quantities (friction velocity, TKE, Monin-Obukhov length, stability classification, anisotropy), and chain calculations fluently.

```python
from hera import toolkitHome

# Tip: if you created the project with `hera-project project create`, you can omit projectName
hf = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_HIGHFREQ, projectName="MY_PROJECT")

# Create a turbulence calculator — pass a data source name or a DataFrame
turb = hf.analysis.singlePointTurbulenceStatistics(
    sonicData="sonic_station_A",   # loads from project automatically
    samplingWindow="30min",
    start="2024-03-15 00:00", end="2024-03-16 00:00",
    height=10, buildingHeight=5, averagedHeight=7,
)

# Compute fluctuations and second moments
turb.fluctuations().secondMoments()
result = turb.compute()

# Build mean-field statistics with method chaining
mdc = hf.analysis.MeanDataCalculator(TurbCalcOrData=turb)
mdc.horizontalSpeed().sigma().Ustar().TKE()
# MOLength()/StabilityMOLength() also need TRH mean temperature — see "Mean data calculator"
mean_data = mdc.compute()
```

For the full API, see the [API Reference](../../../developer_guide/api/measurements.md). For implementation details, see the [Developer Guide](../../../developer_guide/measurements/meteorology.md).

---

## Data source format

| Property | Value |
|----------|-------|
| File format | Parquet (`dataFormat: "parquet"`) |
| Expected columns | `u`, `v`, `w`, `T` (sonic); `TC_T`, `TRH`, `RH` (TRH sensors) |
| Index | `DatetimeIndex` (time-series) |
| Frequency | 10–20 Hz |

---

## Initialising the toolkit

```python
from hera import toolkitHome

hf = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_HIGHFREQ, projectName="MY_PROJECT")

# List available data sources
hf.getDataSourceList()
# ['sonic_station_A', 'TRH_station_A', ...]

# Load data (returns dask DataFrame by default)
sonic = hf.getDataSourceData("sonic_station_A")

# Convert to pandas if needed
df = sonic.compute()
```

---

## Loading raw data

The toolkit provides two ways to ingest raw sensor files (Campbell binary TOB1 or TOA5 ASCII). Both automatically detect the file format, rename the sonic columns `U`/`V`/`W` to `u`/`v`/`w` (and `TcT` to `TC_T`), set a `DatetimeIndex` named `Time`, convert values to float, and drop every other column (`RECORD`, `height`, `station`, `instrument` are moved into the metadata dict).

### `loadData` — parse, save, and register (recommended)

Parse a raw file, save the normalised output as Parquet, and register it as a versioned data source in the project. Once loaded, the data is available by name everywhere.

```python
# Parse + normalise + save as data source
doc = hf.loadData(
    name="sonic_10m",
    path="/raw_data/2024_03_15.dat",
    outputDirectory="/data/project/highfreq",   # where to store the parquet
    fromTime="2024-03-15 00:00",
    toTime="2024-03-16 00:00",
    parser="auto",          # auto-detect, or "campbell" / "toa5"
    version=(1, 0, 0),
    metadata={"station": "Haifa", "height": 10, "campaign": "March2024"},
)

# Now accessible by name everywhere:
turb = hf.analysis.singlePointTurbulenceStatistics(
    sonicData="sonic_10m",
    samplingWindow="30min",
    start="2024-03-15 08:00", end="2024-03-15 12:00",
    height=10, buildingHeight=5, averagedHeight=7,
)
```

#### Append vs overwrite

If a data source with the same name already exists:

```python
# Append new data to the existing parquet file
hf.loadData("sonic_10m", "/raw_data/2024_03_16.dat",
    outputDirectory="/data/project/highfreq",
    append=True)

# Or replace entirely
hf.loadData("sonic_10m", "/raw_data/2024_03_16.dat",
    outputDirectory="/data/project/highfreq",
    overwrite=True)
```

- **Neither flag** (default): raises ``ValueError`` if the data source exists
- **append=True**: loads the existing parquet, concatenates new data, deduplicates by timestamp
- **overwrite=True**: replaces the file and document
- Both ``True``: raises ``ValueError`` (mutually exclusive)

New data sources get a unique filename via the project counter `highfreq_<name>` (e.g. ``sonic_10m_0.parquet``). `loadData` returns one document, or a list of documents for multi-device files.

For multi-device files (TOA5 with multiple sonics), device names are appended automatically:

```python
docs = hf.loadData(name="station_A", path="/raw_data/multi_device.dat",
    outputDirectory="/data/project/highfreq")
# Creates: "station_A_Raw_Sonic_1", "station_A_Raw_Sonic_2", etc.
```

### `parseData` — parse and normalise without saving

For previewing data or one-off analysis without registering a data source. Returns parser-extracted metadata for inspection — useful for deciding what to store with `loadData`:

```python
results = hf.parseData(
    path="/raw_data/2024_03_15.dat",
    parser="auto",
)
# Returns list of (normalised_dataframe, parser_metadata_dict)

df, parser_meta = results[0]
print(df.head())
#                         u     v     w      T
# Time
# 2024-03-15 00:00:00  -2.65  3.05  1.96  24.93
# 2024-03-15 00:00:00  -2.70  3.10  1.90  24.95

print(parser_meta)
# {'station': 'CR3000_1', 'instrument': 'CSAT3', 'height': 10, 'deviceType': 'sonic'}
# deviceType is always set ('sonic' or 'TRH'); station/instrument/height only when the
# parser output had those columns. TOA5 results also carry 'deviceName' (e.g. 'Raw_Sonic_1').

# Pass directly to analysis:
turb = hf.analysis.singlePointTurbulenceStatistics(sonicData=df, ...)
```

### Supported formats

| Format | Header marker | Auto-detected | Parser name |
|--------|--------------|---------------|-------------|
| Campbell TOB1 binary | `TOB1` in first line | Yes | `"campbell"` |
| Campbell TOA5 ASCII | `TOA5` or CSV with device metadata | Yes | `"toa5"` |

### Normalised output

Both parsers normalise output to a consistent format:

| Device type | Columns | Index |
|-------------|---------|-------|
| Sonic | `u`, `v`, `w`, `T` (float) | `DatetimeIndex` |
| TRH | `TC_T`, `TRH`, `RH` (float) | `DatetimeIndex` |

---

## Analysis pipeline

The analysis layer follows a **calculator pattern**: create a calculator object, call methods to queue computations, then call `compute()` to execute them.

### Step 1: Create a turbulence calculator

You can pass either a **data source name** (string) or a **DataFrame** directly:

```python
# Option A: pass a data source name — the toolkit loads it automatically
turb = hf.analysis.singlePointTurbulenceStatistics(
    sonicData="sonic_station_A",  # data source name from the project
    samplingWindow="30min",
    start="2024-03-15 08:00",
    end="2024-03-15 12:00",
    height=10,
    buildingHeight=5,
    averagedHeight=7,
)

# Option B: pass a DataFrame directly
sonic = hf.getDataSourceData("sonic_station_A")
turb = hf.analysis.singlePointTurbulenceStatistics(
    sonicData=sonic,              # pandas or dask DataFrame with u, v, w, T columns
    samplingWindow="30min",       # averaging/resampling window
    start="2024-03-15 08:00",
    end="2024-03-15 12:00",
    height=10,                    # instrument height (m)
    buildingHeight=5,             # building height (m)
    averagedHeight=7,           # area-averaged building height (m)
    isMissingData=False,        # True if gaps exist in the data
)
```

`start` and `end` are stored as metadata only — the turbulence calculator does **not**
trim the raw data to them. They are applied later by `MeanDataCalculator`, which keeps
rows in `[start, end)`. Slice the DataFrame yourself if you want to limit the raw data.
`inmemory` is accepted but unused.

### Step 2: Compute fluctuations and second moments

```python
# Fluctuations: u', v', w', T', wind direction
turb.fluctuations()
# Result columns: u_bar, v_bar, w_bar, T_bar, wind_dir_bar
# (the raw data also gains up, vp, wp, Tp, wind_dir, wind_dir_p)

# Standard deviations
turb.sigma()
# Adds: sigmaU, sigmaV, sigmaW

# Horizontal sigma
turb.sigmaH()
# Adds: sigmaH = hypot(sigmaU, sigmaV) / sqrt(2)

# Horizontal wind speed
turb.horizontalSpeed()
# Adds: horizontal_speed_bar

# Wind direction standard deviation
turb.wind_dir_std()
# Adds: wind_dir_std

# Second-order moments (Reynolds stresses and fluxes)
turb.secondMoments()
# Adds: uu, vv, ww, uv, uw, vw, uT, vT, wT, TT

# Individual moments can also be called separately:
turb.uu().vv().ww().uw().vw().uv()
turb.uT().vT().wT().TT()

# Higher-order moments
turb.w3()   # w'^3 (skewness)
turb.w4()   # w'^4 (kurtosis)

# Friction velocity
turb.Ustar()
# Adds: Ustar = (uw² + vw²)^0.25

# Turbulent kinetic energy
turb.TKE()
# Adds: TKE = 0.5 * (uu + vv + ww)

# Monin-Obukhov length (from sonic temperature)
turb.MOLength_Sonic()
# Adds: L_Sonic

# Dimensionless ratios
turb.sigmaHOverUstar()
turb.sigmaWOverUstar()
turb.sigmaHOverWindSpeed()
turb.sigmaWOverWindSpeed()
turb.uStarOverWindSpeed()
turb.w3OverSigmaW3()

# Stability classification
turb.StabilityMOLength_Sonic()
# Adds: StabilityMOLength_Sonic (categorical: very unstable → very stable)

# Structure functions (dissipation rate estimation)
# Adds D11_<tau>s / D111_<tau>s and u_mag columns
turb.StrucFun(tau_range=[1, 2, 5, 10])
turb.ThirdStrucFun(tau_range=[1, 2, 5, 10])
```

Also available: `wTKE()`, `zoL_Sonic(zmd)`, `zOverL_Sonic()` (uses the effective height
`height + buildingHeight - 0.7 * averagedHeight`), `Lminus1_masked_Sonic()` and
`StrucFunDir(...)`.

### Step 3: Execute and retrieve

```python
# Execute all queued calculations
result = turb.compute()

# result is an InMemoryAvgData (DataFrame-like) with all computed columns
print(result.columns.tolist())
# ['u_bar', 'v_bar', 'w_bar', 'T_bar', 'sigmaU', 'sigmaV', 'sigmaW',
#  'uu', 'vv', 'ww', 'uw', 'vw', 'Ustar', 'TKE', ...]
```

---

## Mean data calculator

The `MeanDataCalculator` builds on the turbulence calculator output to derive higher-level mean-field statistics. All methods return `self` for **fluent chaining**.

On construction it calls `secondMoments().compute(mode=compute_mode_turb)` on the turbulence
calculator, joins the averaging-calculator output if given, and keeps only rows in
`[start, end)` from the calculator metadata.

!!! warning "Known issue in master"
    Passing a pandas/dask DataFrame as `TurbCalcOrData` or `AverageCalcOrData` raises
    `ValueError`: the type check compares against the `pandas`/`dask` *modules*
    (`meandatacalculator.py`, `__init__`). Pass calculator objects.

```python
# Create from a turbulence calculator
mdc = hf.analysis.MeanDataCalculator(TurbCalcOrData=turb)

# Or combine with an averaging calculator (for TRH mean fields)
avg_calc = hf.analysis.AveragingCalculator(
    deviceNameOrData=trh_df,
    samplingWindow="30min",
    start="2024-03-15 08:00", end="2024-03-15 12:00",
    height=10, buildingHeight=5, averagedHeight=7,
)
mdc = hf.analysis.MeanDataCalculator(
    TurbCalcOrData=turb,
    AverageCalcOrData=avg_calc,
)
```

### Available methods (all chainable)

**Time columns:**

```python
mdc.hour()           # Add "hour" column from index
mdc.timeWithinDay()  # Add fractional hour (hour + min/60 + sec/3600)
```

**Wind statistics:**

```python
mdc.horizontalSpeed()       # horizontal_speed_bar = hypot(u_bar, v_bar)
mdc.sigma()                 # sigmaU, sigmaV, sigmaW
mdc.sigmaAligned()          # sigmaU_aligned, sigmaV_aligned (rotated to wind direction)
mdc.sigmaH()                # horizontal std dev
mdc.Ustar()                 # friction velocity = (uw² + vw²)^0.25
```

**Dimensionless ratios:**

```python
mdc.sigmaHOverUstar()       # sigmaH / u*
mdc.sigmaUOverUstar()       # sigmaU_aligned / u*
mdc.sigmaVOverUstar()       # sigmaV_aligned / u*
mdc.sigmaWOverUstar()       # sigmaW / u*
mdc.sigmaHOverWindSpeed()   # turbulence intensity (horizontal)
mdc.sigmaWOverWindSpeed()   # turbulence intensity (vertical)
mdc.uStarOverWindSpeed()    # u* / wind speed
mdc.absWOverSigmaW()        # |w_bar| / sigmaW
```

**Reynolds stress alignment:**

```python
mdc.alignedStress()   # Rotate stress tensor to mean wind direction
# Adds: uu_aligned, uv_aligned, vv_aligned, uw_aligned, vw_aligned
```

**Energetics:**

```python
mdc.TKE()             # Turbulent kinetic energy = 0.5 * (uu + vv + ww)
mdc.Rvw()             # Correlation coefficient vw / sqrt(vv * ww)
mdc.Ruw()             # Correlation coefficient uw / sqrt(uu * ww)
```

**Stability:**

```python
mdc.MOLength()              # Monin-Obukhov length L — uses TC_T_bar, so it needs an
                            # AverageCalcOrData built from TRH data
mdc.effectivez()            # z_eff = height + buildingHeight - 0.7 * averagedHeight,
                            # stored in mdc.metaData['effectivez'] (not a column)
mdc.zOverL()                # Dimensionless stability z/L
mdc.StabilityMOLength()     # Classify: very unstable → very stable
```

**Anisotropy:**

```python
mdc.anisotropyEigs()   # Eigenvalues of anisotropy tensor + Lumley triangle coords
mdc.anisotropyCats()   # Classify: isotropic, 2-component axisymmetric, 1-component
```

**Dissipation rate:**

```python
# Needs the D11_<tau>s / D111_<tau>s and u_mag columns from turb.StrucFun / turb.ThirdStrucFun
mdc.StrucFun_eps(tau_range=[1, 2, 5, 10], rmin=0, rmax=10)       # adds eps_D11
mdc.ThirdStrucFun_eps(tau_range=[1, 2, 5, 10], rmin=0, rmax=10)  # adds eps_D111
```

**Filtering:**

```python
# Apply threshold filters (lt, lte, gt, gte, abs_lt, abs_lte, abs_gt, abs_gte, eq, neq)
mdc.thresholds([("horizontal_speed_bar", "gt", 1.0), ("sigmaW", "lt", 5.0)], inplace=True)

# REMOVE the rows inside [start, end) — filterDates keeps what is outside the interval
mdc.filterDates(start="2024-03-15 09:00", end="2024-03-15 11:00", inplace=True)
```

Use `inplace=True`: with the default `inplace=False` both methods try to build a new
`MeanDataCalculator` from a DataFrame, which hits the known issue above.

### Retrieve final result

```python
mean_data = mdc.compute()
# Returns pandas.DataFrame with all computed columns
```

---

## Averaging calculator

For computing simple time-averaged means (useful for TRH or supplementary variables):

```python
avg = hf.analysis.AveragingCalculator(
    deviceNameOrData=trh_df,
    samplingWindow="30min",
    start="2024-03-15 08:00", end="2024-03-15 12:00",
    height=10, buildingHeight=5, averagedHeight=7,
)

# Execute
result = avg.compute()
# Returns DataFrame with columns: TC_T_bar, TRH_bar, RH_bar
```

---

## Compute modes

The turbulence and averaging calculators accept a `mode` in `compute()`:

| Mode | Description |
|------|-------------|
| `'not_from_db_and_not_save'` | Compute locally, don't persist (default) |
| `'from_db_and_save'` | Check DB first; compute and save if missing |
| `'from_db_and_not_save'` | Check DB first; compute if missing, don't save |
| `'not_from_db_and_save'` | Compute locally and save to cache |

!!! warning "Only the default mode works in master"
    The DB modes reference `datalayer.Cache` on the `hera.datalayer.collection` module
    (`abstractcalculator.py`), which has no `Cache` attribute, so they raise
    `AttributeError`. The `*_save` modes also need `set_saveProperties(dataFormat=..., ...)`
    first, and the `parquet` saver calls dask's `repartition(partition_size=...)` on a pandas
    result. Use the default `'not_from_db_and_not_save'` and save the returned DataFrame
    yourself.

---

## Complete example

```python
from hera import toolkitHome
import pandas as pd

# Initialise
# Tip: if you created the project with `hera-project project create`, you can omit projectName
hf = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_HIGHFREQ, projectName="WindStudy")

# Load sonic data
sonic = hf.getDataSourceData("sonic_10m").compute()
trh = hf.getDataSourceData("trh_10m").compute()   # TC_T, TRH, RH — needed for MOLength

# Define analysis window
start = "2024-03-15 08:00"
end = "2024-03-15 12:00"

# Step 1: Turbulence statistics
turb = hf.analysis.singlePointTurbulenceStatistics(
    sonicData=sonic,
    samplingWindow="30min",
    start=start, end=end,
    height=10, buildingHeight=0, averagedHeight=0,
)
turb.fluctuations().secondMoments().sigma().Ustar().TKE()
turb_result = turb.compute()

avg = hf.analysis.AveragingCalculator(
    deviceNameOrData=trh,
    samplingWindow="30min",
    start=start, end=end,
    height=10, buildingHeight=0, averagedHeight=0,
)

# Step 2: Mean-field statistics with chaining
mdc = hf.analysis.MeanDataCalculator(TurbCalcOrData=turb, AverageCalcOrData=avg)
mdc.horizontalSpeed().hour().timeWithinDay()
mdc.sigmaHOverUstar().sigmaWOverUstar()
mdc.MOLength().StabilityMOLength()
mdc.anisotropyEigs().anisotropyCats()

mean_data = mdc.compute()

# Step 3: Inspect results
print(mean_data[["horizontal_speed_bar", "Ustar", "TKE", "L", "StabilityMOLength"]].head())
#                      horizontal_speed_bar  Ustar    TKE      L      StabilityMOLength
# 2024-03-15 08:00:00               4.2      0.31    1.45   -120.5   unstable
# 2024-03-15 08:30:00               3.8      0.28    1.20    -85.3   unstable
# 2024-03-15 09:00:00               5.1      0.42    2.10   -250.1   neutral/near neutral

# Step 4: Filter and analyse
n_total = len(mean_data)
mdc.thresholds([("StabilityMOLength", "eq", "stable")], inplace=True)
print(f"Stable periods: {len(mdc.compute())} out of {n_total}")
```
