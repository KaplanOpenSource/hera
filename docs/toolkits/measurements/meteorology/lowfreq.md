# MeteoLowFreq

**Toolkit constant:** `toolkitHome.METEOROLOGY_LOWFREQ` (`"MeteoLowFreq"`)
**Toolkit name in the database:** `lowFreqMeteorology`
**Implementation:** `hera.measurements.meteorology.lowfreqdata.toolkit.lowFreqToolKit`

Low-frequency meteorological station data (hourly/10-minute records from IMS and
TOA5 loggers). The toolkit enriches a time series with calendar and season
columns, computes hour-of-day distributions, and plots those distributions as
scatter, line and probability-contour figures — per day and per season.

```python
from hera import toolkitHome

lf = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_LOWFREQ, projectName="MY_PROJECT")

# Load station data
df = lf.getDataSourceData("YAVNEEL")

# Analysis: add year/month/day/time/season columns
enriched = lf.analysis.addDatesColumns(df)

# Presentation: scatter of a field against hour-of-day
lf.presentation.dailyPlots.plotScatter(enriched, plotField="TD")

# Presentation: probability contour of a field against hour-of-day
lf.presentation.dailyPlots.plotProbContourf(enriched, plotField="TD")

# Presentation: the same contour, one panel per season
lf.presentation.seasonalPlots.plotProbContourf_bySeason(enriched, plotField="TD")
```

For the full API, see the [API Reference](../../developer_guide/api/measurements.md).

---

## Data source format

| Property | Value |
|----------|-------|
| File format | Parquet (`dataFormat: "parquet"`) |
| Index | `DatetimeIndex` (time series) |
| Frequency | 10-minute to hourly |
| `docType` | `lowFreqMeteorology_LowFreqData` |

Column names come from the station data itself and are not fixed by the
toolkit — every analysis and plotting method takes the column to work on as an
argument. IMS parquet files, for example, carry `TD` (dry-bulb temperature),
`RH`, `WS` and `WD`.

---

## Initialising the toolkit

```python
from hera import toolkitHome

lf = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_LOWFREQ, projectName="MY_PROJECT")

# List available data sources
lf.getDataSourceList()

# Load one (dask or pandas, depending on how the source was registered)
df = lf.getDataSourceData("YAVNEEL")
```

`lowFreqToolKit` also accepts `filesDirectory` and `connectionName`, both
optional and both forwarded to `abstractToolkit`.

---

## Analysis — `lf.analysis`

### `addDatesColumns(data, datecolumn=None, monthcolumn=None)`

Adds calendar and season columns. Accepts a pandas DataFrame, a dask DataFrame,
or a **string path to a parquet file** (read with `pandas.read_parquet`).

| Argument | Meaning |
|----------|---------|
| `data` | DataFrame or path to a parquet file |
| `datecolumn` | Column holding the date. `None` (default) uses the index and stores it in a new `curdate` column |
| `monthcolumn` | Column holding the month. `None` (default) derives it into a new `monthonly` column |

Columns added:

| Column | Contents |
|--------|----------|
| `curdate` | The datetime (only when `datecolumn is None`) |
| `yearonly` | Year, integer |
| `monthonly` | Month 1–12 (only when `monthcolumn is None`) |
| `dayonly` | Day of month |
| `timeonly` | `datetime.time` |
| `Time` | Hour and minute as one integer, `hour*100 + minute` |
| `season` | `Winter` (Dec–Feb), `Spring` (Mar–May), `Summer` (Jun–Aug), `Autumn` (Sep–Nov) |

Timezone-aware inputs are converted to naive datetimes before the columns are
derived.

### `calcHourlyDist(data, Field, bins=30, normalization='density')`

Two-dimensional histogram of `Field` against hour-of-day. Note the capital **F**
in `Field`.

| Argument | Meaning |
|----------|---------|
| `data` | DataFrame with a datetime index, or a path to a parquet file |
| `Field` | Name of the column to build the distribution for |
| `bins` | Number of bins (default 30) |
| `normalization` | `'density'` (default), `'max_normalized'` or `'y_normalized'` — any other value raises `ValueError` |

Returns a 3-tuple `(x_mid, y_mid, M.T)`: the bin centres along the hour axis, the
bin centres along the value axis, and the transposed histogram. Rows where
`Field` is NaN or ≤ −5000 (the station no-data marker) are dropped first.

### `resampleSecondMoments(data, SamplingWindow, fieldsFirstMoments, fieldsSecondMoments)`

Resamples `fieldsFirstMoments` to `SamplingWindow` by mean and adds the
covariance of every pair drawn from `fieldsSecondMoments`. Expects the
`<field>_bar` columns the covariance is computed against to already be present.

---

## Presentation — `lf.presentation`

Two groups of plots hang off the presentation layer:

| Accessor | Class | Plots |
|----------|-------|-------|
| `lf.presentation.dailyPlots` | `DailyPlots` | `plotScatter`, `dateLinePlot`, `plotProbContourf` |
| `lf.presentation.seasonalPlots` | `SeasonalPlots` | `plotProbContourf_bySeason` |

All of them plot the chosen field against hour-of-day (0–24).

### `dailyPlots.plotScatter(data, plotField, ax=None, scatter_properties={}, ax_functions_properties={})`

Scatter of `plotField` against time of day, via `seaborn.scatterplot`.
`scatter_properties` updates the keyword arguments passed to seaborn;
`ax_functions_properties` adds or replaces the functions applied to the axes.

### `dailyPlots.dateLinePlot(data, plotField, date, legend=True, ax=None, line_properties={}, ax_functions_properties={})`

Line plot of `plotField` across a single day. `date` is a `'YYYY-MM-DD'` string.

### `dailyPlots.plotProbContourf(data, plotField, levels=None, scatter=True, withLabels=True, colorbar=True, Cmapname='jet', ax=None, ..., normalization='max_normalized')`

Filled probability contour of `plotField` against hour-of-day. Note that the
default normalization here is `'max_normalized'`, unlike `calcHourlyDist`.

| Argument | Meaning |
|----------|---------|
| `plotField` | The column to plot — **a single string**, not a list |
| `levels` | Overrides the default contour and contourf levels |
| `scatter` | Overlay the raw points (default `True`) |
| `withLabels` | Label the contour lines (default `True`) |
| `colorbar` | Draw a colorbar (default `True`) |
| `Cmapname` | Colormap name (default `'jet'`) |
| `scatter_properties`, `contour_values`, `contour_properties`, `contourf_properties`, `labels_properties`, `ax_functions_properties` | Dicts merged into the matplotlib/seaborn defaults |

### `seasonalPlots.plotProbContourf_bySeason(data, plotField, ..., figsize=[15, 10])`

Applies `plotProbContourf` once per season on a shared figure. Takes the same
arguments plus `figsize`, and requires the `season` column produced by
`addDatesColumns`.

---

## Complete example

```python
from hera import toolkitHome

lf = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_LOWFREQ, projectName="MY_PROJECT")

df = lf.getDataSourceData("YAVNEEL")
enriched = lf.analysis.addDatesColumns(df)

# Numbers
x_mid, y_mid, hist = lf.analysis.calcHourlyDist(enriched, Field="TD", bins=40)

# Figures
lf.presentation.dailyPlots.plotScatter(enriched, plotField="TD")
lf.presentation.dailyPlots.dateLinePlot(enriched, plotField="TD", date="2020-07-15")
lf.presentation.seasonalPlots.plotProbContourf_bySeason(enriched, plotField="TD")
```
