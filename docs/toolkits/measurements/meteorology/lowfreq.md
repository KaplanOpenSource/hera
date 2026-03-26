# MeteoLowFreq

**Toolkit name:** `MeteoLowFreq`

Hourly and daily meteorological station data with analysis and visualization.

```python
lf = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_LOWFREQ, projectName="MY_PROJECT")

# Load station data
df = lf.getDataSourceData("YAVNEEL")

# Analysis: enrich with date columns (year, month, season, etc.)
enriched = lf.analysis.addDatesColumns(df, datecolumn="datetime")

# Analysis: hourly distribution of a variable
hourly = lf.analysis.calcHourlyDist(enriched, field="wind_speed", density=True)

# Presentation: daily scatter plot
lf.presentation.dailyPlots.plotScatter(enriched, plotField="temperature")

# Presentation: seasonal hourly distribution
lf.presentation.seasonalPlots.plotSeasonalHourly(enriched, field="wind_speed")

# Presentation: probability contour by season
lf.presentation.seasonalPlots.plotProbContourf_bySeason(enriched, fields=["T", "RH"])
```

For the full API, see the [API Reference](../../developer_guide/api/measurements.md).
