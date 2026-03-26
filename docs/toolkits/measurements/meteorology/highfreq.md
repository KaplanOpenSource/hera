# MeteoHighFreq

**Toolkit name:** `MeteoHighFreq`

High-frequency (10-20 Hz) sonic anemometer and TRH sensor data for turbulence analysis.

```python
hf = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_HIGHFREQ, projectName="MY_PROJECT")

# Load sonic anemometer data
sonic = hf.getDataSourceData("sonic_station_A")

# Analysis: calculate mean statistics
stats = hf.analysis.calculateMeanData(sonic)
```

For the full API, see the [API Reference](../../developer_guide/api/measurements.md).
