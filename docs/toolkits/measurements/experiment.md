# Experiment Toolkit

**Toolkit name:** `experiment`

Manages experimental workflows — organizing raw data files into structured experiments with metadata tracking.

```python
exp = toolkitHome.getToolkit(toolkitHome.EXPERIMENT, projectName="MY_PROJECT")

# List available experiments
exp.getExperimentsMap()

# Get a specific experiment
my_exp = exp.getExperiment(experimentName="wind_tunnel_march_2024")
```

For the full API, see the [API Reference](../../developer_guide/api/measurements.md).
