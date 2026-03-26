# Experiment toolkit

The experiment system organizes raw data files into structured experiments:

| Component | Class | Purpose |
|-----------|-------|---------|
| Toolkit | `experimentHome` | Create/list/get experiments |
| Analysis | `experimentAnalysis` | Experiment-specific analysis |
| Data engine | `parquetDataEngineHera` | Parquet-based data storage |
| Presentation | Experiment presentation | Visualization layer |

Experiments are stored as measurement documents and can be exposed as toolkits through `ToolkitHome.getExperimentToolkitDocuments()`.
