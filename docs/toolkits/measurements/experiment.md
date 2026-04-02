# Experiment Toolkit

**Toolkit name:** `experiment`

Manages experimental workflows — organizing raw data files into structured experiments with trials, device types, and sensors. Provides analysis (transmission frequency, turbulence, metadata enrichment) and presentation (device maps, heatmaps, LaTeX reports).

```python
from hera import toolkitHome

# Tip: if you created the project with `hera-project project create`, you can omit projectName
home = toolkitHome.getToolkit(toolkitHome.EXPERIMENT, projectName="MY_PROJECT")

# List available experiments
print(home.keys())  # ['IMS_experiment', 'Haifa2014']

# Get a specific experiment
exp = home.getExperiment("Haifa2014")
# or dictionary-style:
exp = home["Haifa2014"]

# Access trial data
trial = exp.trialSet["Measurements"]["Trial_01"]
df = trial.getData(deviceType="Sonic")

# Analyze transmission health
freq = exp.analysis.getDeviceTypeTransmissionFrequencyOfTrial(
    deviceType="Sonic", trialName="Trial_01"
)

# Visualize device functionality heatmap
exp.presentation.plotDeviceTypeFunctionality(
    deviceType="Sonic", trialName="Trial_01"
)
```

For the full API, see the [API Reference](../../developer_guide/api/measurements.md). For implementation details, see the [Developer Guide](../../developer_guide/measurements/experiment.md).

---

## Concepts

### The Argos data model

An experiment is defined in the **Argos** experiment management system ([ArgosWEB](https://github.com/KaplanOpenSource/argos)) and exported as a ZIP file. The data model has four core objects:

**Entity Types and Entities** — devices and sensors:

- An **Entity Type** is a class of device (e.g., "Sonic", "TRH", "Gateway"). It defines the attribute schema — which properties every device of this type has.
- An **Entity** is a specific device instance (e.g., "sonic01", "TRH_North"). It has its own attribute values.

**Trial Sets and Trials** — experimental configurations:

- A **Trial Set** groups related trials (e.g., "Measurements", "Calibration"). It defines the trial-level property schema.
- A **Trial** is a specific time-bounded experimental run. It assigns entities to locations, sets per-trial attribute values, and defines `TrialStart`/`TrialEnd` timestamps.

### Property scopes

Each attribute has a **scope** that determines where its value is set:

| Scope | Level | Changes per trial? | Example |
|-------|-------|-------------------|---------|
| **Constant** | Entity type | No — same for all devices of this type | `StoreDataPerDevice=false` |
| **Device** | Entity instance | No — fixed per device | `stationName="Check_Post"`, `height=9` |
| **Trial** | Per-device-per-trial | Yes — different in each trial | `location`, calibration values, thresholds |

### Containment hierarchy

Entities can be nested — a TRH sensor can be "contained in" a sonic anemometer station. Child entities **inherit** location and attributes from their parents. For example, if `TRH01` is contained in `sonic01`, and `TRH01` has no location set, it inherits `sonic01`'s location.

### Hera class hierarchy

In Hera, the Argos data model is extended with data-engine awareness:

| Level | Class | Description |
|-------|-------|-------------|
| **Experiment Home** | `experimentHome` | Factory — lists and retrieves experiments in a project |
| **Experiment** | `experimentSetupWithData` | A single experiment with its configuration, trials, and devices |
| **Trial Set** | `TrialSetWithData` | A named group of trials (e.g., "Measurements", "Calibration") |
| **Trial** | `TrialWithdata` | A single trial with start/end times and data access |
| **Entity Type** | `EntityTypeWithData` | A device type (e.g., "Sonic", "TRH") — all sensors of that kind |
| **Entity** | `EntityWithData` | A single sensor/device (e.g., "S01", "TRH_North") |

Each experiment has a **data engine** that handles the actual data retrieval — Parquet files, MongoDB via Pandas, or MongoDB via Dask. All trial and entity objects share the same engine instance.

---

---

## ArgosWEB — experiment authoring

[ArgosWEB](https://github.com/KaplanOpenSource/argos) is the web-based experiment management UI where experiments are designed before any code is written or data is collected. It is the **entry point** of the experiment pipeline.

### What you do in ArgosWEB

1. **Create an experiment** — give it a name, dates, and description
2. **Define entity types** — create device categories (e.g. "Sonic", "TRH", "Gateway") and define their attribute schema (which properties each device has, with types like String, Number, Boolean, Date)
3. **Create entity instances** — add individual devices (e.g. "sonic01", "sonic02") with their specific attribute values (station name, height, serial number)
4. **Set up trial sets and trials** — define experimental configurations with `TrialStart`/`TrialEnd` dates and per-trial properties
5. **Place devices on maps** — upload site images and position devices with lat/lon coordinates
6. **Configure containment** — nest devices (e.g. TRH sensor is "contained in" a sonic anemometer station — the TRH inherits the sonic's location)
7. **Export as ZIP** — download the experiment definition as a ZIP file containing `data.json` + map images

### How ArgosWEB connects to Hera

```
ArgosWEB (browser)
    │
    ▼ Export ZIP
experiment.zip
    │  └── data.json (trial sets, entity types, trials, entities, maps)
    │  └── images/ (site maps)
    │
    ▼ hera-experiment create
Hera experiment directory
    │  ├── code/           (experiment class)
    │  ├── data/           (parquet files — filled during data collection)
    │  ├── runtimeExperimentData/
    │  │   └── experiment.zip
    │  └── repository.json (for loading into Hera projects)
    │
    ▼ hera-project repository add + project create
Hera Project (MongoDB)
    │
    ▼ toolkitHome.getToolkit(EXPERIMENT)
Python analysis
```

### Alternative: live GraphQL access

Instead of exporting a ZIP, you can access ArgosWEB experiments directly via its GraphQL API:

```python
from argos.experimentSetup import webExperimentFactory

factory = webExperimentFactory(
    url="http://argos-server/graphql",
    token="your-auth-token",
)
experiment = factory.getExperiment("MyExperiment")
# Same API as file-based experiments — trial sets, entities, etc.
```

This is useful for automation or when you want to avoid the export step.

---

## Experiment lifecycle

### 1. Define in ArgosWEB

Create the experiment in the Argos web UI (see [ArgosWEB section above](#argosweb--experiment-authoring)):
- Define entity types and their attribute schemas
- Create entity instances (devices/sensors)
- Create trial sets and trials with `TrialStart`/`TrialEnd` dates
- Place devices on map images with coordinates
- Set up containment hierarchy (which sensor is on which station)
- Export as ZIP file

### 2. Create experiment directory

```bash
hera-experiment create MyExperiment --zip /path/to/exported.zip --path /experiments/
```

This creates the standard directory structure:

```
MyExperiment/
├── code/
│   └── MyExperiment.py              # Experiment class (customisable)
├── data/                            # Parquet files (one per device type)
│   ├── Sonic.parquet
│   └── TRH.parquet
├── runtimeExperimentData/
│   ├── Datasources_Configurations.json
│   └── MyExperiment.zip             # Argos metadata
└── MyExperiment_repository.json     # For loading into Hera projects
```

### 3. Collect data

During the experiment, data flows from sensors to Parquet files:

```
Devices → Node-RED → Kafka → pyArgos consumer → Parquet files
         (normalise)  (1 topic   (batch consume)   (data/ dir)
                      per type)
```

Or data can be loaded from Campbell binary/TOA5 files after the fact.

### 4. Load into Hera project

```bash
# Register repository (one-time)
hera-project repository add MyExperiment/MyExperiment_repository.json

# Create project (loads all registered repositories)
hera-project project create MY_PROJECT

# Or update existing project
hera-project project updateRepositories MY_PROJECT
```

### 5. Analyse in Python

```python
from hera import toolkitHome

home = toolkitHome.getToolkit(toolkitHome.EXPERIMENT, projectName="MY_PROJECT")
exp = home["MyExperiment"]

# Access trial data
df = exp.trialSet["Measurements"]["Trial_01"].getData(deviceType="Sonic")

# Analyse
exp.analysis.addTrialProperties(df, "Trial_01")
```

---

## Exploring experiment metadata

Before accessing data, you can inspect the experiment's structure:

```python
exp = home["MyExperiment"]

# Experiment configuration
print(exp.name)
print(exp.configuration)

# Entity types and their properties
for name, etype in exp.entityType.items():
    print(f"{name}: {etype.numberOfEntities} entities")
    print(etype.propertiesTable)       # attribute schema
    print(etype.entitiesTable)         # all devices as DataFrame

# Trial sets and trials
for ts_name, ts in exp.trialSet.items():
    print(f"Trial set: {ts_name}")
    print(ts.trialsTable)              # all trials as DataFrame
    for trial_name, trial in ts.items():
        print(f"  Trial: {trial_name}")
        print(f"    Start: {trial.properties['TrialStart']}")
        print(f"    End: {trial.properties['TrialEnd']}")
        print(trial.entitiesTable())   # devices in this trial with locations
```

---

### Data storage: `StoreDataPerDevice`

Each entity type (device type) has a `StoreDataPerDevice` flag that controls how measurement data is organized on disk:

| `StoreDataPerDevice` | Parquet file layout | Example |
|----------------------|--------------------|---------|
| `false` (default) | **One file per entity type** — all devices of that type in a single parquet file, with a `deviceName` column to distinguish them | `data/Sonic.parquet` contains data from sonic01, sonic02, ... |
| `true` | **One file per device** — each device has its own parquet file | `data/sonic01.parquet`, `data/sonic02.parquet`, ... |

This flag is defined in the experiment metadata (Argos zip file) as a `Constant`-scope property on the entity type. It affects:

- **How data is stored**: the repository JSON creates one `Experiment_rawData` document per type (if `false`) or per device (if `true`)
- **How data is queried**: when `StoreDataPerDevice=false`, the engine loads the single file and filters by `deviceName`; when `true`, it loads the specific device's file directly
- **CLI usage**: when using `hera-experiment data`, pass `--perDevice True` if the entity type stores data per device

```python
# StoreDataPerDevice=false (default): one file, filter by device name
df = trial.getData(deviceType="Sonic", deviceName="sonic01")
# Loads Sonic.parquet, filters to sonic01 rows

# StoreDataPerDevice=true: separate files per device
df = trial.getData(deviceType="PID", deviceName="PID_01")
# Loads PID_01.parquet directly
```

---

## Listing experiments

```python
home = toolkitHome.getToolkit(toolkitHome.EXPERIMENT, projectName="MY_PROJECT")

# List experiment names
home.keys()
# ['IMS_experiment', 'Haifa2014']

# Get a map of experiment names → datasource documents
home.getExperimentsMap()

# Get a formatted table of all experiments
home.getExperimentsTable()
```

---

## Loading an experiment

```python
exp = home.getExperiment("Haifa2014")

# Experiment properties
print(exp.name)                    # 'Haifa2014'
print(exp.configuration)           # full config dict
print(exp.defaultTrialSet)         # name of the default trial set

# Available trial sets and device types
print(list(exp.trialSet.keys()))   # ['Measurements', 'Calibration']
print(list(exp.entityType.keys())) # ['Sonic', 'TRH', 'PID']
```

---

## Accessing trial data

Trials are time-bounded segments of an experiment. Each trial has `TRIALSTART` and `TRIALEND` properties that are automatically used when you call `getData()` without specifying a time range.

```python
# Navigate: experiment → trial set → trial
trial = exp.trialSet["Measurements"]["Trial_01"]

# Get all Sonic data for this trial
df = trial.getData(deviceType="Sonic")

# Get data for a specific device
df = trial.getData(deviceType="Sonic", deviceName="S01")

# Get data with device metadata merged in
df = trial.getData(deviceType="Sonic", withMetadata=True)

# Override time range
df = trial.getData(
    deviceType="TRH",
    startTime="2024-03-15 08:00",
    endTime="2024-03-15 12:00"
)
```

### Shortcut: default trial set

```python
# Access trials from the default trial set directly
trials = exp.trialsOfDefaultTrialSet
```

---

## Accessing device data

Entity types (device types) and entities (individual devices) also provide data access:

```python
# All data for a device type
sonic_type = exp.entityType["Sonic"]
df_all = sonic_type.getData()

# Data for a device type during a specific trial
df_trial = sonic_type.getDataTrial(trialSetName="Measurements", trialName="Trial_01")

# Data for a single device
device = sonic_type["S01"]
df_device = device.getData()

# With time filtering
df_device = device.getData(startTime="2024-03-15 08:00", endTime="2024-03-15 12:00")
```

---

## Time-range queries

For queries not tied to a specific trial, use `getDataFromDateRange` on the experiment:

```python
df = exp.getDataFromDateRange(
    deviceType="TRH",
    startTime="2024-03-15 00:00",
    endTime="2024-03-16 00:00",
    deviceName="TRH_North",    # optional — all devices if omitted
    withMetadata=True           # merge device metadata
)
```

---

## Direct data engine access

For advanced use, access the data engine directly:

```python
engine = exp.getExperimentData()

# Parquet engine: lazy Dask DataFrame
dask_df = engine.getData(deviceType="Sonic", autoCompute=False)
pandas_df = dask_df.compute()

# Or compute immediately
pandas_df = engine.getData(deviceType="Sonic", autoCompute=True)

# Per-device organization
df = engine.getData(deviceType="Sonic", perDevice=True)
```

---

## Analysis

The analysis layer provides methods for device diagnostics, metadata enrichment, and turbulence calculations.

```python
analysis = exp.analysis
```

### Device locations

```python
locations = analysis.getDeviceLocations(
    entityTypeName="Sonic",
    trialName="Trial_01",
    trialSetName="Measurements"   # uses default trial set if omitted
)
# Returns DataFrame with device positions and metadata
```

### Transmission frequency

Analyze how reliably each device transmitted data during a trial:

```python
freq = analysis.getDeviceTypeTransmissionFrequencyOfTrial(
    deviceType="Sonic",
    trialName="Trial_01",
    trialSetName="Measurements",   # uses default trial set if omitted
    samplingWindow="1min",         # time bin size (default: "1min")
    normalize=True,                # normalize to planned message rate
    completeTimeSeries=True,       # fill gaps with zeros
    completeDevices=True,          # include non-transmitting devices
    wideFormat=True,               # pivot table (devices × time)
    recalculate=False              # use cached result if available
)
```

When `normalize=True`, values represent fraction of expected messages (1.0 = perfect). Results are cached in the data layer — set `recalculate=True` to force recomputation.

### Planned message count

```python
expected = analysis.getDeviceTypePlannedMessageCount(
    deviceType="Sonic",
    samplingWindow="1min"
)
# Returns float: expected messages per window
```

### Adding metadata to data

```python
# Merge device metadata (location, properties) into a DataFrame
df_with_meta = analysis.addMetadata(
    dataset=raw_df,
    trialName="Trial_01",
    trialSetName="Measurements"
)

# Add time-from-start and time-from-release columns
df_enriched = analysis.addTrialProperties(
    data=df_with_meta,
    trialName="Trial_01",
    trialSetName="Measurements"
)
# Adds columns: fromStart, fromRelease, fromStartSeconds, fromReleaseSeconds
```

### Turbulence statistics

For sonic anemometer data:

```python
stats = analysis.getTurbulenceStatistics(
    sonicData=sonic_df,
    samplingWindow="30min",
    height=10   # measurement height in meters
)
```

---

## Presentation

The presentation layer provides visualizations for experiment setup, device diagnostics, and reporting.

```python
pres = exp.presentation

# Control figure saving
pres.saveFigures = True
pres.savePath = "/path/to/output"
```

### Experiment site image

```python
# Plot an experiment site image with coordinate grid
ax = pres.plotImage(
    imageName="site_overview",
    withGrid=True,
    majorLocator=10   # grid spacing
)
```

### Device locations on map

```python
# Plot devices on an experiment map image
fig, ax = pres.plotDevicesOnImage(
    trialSetName="Measurements",
    trialName="Trial_01",
    deviceType="Sonic",
    mapName="floor_plan"
)

# Plot devices in ITM coordinates
fig, ax = pres.plotDevices(
    trialSetName="Measurements",
    trialName="Trial_01",
    deviceType="Sonic",
    mapName="site_overview"
)
```

### Device functionality heatmap

Visualize transmission health across devices and time. Color-codes each cell: red = no data, orange = partial, green = healthy.

```python
ax, pivot_table = pres.plotDeviceTypeFunctionality(
    deviceType="Sonic",
    trialName="Trial_01",
    trialSetName="Measurements",
    samplingWindow="1min",
    equalSquares=False   # True for square cells
)
```

### LaTeX report generation

Generate a PDF report with device maps and metadata tables:

```python
pres.generateLatexTable(
    latex_template="report_template.tex",   # Jinja2 template
    folder_path="/path/to/output"
)
```

---

## CLI reference

### List experiments

```bash
hera-experiment list --projectName MY_PROJECT
```

### Show experiment table

```bash
hera-experiment table --projectName MY_PROJECT
```

### Retrieve data

```bash
hera-experiment data Haifa2014 Sonic --projectName MY_PROJECT
hera-experiment data Haifa2014 TRH --deviceName TRH_North --perDevice True
```

### Create a new experiment

Scaffolds a complete experiment directory with boilerplate code, config files, and a repository JSON:

```bash
hera-experiment create my_experiment --path /path/to/experiments
hera-experiment create my_experiment --zip /path/to/argos_export.zip --relative
```

This creates:

```
my_experiment/
├── code/
│   └── my_experiment.py              # Experiment class (extends experimentSetupWithData)
├── data/                             # Place data files here
├── runtimeExperimentData/
│   ├── Datasources_Configurations.json
│   └── my_experiment.zip             # Argos metadata (if --zip provided)
└── my_experiment_repository.json     # Repository for loading into projects
```

The generated class provides hooks for custom analysis and presentation:

```python
class my_experiment(experimentSetupWithData):
    def __init__(self, projectName, pathToExperiment, filesDirectory):
        super().__init__(projectName, pathToExperiment, filesDirectory)
        self._analysis = my_experimentAnalysis(self)
        self._presentation = my_experimentPresentation(self, self.analysis)

class my_experimentAnalysis(experimentAnalysis):
    pass  # Add custom analysis methods here

class my_experimentPresentation(experimentPresentation):
    pass  # Add custom presentation methods here
```

### Load experiment into project

```bash
# Method 1: Register repository, then create/update project
hera-project repository add my_experiment/my_experiment_repository.json
hera-project project create MY_PROJECT
# or: hera-project project updateRepositories MY_PROJECT

# Method 2: Direct load
hera-experiment load --experiment /path/to/my_experiment MY_PROJECT
```

---

## Data engine types

The experiment toolkit supports three data backends, selected at initialization:

| Engine | Constant | Backend | Returns |
|--------|----------|---------|---------|
| Parquet (default) | `PARQUETHERA` | Hera data layer + Parquet files | `dask.DataFrame` (lazy) or `pandas.DataFrame` |
| Pandas/MongoDB | `PANDASDB` | Direct MongoDB queries | `pandas.DataFrame` |
| Dask/MongoDB | `DASKDB` | MongoDB via Dask | `dask.DataFrame` (lazy) |

To use a non-default engine when loading an experiment programmatically:

```python
from hera.measurements.experiment.dataEngine import PANDASDB

exp = experimentSetupWithData(
    projectName="MY_PROJECT",
    pathToExperiment="/path/to/experiment",
    dataType=PANDASDB
)
```

---

## Complete example

```python
from hera import toolkitHome

# Load experiment
# Tip: if you created the project with `hera-project project create`, you can omit projectName
home = toolkitHome.getToolkit(toolkitHome.EXPERIMENT, projectName="WindTunnel")
exp = home["march_2024"]

# Explore structure
print(f"Experiment: {exp.name}")
print(f"Trial sets: {list(exp.trialSet.keys())}")
print(f"Device types: {list(exp.entityType.keys())}")

# Get trial data with metadata
trial = exp.trialSet["Measurements"]["Release_01"]
df = trial.getData(deviceType="Sonic", withMetadata=True)

# Enrich with trial timing
df = exp.analysis.addTrialProperties(df, trialName="Release_01")
print(df[["deviceName", "fromReleaseSeconds", "wind_speed"]].head())

# Check device health
ax, freq = exp.presentation.plotDeviceTypeFunctionality(
    deviceType="Sonic",
    trialName="Release_01",
    samplingWindow="1min"
)

# Get device locations
locations = exp.analysis.getDeviceLocations(
    entityTypeName="Sonic",
    trialName="Release_01"
)

# Plot devices on site map
fig, ax = exp.presentation.plotDevicesOnImage(
    trialSetName="Measurements",
    trialName="Release_01",
    deviceType="Sonic",
    mapName="site_plan"
)
```
