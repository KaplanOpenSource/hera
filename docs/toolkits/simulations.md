# Simulations Toolkits

Simulation toolkits manage the lifecycle of computational models — setting up cases, running simulations, and post-processing results.

All simulation toolkits are accessed via `toolkitHome.getToolkit()` and bound to a project:

```python
from hera import toolkitHome

toolkit = toolkitHome.getToolkit(toolkitHome.SIMULATIONS_OPENFOAM, projectName="MY_PROJECT")
```

---

## OpenFOAM

**Toolkit name:** `OpenFOAM`

Full lifecycle management for OpenFOAM CFD simulations: templates, case setup, mesh generation, running, and post-processing.

### Templates

Templates are saved OpenFOAM case configurations that can be reused across simulations.

```python
of = toolkitHome.getToolkit(toolkitHome.SIMULATIONS_OPENFOAM, projectName="MY_PROJECT")

# Save a case directory as a template
of.saveTemplate(templateName="simpleFoam_base", caseDirectory="/path/to/case")

# Load a template to a new location
of.loadTemplate(
    templateName="simpleFoam_base",
    toDirectory="/data/simulations",
    caseName="run_001"
)
```

### Running simulations

```python
# Run an OpenFOAM case
of.runCase(casePath="/data/simulations/run_001", nproc=4)
```

### Post-processing

```python
# Access the analysis layer for data extraction
of.analysis.extractData(casePath="/data/simulations/run_001")
```

From the CLI:

```bash
# List templates
hera-openFoam simpleFoam templates list --projectName MY_PROJECT

# Save a template
hera-openFoam simpleFoam templates save myTemplate \
    --projectName MY_PROJECT --directory /path/to/case

# Load a template
hera-openFoam simpleFoam templates load myTemplate \
    --projectName MY_PROJECT --toDirectory /path/to/output --caseName run_001
```

---

## LSM (Lagrangian Stochastic Model)

**Toolkit name:** `LSM`

The LSM toolkit manages atmospheric dispersion simulations using Lagrangian particle tracking. It handles the full lifecycle: defining simulation templates, running simulations (including batch runs on Slurm clusters), and analyzing results (concentration fields, dosage calculations).

### Initializing the LSM toolkit

```python
lsm = toolkitHome.getToolkit(toolkitHome.LSM, projectName="MY_PROJECT")

# With options:
lsm = toolkitHome.getToolkit(
    toolkitHome.LSM,
    projectName="MY_PROJECT",
    to_xarray=True,      # save results as xarray (default: True)
    to_database=False,    # store simulation runs in the database (default: False)
    forceKeep=False       # keep original Lagrangian files after converting to xarray (default: False)
)
```

### Templates

LSM simulations are defined by **templates** — JSON descriptors that specify the model parameters (domain, meteorology, release conditions, etc.). Templates are stored as data sources:

```python
# Load a template from a JSON file
template = lsm.loadData(
    "/path/to/lsm_template.json",
    saveMode=lsm.TOOLKIT_SAVEMODE_FILEANDDB  # save to DB for reuse
)

# List all templates
templates = lsm.getTemplates()

# Get a template by name
template = lsm.getTemplateByName("urban_dispersion")

# View templates as a table
lsm.getTemplatesTable()
```

### Running simulations

Once you have a template, you can run simulations with specific parameters:

```python
# Run a simulation from a template
template.run(
    topography=topography_data,     # from GIS toolkit
    stations=weather_stations,      # meteorological stations DataFrame
    simulationName="run_001",
    windSpeed=5.0,
    windDirection=270,
    releaseRate=1.0
)
```

### Querying simulation results

```python
# Get simulations matching specific parameters
simulations = lsm.getSimulations(windSpeed=5.0, windDirection=270)

# List all simulations as a table
lsm.getSimulationsList()
```

### Working with a single simulation

Each simulation result is a `SingleSimulation` object that provides concentration and dosage calculations:

```python
sim = lsm.getSimulations(windSpeed=5.0)[0]

# Get concentration field (xarray Dataset with 'C' field)
# Q is the total released mass (default: 1 kg)
concentration = sim.getConcentration(Q=1*kg)

# The result is an xarray Dataset with:
# - 'C' field: concentration in mg/m³ at each (x, y, z, datetime)
# - attrs: dt (time step), Q (released mass), C (concentration units)

# Get dosage (cumulative exposure)
dosage = sim.getDosage(Q=1*kg)
# Returns xarray Dataset with 'Dosage' field

# Get concentration at a specific point and time
c_value = sim.getConcentrationAtPoint(x=100, y=200, datetime=target_time)
```

### Integration with Risk Assessment

LSM results feed directly into the risk assessment toolkit:

```python
# Get concentration from LSM
concentration = sim.getConcentration(Q=1e6*mg)

# Use in risk assessment
risk = toolkitHome.getToolkit(toolkitHome.RISKASSESSMENT, projectName="MY_PROJECT")
agent = risk.getAgent("Chlorine")
toxic_loads = agent["RegularPopulation"].calculateToxicLoads(concentration, field="C")
```

### Batch simulations with Slurm

For running many simulations with parameter variations on a cluster:

```python
lsm.prepareSlurmLSMExecution(
    baseParameters={"windSpeed": 5.0, "releaseRate": 1.0},
    jsonVariations={"windDirection": [0, 90, 180, 270]},
    templateName="urban_dispersion",
    stations=stations_df,
    topography="/path/to/topography",
    jobName="dispersion_sweep"
)
# Creates Slurm submission scripts for all parameter combinations
```

### CLI

```bash
# Load a template
hera-LSM load MY_PROJECT /path/to/template.json modelFolder "0,0,1" params.json

# List simulations
hera-LSM list MY_PROJECT
```

---

## Gaussian Dispersion

**Toolkit name:** `GaussianDispersion`

Gaussian puff and plume models for atmospheric pollutant transport.

```python
gauss = toolkitHome.getToolkit(toolkitHome.GAUSSIANDISPERSION, projectName="MY_PROJECT")

# The Gaussian toolkit provides methods for:
# - Cloud dispersion modeling
# - Meteorology integration (wind field, stability class)
# - Downwind concentration calculations
```

---

## Wind Profile

**Toolkit name:** `WindProfile`

Vertical wind profile modeling and analysis.

```python
wp = toolkitHome.getToolkit(toolkitHome.WINDPROFILE, projectName="MY_PROJECT")

# Wind profile analysis and modeling
```

---

## Hermes Workflows

**Toolkit name:** `hermesWorkflows`

Manages simulation groups and workflow pipelines that chain multiple simulation steps together.

```bash
# List workflow groups
hera-workflows list groups --projectName MY_PROJECT
```

For the full API details of each toolkit, see the [Toolkit Catalog](overview.md) and the [API Reference](../developer_guide/api/simulations.md).
