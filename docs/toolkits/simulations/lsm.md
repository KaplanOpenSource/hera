# LSM (Lagrangian Stochastic Model)

**Toolkit name:** `LSM`

The LSM toolkit manages atmospheric dispersion simulations using Lagrangian particle tracking. It handles the full lifecycle: defining simulation templates, running simulations (including batch runs on Slurm clusters), and analyzing results (concentration fields, dosage calculations).

## Initializing the LSM toolkit

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

## Templates

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

## Running simulations

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

## Querying simulation results

```python
# Get simulations matching specific parameters
simulations = lsm.getSimulations(windSpeed=5.0, windDirection=270)

# List all simulations as a table
lsm.getSimulationsList()
```

## Working with a single simulation

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

## Integration with Risk Assessment

LSM results feed directly into the risk assessment toolkit:

```python
# Get concentration from LSM
concentration = sim.getConcentration(Q=1e6*mg)

# Use in risk assessment
risk = toolkitHome.getToolkit(toolkitHome.RISKASSESSMENT, projectName="MY_PROJECT")
agent = risk.getAgent("Chlorine")
toxic_loads = agent["RegularPopulation"].calculateToxicLoads(concentration, field="C")
```

## Batch simulations with Slurm

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

## CLI

```bash
# Load a template
hera-LSM load MY_PROJECT /path/to/template.json modelFolder "0,0,1" params.json

# List simulations
hera-LSM list MY_PROJECT
```

For the full API, see the [API Reference](../../developer_guide/api/simulations.md).
