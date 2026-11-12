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
# Load a template from a JSON file (the JSON must contain a "name" key)
template = lsm.loadData(
    "/path/to/lsm_template.json",
    saveMode=lsm.TOOLKIT_SAVEMODE_FILEANDDB,  # save to DB for reuse (ValueError if it exists;
                                              # use TOOLKIT_SAVEMODE_FILEANDDB_REPLACE to replace)
    version=[0, 0, 1],                        # optional; other keyword arguments are stored in the description
)

# List all templates
templates = lsm.getTemplates()

# Get a template by name (None if not found)
template = lsm.getTemplateByName("urban_dispersion", templateVersion=None)

# View templates as a table
lsm.getTemplatesTable()
```

## Running simulations

Once you have a template, you can run simulations with specific parameters:

```python
# Run a simulation from a template
sim = template.run(
    topography=topography_data,     # topography file content (string); None = flat, homogeneous wind
    stations=weather_stations,      # meteorological stations DataFrame (optional)
    canopy=None,                    # optional canopy file
    depositionRates=None,           # optional deposition velocity (or list)
    simulationName="run_001",       # default: LSM_Simulation_<counter>
    saveMode=lsm.TOOLKIT_SAVEMODE_FILEANDDB,
    params={},                      # parameters overriding the template defaults
    windSpeed=5.0,                  # extra keyword arguments are merged into params
    windDirection=270,
    releaseRate=1.0
)
# Returns a SingleSimulation when the toolkit was created with to_xarray=True, otherwise None
```

The model runs in `<filesDirectory>/<simulationName>`. If the template defines `units` for a parameter, pass that parameter as a pint quantity (e.g. `windSpeed=5*ureg.m/ureg.s`) or a unit string (`"5 m/s"`).

## Querying simulation results

```python
from hera.utils import ureg

# Get simulations (list of SingleSimulation) matching specific parameters
simulations = lsm.getSimulations(windSpeed=5*ureg.m/ureg.s, windDirection=270*ureg.degree)
simulations = lsm.getSimulations(simulationName="run_001")

# List all simulations as a table (long format; wideFormat=True pivots it)
lsm.getSimulationsList(wideFormat=True)
```

`getSimulations` converts the query values with the `units` of the template named by `unitsTemplateVersion` (default `"v4-general"`), so that template must be loaded in the project; when it defines units, query values must be pint (or unum) quantities, otherwise a `ValueError` is raised.

## Working with a single simulation

Each simulation result is a `SingleSimulation` object that provides concentration and dosage calculations:

```python
sim = lsm.getSimulations(simulationName="run_001")[0]

# Get concentration field (xarray Dataset with 'C' field)
# Q is the total released mass (default: 1 kg)
concentration = sim.getConcentration(Q=1*ureg.kg)

# The result is an xarray Dataset with:
# - 'C' field: concentration in mg/m³ at each (x, y, z, datetime)
# - attrs: dt (time step), Q (released mass), C (concentration units)

# Get dosage (cumulative exposure)
dosage = sim.getDosage(Q=1*ureg.kg)
# Returns xarray Dataset with 'Dosage' field

# Get concentration at a specific point and time
c_value = sim.getConcentrationAtPoint(x=100, y=200, datetime=target_time)

# All three accept Q, time_units (default min) and q_units (default mg)
```

## Integration with Risk Assessment

LSM results feed directly into the risk assessment toolkit:

```python
# Get concentration from LSM
concentration = sim.getConcentration(Q=1e6*ureg.mg)

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
    jsonVariations=[{"windDirection": [0, 90, 180, 270]}],   # list of variation groups (or a path to a JSON file)
    templateName="urban_dispersion",
    stations=stations_df,
    topography="/path/to/topography",
    jobName="dispersion_sweep"
)
# Writes <filesDirectory>/simulationsScripts/LSM_Simulation_<i>/run_sim.py for every combination,
# plus cases.txt and submit_all.sh (a Slurm array job) in the same directory
```

!!! bug "Known issue"
    `prepareSlurmLSMExecution` currently fails with `TypeError`: it calls `LSMTemplate.prepareParams(desc=None, ...)` but the parameter is named `template_desc` (`hera/simulations/LSM/toolkit.py:366`).

## CLI

```bash
# Set up a template that was loaded from a repository (creates the symbolic links to the model code)
hera-LSM setup <templateName> <codeDir> [--templateVersion 1] [--projectName MY_PROJECT]

# List the templates in the project (name, version, path, model folder)
hera-LSM list [--projectName MY_PROJECT]
```

Without `--projectName` the project name is read from `caseConfiguration.json`. Templates are added to a project by loading the repository that contains them (`hera-project project updateRepositories`).

!!! bug "Known issue"
    In master `hera-LSM` fails on every command: `hera/bin/hera-LSM:41` calls `parser.parse_args()()` (the parsed namespace is called), which raises `TypeError`.

For the full API, see the [API Reference](../../developer_guide/api/simulations.md).
