# Simulations Toolkits

Simulation toolkits manage the lifecycle of computational models — setting up cases, running simulations, and post-processing results.

All simulation toolkits are accessed via `toolkitHome.getToolkit()` and bound to a project:

```python
from hera import toolkitHome

# Tip: if you created the project with `hera-project project create`, you can omit projectName
toolkit = toolkitHome.getToolkit(toolkitHome.SIMULATIONS_OPENFOAM, projectName="MY_PROJECT")
```

---

## OpenFOAM

**Toolkit name:** `OpenFOAM` | **Constant:** `toolkitHome.SIMULATIONS_OPENFOAM`

The OpenFOAM toolkit manages the full CFD simulation lifecycle: creating cases, configuring meshes and boundary conditions, running solvers, and extracting results via VTK post-processing pipelines.

### Getting started

```python
from hera import toolkitHome

of = toolkitHome.getToolkit(toolkitHome.SIMULATIONS_OPENFOAM, projectName="MY_PROJECT")
```

The toolkit provides access to:

| Attribute | What it does |
|-----------|-------------|
| `of.OFObjectHome` | Field catalog — create empty fields with correct dimensions |
| `of.analysis` | VTK post-processing pipeline |
| `of.presentation` | Export to CSV, VTK, parquet formats |
| `of.stochasticLagrangian` | Lagrangian particle dispersion solver |
| `of.buoyantReactingFoam` | Compressible reactive flow solver |

### Creating an empty case

```python
# Create a new OpenFOAM case with specified fields
of.createEmptyCase(
    caseDirectory="/data/simulations/wind_study",
    fieldList=["U", "p", "k", "epsilon"],
    flowType=of.FLOWTYPE_INCOMPRESSIBLE
)

# Write an empty field file (with correct dimensions for the flow type)
of.writeEmptyField(
    fieldName="U",
    flowType=of.FLOWTYPE_INCOMPRESSIBLE,
    caseDirectory="/data/simulations/wind_study",
    timeOrLocation="0"
)
```

**Flow types:**

| Constant | Use case |
|----------|---------|
| `of.FLOWTYPE_INCOMPRESSIBLE` | simpleFoam, wind simulations |
| `of.FLOWTYPE_COMPRESSIBLE` | buoyantReactingFoam, thermal flows |
| `of.FLOWTYPE_DISPERSION` | Lagrangian particle tracking |

### Working with fields (OFObjectHome)

The field catalog knows the correct dimensions for standard OpenFOAM fields across flow types:

```python
# Get an empty velocity field for an incompressible case
U_field = of.OFObjectHome.getEmptyField(
    fieldName="U",
    flowType=of.FLOWTYPE_INCOMPRESSIBLE,
    noOfProc=4  # number of processors for parallel cases
)

# Set a uniform internal field value
U_field.setInternalUniformFieldValue([5, 0, 0])

# Add boundary conditions
U_field.addBoundaryField("inlet", type="fixedValue", value="uniform (5 0 0)")
U_field.addBoundaryField("outlet", type="zeroGradient")
U_field.addBoundaryField("walls", type="noSlip")

# Write to the case directory
U_field.writeToCase("/data/simulations/wind_study", timeOrLocation="0")

# Read a field from an existing case
U_existing = of.OFObjectHome.readFieldFromCase(
    fieldName="U",
    flowType=of.FLOWTYPE_INCOMPRESSIBLE,
    caseDirectory="/data/simulations/wind_study",
    timeStep="100"
)

# Convert to pandas DataFrame for analysis
df = of.OFObjectHome.readFieldAsDataFrame(
    fieldName="U",
    caseDirectory="/data/simulations/wind_study",
    times=["100", "200", "300"]
)
```

**Predefined fields:** U, p, p_rgh, epsilon, omega, alphat, nut, k, T, Tbackground, cellCenters, Hmix, ustar, distanceFromWalls

### Configuring mesh with workflows

Use workflows to set up the mesh and solver configuration:

```python
# Configure block mesh boundaries from coordinates
of.simpleFoam.blockMesh_setBoundFromBounds(
    eulerianWF=workflow,
    minx=0, maxx=1000,
    miny=0, maxy=500,
    minz=0, maxz=200,
    dx=10, dy=10, dz=5
)

# Or from an OBJ geometry file
of.simpleFoam.blockMesh_setBoundFromFile(
    eulerianWF=workflow,
    fileName="/data/terrain.obj",
    dx=10, dy=10, dz=5
)

# Adjust domain height
of.simpleFoam.blockMesh_setDomainHeight(
    eulerianWF=workflow,
    Z=300, dz=5
)
```

### Running simulations

#### Via workflow (recommended)

```python
# Run a full simulation defined by a workflow JSON
of.runOFSimulation("wind_study_simpleFoam")

# The workflow handles: mesh generation, IC/BC setup, solver execution
```

#### Batch runs with Slurm

```python
# Generate Slurm scripts for parameter sweep
of.prepareSlurmWorkflowExecution(
    workflowName="simpleFoam_sweep",
    variations={"windSpeed": [3, 5, 8, 12]},
    jobName="wind_sweep"
)
```

### Reading mesh data

```python
# Get cell centers (runs foamJob postProcess internally)
mesh = of.getMesh(caseDirectory="/data/simulations/wind_study")

# Get mesh bounding box
extent = of.getMeshExtent(caseDirectory="/data/simulations/wind_study")

# Get list of computed time steps
times = of.getTimeList("wind_study_simpleFoam")
```

### Post-processing with VTK pipeline

The analysis layer provides a VTK filter pipeline for extracting data from simulation results:

```python
# Create a VTK pipeline
pipeline = of.analysis.getVTKPipeline()

# Add filters to extract data
pipeline.addFilter("ground_slice", filterType="Slice",
    write=True,
    params={"origin": [500, 250, 2], "normal": [0, 0, 1]}
)

pipeline.addFilter("centerline", filterType="PlotOverLine",
    write=True,
    params={"point1": [0, 250, 10], "point2": [1000, 250, 10]}
)

pipeline.addFilter("cell_data", filterType="CellCenters",
    write=True
)

# Register the pipeline with a specific simulation case
registered = pipeline.registerPipeline(
    nameOrWorkflowFileOrJSONOrResource="wind_study_simpleFoam",
    serverName="local",
    caseType=of.CASETYPE_RECONSTRUCTED
)

# Execute and get results (with DB caching)
results = registered.getData(
    regularMesh=True,      # True: xarray (regular grid), False: pandas (unstructured)
    filterName="ground_slice",
    timeList=["100", "200"],
    fieldNames=["U", "p"]
)

# Results are cached — subsequent calls load from DB
cached = registered.getData(filterName="ground_slice", timeList=["100"])
```

**Available filter types:**

| Filter | What it does | Key parameters |
|--------|-------------|---------------|
| `Slice` | Planar cut through domain | `origin`, `normal` |
| `PlotOverLine` | Sample values along a line | `point1`, `point2`, sample pattern |
| `CellCenters` | Extract cell center values | — |
| `ExtractBlock` | Extract specific mesh patches | `patchList`, `internalMesh` |
| `DescriptiveStatistics` | Summary statistics | `variables` |
| `IntegrateVariables` | Integrate fields | — |

### Exporting results

```python
# Export as ParaView-compatible CSV (one file per time step)
of.presentation.to_paraview_CSV(
    data=results,
    outputdirectory="/data/output",
    filename="wind_study"
)

# Export particle data as VTK (unstructured)
of.presentation.toUnstructuredVTK(
    data=particle_data,
    outputdirectory="/data/output",
    filename="particles"
)

# Export regular grid as structured VTK
of.presentation.toStructuredVTK(
    data=grid_data,
    outputdirectory="/data/output",
    filename="concentration"
)
```

### Lagrangian dispersion (stochastic solver)

For particle dispersion simulations coupled with an existing flow field:

```python
# Create a dispersion flow field linked to an existing flow
of.stochasticLagrangian.createDispersionFlowField(
    flowName="dispersion_001",
    flowData=flow_data,
    OriginalFlowField="wind_study_simpleFoam",
    dispersionDuration=3600,
    flowType=of.FLOWTYPE_DISPERSION
)

# Define particle sources
of.stochasticLagrangian.makeSource_Cylinder(
    x=500, y=250, z=10,
    nParticles=10000,
    radius=5, height=20,
    dispersionName="dispersion_001",
    fileName="sources.txt"
)

# Also available: makeSource_Point, makeSource_Circle,
# makeSource_Sphere, makeSource_Rectangle, makeSource_Cube
```

### CLI commands

```bash
# List templates
hera-openFoam simpleFoam templates list --projectName MY_PROJECT

# Save a template
hera-openFoam simpleFoam templates save myTemplate \
    --projectName MY_PROJECT --directory /path/to/case

# Load a template
hera-openFoam simpleFoam templates load myTemplate \
    --projectName MY_PROJECT --toDirectory /path/to/output --caseName run_001

# List workflow groups
hera-openFoam simpleFoam workflows list groups --projectName MY_PROJECT
```

### Further reading

- **[Simulations API Reference](../developer_guide/api/simulations.md)** — auto-generated API docs for all OpenFOAM classes
- **[Simulations Implementation](../developer_guide/simulations.md)** — developer guide with architecture, solver hierarchy, VTK pipeline internals
- **[CLI Reference](../cli/reference.md#hera-openfoam)** — full CLI command reference

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
