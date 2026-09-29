# OpenFOAM

**Toolkit name:** `OpenFOAM` | **Constant:** `toolkitHome.SIMULATIONS_OPENFOAM`

The OpenFOAM toolkit manages the full CFD simulation lifecycle: creating cases, configuring meshes and boundary conditions, running solvers, extracting results via VTK post-processing pipelines, and computing Lagrangian dispersion concentrations.

```python
from hera import toolkitHome

# Tip: if you created the project with `hera-project project create`, you can omit projectName
of = toolkitHome.getToolkit(toolkitHome.SIMULATIONS_OPENFOAM, projectName="MY_PROJECT")
```

The toolkit provides access to:

| Attribute | What it does |
|-----------|-------------|
| `of.OFObjectHome` | Field catalog — create empty fields with correct dimensions |
| `of.analysis` | VTK post-processing pipeline |
| `of.presentation` | Export to CSV, VTK, parquet formats |
| `of.stochasticLagrangian` | Lagrangian particle dispersion solver extension |
| `of.buoyantReactingFoam` | Compressible reactive flow solver extension |

For implementation details, see the [Developer Guide](../../developer_guide/simulations/openfoam.md). For workflow management, see [Hermes Workflows](workflows.md).

---

## Creating an empty case

An OpenFOAM case is a directory containing `0/` (initial conditions), `constant/` (mesh and physical properties), and `system/` (solver settings). The toolkit creates this structure and populates the field files with correct dimensions.

```python
# Create a case with velocity, pressure, and turbulence fields
of.createEmptyCase(
    caseDirectory="/data/simulations/wind_study",
    fieldList=["U", "p", "k", "epsilon"],
    flowType=of.FLOWTYPE_INCOMPRESSIBLE,
)
# Creates:
#   wind_study/
#   ├── 0/
#   │   ├── U        (vector, dimensions [0 1 -1 0 0 0 0])
#   │   ├── p        (scalar, dimensions [0 2 -2 0 0 0 0] for incompressible)
#   │   ├── k        (scalar, turbulent kinetic energy)
#   │   └── epsilon   (scalar, dissipation rate)
#   ├── constant/
#   └── system/
```

**Flow types** control which dimensions are used for pressure and other fields:

| Constant | Pressure dimensions | Use case |
|----------|-------------------|---------|
| `of.FLOWTYPE_INCOMPRESSIBLE` | m²/s² (kinematic) | simpleFoam, wind simulations |
| `of.FLOWTYPE_COMPRESSIBLE` | kg/(m·s²) (dynamic) | buoyantReactingFoam, thermal flows |
| `of.FLOWTYPE_DISPERSION` | — | Lagrangian particle tracking |

---

## Working with fields (OFObjectHome)

The field catalog knows the correct dimensions for standard OpenFOAM fields. You can create, read, modify, and write fields programmatically.

### Creating a field from scratch

```python
# Get an empty velocity field with correct dimensions
U_field = of.OFObjectHome.getEmptyField(
    fieldName="U",
    flowType=of.FLOWTYPE_INCOMPRESSIBLE,
    noOfProc=4,   # number of processors (for parallel decomposed cases)
)

# Set the initial value for all cells
U_field.setInternalUniformFieldValue([5, 0, 0])  # 5 m/s in x direction

# Set boundary conditions
U_field.addBoundaryField("inlet", type="fixedValue", value="uniform (5 0 0)")
U_field.addBoundaryField("outlet", type="zeroGradient")
U_field.addBoundaryField("ground", type="noSlip")
U_field.addBoundaryField("top", type="slip")

# Write to the case directory at time 0
U_field.writeToCase("/data/simulations/wind_study", timeOrLocation="0")
```

### Creating a field pre-loaded with boundaries from an existing case

If you have an existing case and want to create a new field that matches its boundary structure:

```python
# Creates a field with the same patches as the case, plus a custom internal value
ustar_field = of.OFObjectHome.getEmptyFieldFromCase(
    fieldName="ustar",
    flowType=of.FLOWTYPE_INCOMPRESSIBLE,
    internalValue=0.3,      # uniform initial value
    caseDirectory="/data/simulations/wind_study",
)
ustar_field.writeToCase("/data/simulations/wind_study", timeOrLocation="100")
```

### Reading fields from a completed simulation

```python
# Read velocity at a specific timestep
U_result = of.OFObjectHome.readFieldFromCase(
    fieldName="U",
    flowType=of.FLOWTYPE_INCOMPRESSIBLE,
    caseDirectory="/data/simulations/wind_study",
    timeStep="100",
)

# Access data
internal_field = U_result.getDataFrame()
print(internal_field.head())
#     Ux    Uy    Uz    region      processor  index
# 0   4.92  0.12  0.01  internalField  0       0
# 1   5.01  0.08  0.02  internalField  0       1
```

### Reading fields across multiple timesteps as a DataFrame

```python
# Read velocity at multiple timesteps — useful for time-series analysis
df = of.OFObjectHome.readFieldAsDataFrame(
    fieldName="U",
    caseDirectory="/data/simulations/wind_study",
    times=["100", "200", "300"],
    readParallel=True,   # reads from processor* directories
)
# Returns DataFrame with columns: Ux, Uy, Uz, processor, time, index
```

### Available predefined fields

| Field | Type | Description |
|-------|------|-------------|
| `U` | vector | Velocity |
| `p`, `p_rgh` | scalar | Pressure (kinematic or with gravity) |
| `k` | scalar | Turbulent kinetic energy |
| `epsilon`, `omega` | scalar | Dissipation rate / specific dissipation |
| `nut`, `alphat` | scalar | Turbulent viscosity / thermal diffusivity |
| `T`, `Tbackground` | scalar | Temperature |
| `cellCenters` | vector | Cell center coordinates (computed by postProcess) |
| `Hmix`, `ustar` | scalar | Mixing height / friction velocity (for dispersion) |
| `distanceFromWalls` | scalar | Wall distance field |

You can register custom fields:

```python
of.OFObjectHome.addFieldDefinitions(
    fieldName="myScalar",
    dimensions={"default": [0, 0, 0, 0, 0, 0, 0]},
    fieldType="scalar",
)
```

---

## Configuring mesh with workflows

Workflows define the complete mesh and solver setup as a JSON-based DAG. Use workflow methods to configure the mesh programmatically.

### Block mesh from coordinates

```python
# Get a workflow for simpleFoam
wf = of.getHermesWorkflowFromDB("wind_study_simpleFoam")

# Configure a rectangular block mesh
of.buoyantReactingFoam.blockMesh_setBoundFromBounds(
    eulerianWF=wf,
    minx=0, maxx=1000,    # domain extent in x (meters)
    miny=0, maxy=500,     # domain extent in y
    minz=0, maxz=200,     # domain extent in z (height)
    dx=10, dy=10, dz=5,   # cell spacing
)
# This creates a mesh of 100×50×40 = 200,000 cells
```

### Block mesh from geometry file

```python
# Configure mesh from a terrain OBJ file (the mesh follows the terrain surface)
of.buoyantReactingFoam.blockMesh_setBoundFromFile(
    eulerianWF=wf,
    fileName="/data/terrain/haifa_terrain.obj",
    dx=10, dy=10, dz=5,
)
```

### Adjusting domain height

```python
# Change the domain height after initial mesh setup
of.buoyantReactingFoam.blockMesh_setDomainHeight(
    eulerianWF=wf,
    Z=300,    # new domain height
    dz=5,     # vertical cell spacing
)
```

### Hydrostatic pressure initialisation

For compressible or buoyant flows, initialise pressure with a hydrostatic profile:

```python
# Compute p = p_ground - ρg·z for all cells
# Also handles fixedValue boundary patches
p_field = of.buoyantReactingFoam.IC_getHydrostaticPressure(
    caseDirectory="/data/simulations/thermal_study",
    fieldName="p",
    groundPressure=101000,   # Pa
)
p_field.writeToCase("/data/simulations/thermal_study", timeOrLocation="0")
```

---

## Running simulations

### Via workflow (recommended)

A workflow encapsulates the entire simulation: mesh generation, initial/boundary conditions, solver execution.

```python
# Run a complete simulation from a workflow stored in the DB
of.runOFSimulation("wind_study_simpleFoam")
# This:
#   1. Retrieves the workflow JSON from MongoDB
#   2. Builds a Luigi task DAG (mesh → IC/BC → solver → postProcess)
#   3. Writes the generated Python module
#   4. Executes via: python3 -m luigi --module wind_study finalnode_xx_0 --local-scheduler
#   5. Cleans up the generated module
```

### Adding a workflow to the database

```python
# Load a workflow JSON and add it to a named group
doc = of.addWorkflowToGroup(
    workflowJSON="/path/to/workflow.json",
    groupName="wind_study",
)
# Generates name: wind_study_0001 (auto-incremented)
print(f"Workflow stored as: {doc.desc['workflowName']}")
```

### Listing and comparing workflows

```python
# List all workflow groups in the project
of.listGroups()

# List workflows in a specific group
of.listWorkflows("wind_study", listNodes=True, listParameters=True)

# Compare parameters across workflows in a group
diff_table = of.compareWorkflowInGroup("wind_study")
print(diff_table)
# Shows a DataFrame with rows = parameters that differ, columns = workflow names
```

### Batch runs with Slurm

```python
# Generate Slurm submission scripts for a parameter sweep
of.prepareSlurmWorkflowExecution(
    workflowName="simpleFoam_sweep",
    variations={"windSpeed": [3, 5, 8, 12]},   # creates 4 cases
    jobName="wind_sweep",
)
# Generates one Slurm .sh script per variation
```

---

## Reading mesh and simulation data

### Mesh cell centers

```python
# Read mesh cell centers (internally runs foamJob postProcess -func writeCellCentres)
mesh = of.getMesh(
    caseDirectory="/data/simulations/wind_study",
    readParallel=True,   # reads decomposed case if processor* dirs exist
    time=0,
)
# Returns DataFrame with columns: x, y, z (and processorNumber, index for parallel)
print(mesh.getDataFrame().head())
```

### Mesh bounding box

```python
extent = of.getMeshExtent(caseDirectory="/data/simulations/wind_study")
print(extent)
# {'x': (0.0, 1000.0), 'y': (0.0, 500.0), 'z': (0.0, 200.0)}
```

### Available timesteps

```python
times = of.getTimeList("wind_study_simpleFoam")
# ['0', '100', '200', '300', '400', '500']
```

### Setting initial conditions from a meteorological dataset

The output of an atmospheric model (GFS/GDAS NetCDF, WRF, or anything that
loads as an `xarray.Dataset`) can be used as the initial and boundary
condition of a case — the meso-scale-to-CFD nesting step. There are two
routes, and both read and write the case with
[foamlib](https://github.com/gerlero/foamlib).

**Route 1 — write the values straight into the fields of the case.**
The dataset is interpolated onto the cell centres (bilinear in the
horizontal, linear in the vertical) and onto the face centres of the
requested patches:

```python
# The cell centres must exist; getMesh runs writeCellCentres for you.
of.getMesh("/path/to/case")

written = of.netcdfToCaseFields(
    caseDirectory="/path/to/case",
    netcdfFile="forecast.nc",
    fieldMap=dict(U=("ugrd", "vgrd", 0), T="tmp"),  # tuple = vector, string = scalar
    xCoordinate="grid_xt", yCoordinate="grid_yt",   # names of the dataset coordinates
    verticalCoordinate="pfull",
    heightVariable="height",     # height [m] per level, for a terrain-following coordinate
    time="2024-03-15T12:00",     # the nearest time step is selected
    caseCRS=ITM,                 # the coordinate system of the case
    datasetCRS=WSG84,            # ... and of the dataset
    patchList=["inlet", "outlet"],
)
```

The field files must already exist — they carry the type, the dimensions and
the boundary conditions; create them with `createEmptyCase` if needed. A
decomposed case is handled by writing every `processor*` directory.
`datasetToCaseFields` is the same function for a dataset that is already
open.

**Route 2 — emit a `setFieldsDict`.** Every cell of the dataset becomes a
`boxToCell` (and a `boxToFace`) region whose box reaches halfway to the
neighbouring cells, to be applied afterwards by OpenFOAM's `setFields`:

```python
of.datasetToSetFieldsDict(
    dataset=weather,
    fieldMap=dict(U=("Ue", "Ve", "W"), T="Temp"),
    xCoordinate="XLONG", yCoordinate="XLAT",   # 1D or 2D (curvilinear) coordinates
    verticalCoordinate="level",
    heightVariable="gpt_hgt_M",
    extent=dict(minX=735000, maxX=754000, minY=193100, maxY=222800, minZ=-60, maxZ=2690),
    defaultFieldValues=["volScalarFieldValue", "T", 300],
    outputFile="/path/to/case/system/setFieldsDict",
)
```

Cells whose centre falls outside `extent` are dropped and the outermost
boxes are clipped to it.

`xarrayToSetFieldsDictDomain` is the older, narrower variant: it returns the
regions as text only, treats the coordinates of the dataset as cell *edges*
already in the coordinate system of the case, and does not write a file.
Prefer `datasetToSetFieldsDict`.

---

## Post-processing with VTK pipeline

The analysis layer provides a ParaView-integrated VTK filter pipeline for extracting data from simulation results, with automatic DB caching.

### Creating a pipeline

```python
# 1. Create an empty pipeline
pipeline = of.analysis.getVTKPipeline()

# 2. Add filters — params is a LIST OF TUPLES (order matters for ParaView)
pipeline.addFilter("ground_slice", filterType="Slice", write=True,
    params=[
        ("SliceType", "Plane"),
        ("SliceType.Origin", [500, 250, 2]),   # dot notation for nested properties
        ("SliceType.Normal", [0, 0, 1]),
    ]
)

pipeline.addFilter("centerline", filterType="PlotOverLine", write=True,
    params=[
        ("Point1", [0, 250, 10]),
        ("Point2", [1000, 250, 10]),
        ("SamplingPattern", "Sample Uniformly"),
        ("Resolution", 100),
    ]
)
```

### Chaining filters (downstream)

Filters can be chained — each filter receives the output of its parent:

```python
# Add an ExtractBlock filter, then chain CellCenters downstream
extract = pipeline.addFilter("boundary", filterType="ExtractBlock", write=False)
extract.setRegionsToExtract(patchList=["inlet", "outlet"], internalMesh=False)

# Chain: boundary → cell_centers (cell_centers receives ExtractBlock output)
extract.addFilter("boundary_values", filterType="CellCenters", write=True)
```

### Using typed filter constructors

For convenience, filter subclasses provide typed methods:

```python
from hera.simulations.openFoam.postProcess.VTKPipeline import vtkFilter_Slice

slice_filter = vtkFilter_Slice(name="ground_slice", write=True)
slice_filter.setPlaneOrigin([500, 250, 2])
slice_filter.setPlaneNormal([0, 0, 1])
pipeline.addFilterFromObj(slice_filter)
```

### Registering with a simulation case

Bind the pipeline to a specific OpenFOAM case:

```python
registered = pipeline.registerPipeline(
    nameOrWorkflowFileOrJSONOrResource="wind_study_simpleFoam",  # DB name or directory
    caseType=of.CASETYPE_RECONSTRUCTED,   # or of.CASETYPE_DECOMPOSED for parallel
)
```

### Executing and retrieving results

```python
# Execute — computes only missing timesteps (incremental caching)
results = registered.getData(
    regularMesh=True,           # True: xarray/zarr, False: pandas/parquet
    filterName="ground_slice",  # or list, or None for all write=True filters
    timeList=["100", "200"],    # or None (all), or "100:500" (range)
    fieldNames=["U", "p"],      # restrict which OF fields are read
)

# Subsequent calls load from cache — no recomputation
cached = registered.getData(filterName="ground_slice", timeList=["100"])

# Force recomputation
fresh = registered.getData(filterName="ground_slice", overwrite=True)
```

`results` is a dict mapping filter names to data (xarray Dataset or pandas DataFrame).

### Available filter types

| Filter | ParaView type | Typed class | Key parameters |
|--------|--------------|-------------|---------------|
| `Slice` | Planar cut | `vtkFilter_Slice` | `SliceType.Origin`, `SliceType.Normal` |
| `PlotOverLine` | Line sample | `vtkFilter_PlotOverLine` | `Point1`, `Point2`, `SamplingPattern`, `Resolution` |
| `CellCenters` | Cell values | `vtkFilter_CellCenters` | — |
| `ExtractBlock` | Patch extraction | `vtkFilter_ExtractBlock` | `Selectors` (XPath: `/Root/boundary/{name}`) |
| `DescriptiveStatistics` | Statistics | `vtkFilter_DescriptiveStatistics` | `VariablesOfInterest` |
| `IntegrateVariables` | Field integral | `vtkFilter_IntegrateVariables` | — |

!!! note "Parameter order matters"
    The `params` list is applied in order. Some ParaView filters require properties to be set in a specific sequence (e.g., `SliceType` must be set before `SliceType.Origin`). Always use a list of tuples, not a dict.

### Caching behaviour

- Results are cached per filter in the project's Cache collection (type `vtk_filter`)
- Cache key = pipeline JSON + filter name + simulation parameters
- Incremental: only computes timesteps not already in cache
- Changing the pipeline definition invalidates the cache (different JSON = different key)
- `overwrite=True` forces full recomputation

---

## Exporting results

### ParaView CSV (one file per timestep)

```python
of.presentation.to_paraview_CSV(
    data=results["ground_slice"],
    outputdirectory="/data/output/csv",
    filename="wind_study",
)
# Creates: wind_study_100.csv, wind_study_200.csv, ...
```

### Unstructured VTK (for particle data)

```python
of.presentation.toUnstructuredVTK(
    data=particle_data,
    outputdirectory="/data/output/vtk",
    filename="particles",
)
```

### Structured VTK (for regular grid data)

```python
of.presentation.toStructuredVTK(
    data=concentration_grid,
    outputdirectory="/data/output/vtk",
    filename="concentration",
)
```

### Loading Lagrangian data for export

```python
# Load parallel Lagrangian particle data as a Dask DataFrame
particles = of.presentation.loadLagrangianDataParallel(
    caseDirectory="/data/simulations/dispersion_001",
    cloudName="kinematicCloud",
)
# Then export:
of.presentation.toUnstructuredVTK(data=particles, ...)
```

---

## Lagrangian dispersion (stochastic solver)

The stochastic Lagrangian extension manages particle dispersion simulations coupled with an existing flow field.

### Setting up a dispersion flow field

Before running a dispersion simulation, you need to create a dispersion flow field — a copy of the original flow with additional fields (ustar, Hmix) and time mapping:

```python
# flowData describes the base flow and dispersion configuration
flow_data = {
    "originalFlow": {
        "source": "wind_study_simpleFoam",   # name in DB, or directory path
        "time": {
            "temporalType": "steadyState",    # or "dynamic"
            "timestep": 500,                  # which flow timestep to use
        },
        "linkMeshSymbolically": True,         # symlink mesh (saves disk space)
        "linkDataSymbolically": True,
    },
    "dispersionFields": {
        "ustar": 0.3,    # friction velocity (uniform initial value)
        "Hmix": 500,      # mixing height
    },
}

# Create the dispersion flow field
result = of.stochasticLagrangian.createDispersionFlowField(
    flowName="dispersion_001",
    flowData=flow_data,
    OriginalFlowField="wind_study_simpleFoam",
    dispersionDuration=3600,   # seconds
)
# This:
#   1. Finds the original flow (DB or filesystem)
#   2. Detects parallel structure and available timesteps
#   3. Maps flow timesteps to dispersion time (steady-state: freeze; dynamic: shift)
#   4. Copies/links constant, system, and time directories
#   5. Writes ustar and Hmix fields at each dispersion timestep
#   6. Registers in the database
```

### Defining particle sources

```python
# Point source (single release point)
of.stochasticLagrangian.writeParticlePositionFile(
    x=500, y=250, z=10,
    nParticles=10000,
    dispersionName="dispersion_001",
    type="Point",
)

# Cylinder source (e.g. smokestack plume)
of.stochasticLagrangian.makeSource_Cylinder(
    x=500, y=250, z=10,
    nParticles=10000,
    radius=5, height=20,
)

# Other source shapes:
# makeSource_Circle(x, y, z, nParticles, radius)
# makeSource_Sphere(x, y, z, nParticles, radius)
# makeSource_Rectangle(x, y, z, nParticles, lengthX, lengthY, rotateAngle)
# makeSource_Cube(x, y, z, nParticles, lengthX, lengthY, lengthZ, rotateAngle)
```

### Reading dispersion results

```python
# Read Lagrangian particle data (with caching)
particles = of.stochasticLagrangian.getCaseResults(
    caseDescriptor="dispersion_001",
    withVelocity=True,
    withReleaseTimes=True,
    withMass=True,
    cloudName="kinematicCloud",
    cache=True,           # save to parquet cache
)
# Returns Dask DataFrame with columns: x, y, z, Ux, Uy, Uz, mass, age, datetime

# Read Eulerian concentration fields
concentrations = of.stochasticLagrangian.getCaseConcentrationsEulerian(
    caseDescriptor="dispersion_001",
    cloudName="kinematicCloud",
)
# Returns xarray Dataset indexed by (datetime, x, y, z)
```

### Computing concentration on a regular mesh

```python
# Compute concentration field on a Cartesian grid
conc = of.stochasticLagrangian.analysis.calcConcentrationFieldFullMesh(
    caseDescriptor="dispersion_001",
    dxdydz=10,        # grid cell size in meters
    extents=None,      # None = auto from flow field mesh extent
)
# Returns xarray Dataset with variable C (concentration in kg/m³)

# Or compute point-wise from particle positions
conc_pw = of.stochasticLagrangian.analysis.calcConcentrationPointWise(
    data=particles,
    dxdydz=10,
)
```

### Parsing solver log files

```python
# Extract mass and parcel fate information from a dispersion solver log
mass_df = of.stochasticLagrangian.getMassFromLog(
    logFile="/data/simulations/dispersion_001/log.StochasticLagrangianSolver",
)
# Returns DataFrame with columns: time, name, action (release/escape/stick), mass, parcels
```

---

## Templates

Save and reuse case configurations:

```python
# Save a case directory as a reusable template
of.saveTemplate(
    templateName="simpleFoam_base",
    caseDirectory="/data/simulations/wind_study",
)

# List available templates
of.listHermesSolverTemplates(solverName="simpleFoam")

# Load a template into a new directory
of.loadTemplate(
    templateName="simpleFoam_base",
    toDirectory="/data/simulations",
    caseName="wind_study_v2",
)
```

---

## CLI commands

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

# List workflows in a group
hera-openFoam simpleFoam workflows list workflows wind_study --projectName MY_PROJECT

# Compare workflows in a group
hera-openFoam simpleFoam workflows compare wind_study --projectName MY_PROJECT
```

---

## Complete example: wind simulation end-to-end

```python
from hera import toolkitHome

# Tip: if you created the project with `hera-project project create`, you can omit projectName
of = toolkitHome.getToolkit(toolkitHome.SIMULATIONS_OPENFOAM, projectName="WindStudy")

# === 1. Create case ===
case_dir = "/data/simulations/haifa_wind"
of.createEmptyCase(
    caseDirectory=case_dir,
    fieldList=["U", "p", "k", "epsilon"],
    flowType=of.FLOWTYPE_INCOMPRESSIBLE,
)

# === 2. Configure mesh via workflow ===
wf = of.getHermesWorkflowFromDB("haifa_wind_simpleFoam")
of.buoyantReactingFoam.blockMesh_setBoundFromBounds(
    eulerianWF=wf, minx=0, maxx=2000, miny=0, maxy=1000,
    minz=0, maxz=300, dx=20, dy=20, dz=10,
)

# === 3. Set initial conditions ===
U = of.OFObjectHome.getEmptyFieldFromCase(
    fieldName="U", flowType=of.FLOWTYPE_INCOMPRESSIBLE,
    internalValue=[5, 0, 0], caseDirectory=case_dir,
)
U.addBoundaryField("inlet", type="fixedValue", value="uniform (5 0 0)")
U.addBoundaryField("outlet", type="zeroGradient")
U.addBoundaryField("ground", type="noSlip")
U.writeToCase(case_dir, timeOrLocation="0")

# === 4. Run the simulation ===
of.runOFSimulation("haifa_wind_simpleFoam")

# === 5. Post-process: extract a ground-level slice ===
pipeline = of.analysis.getVTKPipeline()
pipeline.addFilter("ground", filterType="Slice", write=True, params=[
    ("SliceType", "Plane"),
    ("SliceType.Origin", [1000, 500, 2]),
    ("SliceType.Normal", [0, 0, 1]),
])

registered = pipeline.registerPipeline("haifa_wind_simpleFoam")
results = registered.getData(
    regularMesh=True, filterName="ground",
    fieldNames=["U", "p"],
)

# === 6. Inspect results ===
ground_data = results["ground"]
print(ground_data)
# xarray Dataset with U (vector), p (scalar) on a regular grid

# === 7. Export to CSV for external tools ===
of.presentation.to_paraview_CSV(
    data=ground_data,
    outputdirectory="/data/output",
    filename="haifa_wind_ground",
)
```

---

## Further reading

- **[Hermes Workflows](workflows.md)** — workflow lifecycle, naming conventions, Luigi integration
- **[Developer Guide > OpenFOAM](../../developer_guide/simulations/openfoam.md)** — architecture, solver hierarchy, VTK pipeline internals, refactored method documentation
- **[Developer Guide > Hermes](../../developer_guide/simulations/workflows.md)** — Hermes workflow toolkit internals
- **[API Reference](../../developer_guide/api/simulations.md)** — auto-generated API docs
- **[CLI Reference](../../cli/reference.md#hera-openfoam)** — full CLI command reference
