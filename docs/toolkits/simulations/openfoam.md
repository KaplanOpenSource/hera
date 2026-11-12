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
| `of.presentation` | Export to CSV and VTK formats |
| `of.stochasticLagrangian` | Lagrangian particle dispersion solver extension |
| `of.buoyantReactingFoam` | Compressible reactive flow solver extension (block-mesh and hydrostatic-pressure helpers) |

`OFToolkit` inherits from the [Hermes workflow toolkit](workflows.md), so every workflow method (`addWorkflowToGroup`, `listGroups`, `listWorkflows`, `compareWorkflowInGroup`, `executeWorkflowFromDB`, ...) is also available on `of`.

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
#   ├── 0.orig/       (same fields)
#   ├── 0.parallel/   (same fields, with a proc.* processor boundary)
#   ├── constant/
#   │   └── triSurface/
#   └── system/
```

To add a single field to an existing case:

```python
of.writeEmptyField(
    fieldName="U",
    flowType=of.FLOWTYPE_INCOMPRESSIBLE,
    caseDirectory="/data/simulations/wind_study",
    timeOrLocation=0,              # time directory (or e.g. "0.orig")
    readBoundaryFromCase=False,    # True: take the patches from the case mesh
    writeProcBoundary=False,       # True: add the proc.* processor boundary
)
```

Fields that are not in the catalog can be described with `additionalFieldsDescription` — a dict (or JSON file) of `{fieldName: {dimensions: {...}, fieldType: ...}}`, passed to `OFObjectHome.addFieldDefinitions`.

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
#    processorIndex    Ux    Uy    Uz         region  processor
# 0               0  4.92  0.12  0.01  internalField          0
# 1               1  5.01  0.08  0.02  internalField          0
# (boundary rows also carry 'boundary' and 'type'; getDataFrame(filterInternalPatches=True) drops proc* patches)
```

### Reading fields across multiple timesteps as a DataFrame

!!! warning "Deprecated"
    `readFieldAsDataFrame` is marked deprecated in the code — prefer `readFieldFromCase(...).getDataFrame()` per timestep.

```python
# Read velocity at multiple timesteps — useful for time-series analysis
df = of.OFObjectHome.readFieldAsDataFrame(
    fieldName="U",
    caseDirectory="/data/simulations/wind_study",
    times=["100", "200", "300"],   # None = all time directories
    readParallel=True,   # reads from processor* directories (raises ValueError if there are none)
)
# Returns a DataFrame of the internal field with the field components, time (and processor for parallel cases)
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
    fieldName="myConcentration",
    dimensions={"kg": 1, "m": -3},   # exponents of kg, m, s, K, mol, A, cd; {} = dimensionless
    fieldType=of.FIELDTYPE_SCALAR,   # scalar, vector or tensor
    # compressible_dimensions={...}  # if the units differ for compressible flows
    # overwrite=True                 # replace an existing definition (otherwise ValueError)
)
```

---

## Configuring mesh with workflows

Workflows define the complete mesh and solver setup as a JSON-based DAG. Use workflow methods to configure the mesh programmatically.

`getHermesWorkflowFromDB` returns a workflow object whose class is chosen from the `solver` field of the workflow (e.g. `workflow_simpleFoam`, `workflow_buoyantReactingFoam` in `hera.simulations.openFoam.OFWorkflow`). The block-mesh helpers below only accept Eulerian workflows (subclasses of `workflow_Eulerian`) and raise `ValueError` otherwise. They change the workflow object **in memory** — store it back with `of.updateDocumentWorkflow(document, wf)` before running it from the DB.

### Block mesh from coordinates

```python
# Get a workflow for simpleFoam
wf = of.getHermesWorkflowFromDB("wind_study_0001")

# Configure a rectangular single-block mesh directly on the workflow
wf.set_blockMesh_blockBoundaries(
    minx=0, maxx=1000,    # domain extent in x (meters)
    miny=0, maxy=500,     # domain extent in y
    minz=0, maxz=200,     # domain extent in z (height)
    dx=10, dy=10, dz=5,   # cell spacing
)
# This creates a mesh of 100×50×40 = 200,000 cells
```

!!! bug "Known issue"
    The toolkit wrapper `of.buoyantReactingFoam.blockMesh_setBoundFromBounds(eulerianWF, minx, maxx, miny, maxy, minz, maxz, dx, dy, dz)` currently fails with a `NameError` (it passes an undefined name to the logger — `hera/simulations/openFoam/eulerian/abstractEulerianSolver.py:65`). Call `wf.set_blockMesh_blockBoundaries(...)` directly as shown above.

### Block mesh from geometry file

```python
# Set the block-mesh bounds from the bounding box of an OBJ file (requires FreeCAD)
of.buoyantReactingFoam.blockMesh_setBoundFromFile(
    eulerianWorkFlow=wf,
    fileName="/data/terrain/haifa_terrain.obj",
    dx=10, dy=10, dz=5,
)
```

### Adjusting domain height

```python
# Change the domain height after initial mesh setup
of.buoyantReactingFoam.blockMesh_setDomainHeight(
    eulerianWorkFlow=wf,
    Z=300,    # new domain height
    dz=5,     # vertical cell spacing
)
```

### Hydrostatic pressure initialisation

For compressible or buoyant flows, initialise pressure with a hydrostatic profile:

```python
# Compute p = p_ground - 9.81·z for all cells (z = cell-centre height)
# Also handles fixedValue boundary patches
p_field = of.buoyantReactingFoam.IC_getHydrostaticPressure(
    caseDirectory="/data/simulations/thermal_study",
    fieldName="p",
    groundPressure=101000,   # Pa
)
p_field.writeToCase("/data/simulations/thermal_study", timeOrLocation="0")
```

!!! bug "Known issue"
    `IC_getHydrostaticPressure` (and the CLI `hera-openFoam buoyantReactingFoam IC hydrostaticPressure` that wraps it) currently fails with `AttributeError`: it calls `OFObjectHome.getFieldFromCase`, which does not exist (`hera/simulations/openFoam/eulerian/abstractEulerianSolver.py:177`). It also always reads `p` (the `fieldName` argument is ignored) and uses 101000 Pa for the fixedValue patches regardless of `groundPressure`.

---

## Running simulations

### Via workflow (recommended)

A workflow encapsulates the entire simulation: mesh generation, initial/boundary conditions, solver execution.

```python
# Run a complete simulation from a workflow stored in the DB
of.runOFSimulation("wind_study_0001")
# This:
#   1. Retrieves the workflow JSON from MongoDB
#   2. Builds a Luigi task DAG (mesh → IC/BC → solver → postProcess)
#   3. Writes the workflow JSON and the generated Python module
#   4. Executes via: python3 -m luigi --module wind_study_0001 finalnode_xx_0 --local-scheduler
#   5. Cleans up the generated module
#   6. Runs ./Allrun inside the case directory (the document resource)
```

`runOFSimulation` also accepts `scheduler="central"` (with optional `schedulerHost`, `schedulerPort` and `dispatch_id`) to run the Luigi DAG against a running `luigid` instead of the local scheduler.

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
# Generate a Slurm array submission for a parameter sweep
of.prepareSlurmWorkflowExecution(
    baseConfiguration="wind_study_0001",        # name of a workflow already in the DB
    jsonVariations="/path/to/variations.json",  # path or list of variation groups
    jobName="wind_sweep",
    allocateProcessorsPerRun=8,                 # optional
)
# Each variation is added to the group of the base workflow (wind_study_0002, ...),
# and <filesDirectory>/cases.txt + <filesDirectory>/submit_all.sh are written.
# Each Slurm job runs `hera-workflows sync --force` + `hera-workflows buildExecute`
# for its case and then `bash Allrun` (addAllRun=True).
```

The variations file is a **list of variation groups** (see `hera.utils.jsonutils.JSONVariations`): each group maps JSON paths in the workflow to lists of values. Parameters within a group change together (equal-length lists); different groups are combined as a Cartesian product.

```json
[
    {"<json.path.to.parameter>": [3, 5, 8, 12]}
]
```

!!! note
    In master `baseConfiguration` must be a workflow **name**; passing a workflow dict fails because the group name is read from the DB document.

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
# Returns an OFField of the cell centres (field "cellCenters", file "C");
# getDataFrame() gives columns Cx, Cy, Cz, region, processorIndex (plus processor for parallel cases)
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
times = of.getTimeList("wind_study_0001")   # DB name or case directory
# ['0', '100', '200', '300', '400', '500']
```

`getMeshFromName(name)` and `getMeshExtentFromName(name)` do the same as `getMesh` / `getMeshExtent` for a workflow stored in the DB.

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
open. `ITM` and `WSG84` come from `hera.measurements.GIS`.

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

`netcdfToSetFieldsDict(netcdfFile, ...)` is the same for a NetCDF file on disk.

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
    nameOrWorkflowFileOrJSONOrResource="wind_study_0001",  # DB name or directory
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
cached = registered.getData(regularMesh=True, filterName="ground_slice", timeList=["100"])

# Only the last available timestep
latest = registered.getData(regularMesh=True, filterName="ground_slice", latestTime=True)

# Force recomputation
fresh = registered.getData(regularMesh=True, filterName="ground_slice", overwrite=True)
```

`results` is a dict mapping filter names to data (xarray Dataset or pandas DataFrame). `regularMesh` is required. `registered.getRegularData(...)` and `registered.getNonRegularData(...)` are shortcuts for `regularMesh=True` / `False`.

### Available filter types

| Filter | ParaView type | Typed class | Key parameters |
|--------|--------------|-------------|---------------|
| `Slice` | Planar cut | `vtkFilter_Slice` | `SliceType.Origin`, `SliceType.Normal` |
| `PlotOverLine` | Line sample | `vtkFilter_PlotOverLine` | `Point1`, `Point2`, `SamplingPattern`, `Resolution` |
| `CellCenters` | Cell values | `vtkFilter_CellCenters` | — |
| `ExtractBlock` | Patch extraction | `vtkFilter_ExtractBlock` | `Selectors` (XPath: `/Root/boundary/{name}`) |
| `DescriptiveStatistics` | Statistics | `vtkFilter_DescriptiveStatistics` | `VariablesofInterest` |
| `IntegrateVariables` | Field integral | `vtkFilter_IntegrateVariables` | — |

!!! note "Parameter order matters"
    The `params` list is applied in order. Some ParaView filters require properties to be set in a specific sequence (e.g., `SliceType` must be set before `SliceType.Origin`). Always use a list of tuples, not a dict.

### Caching behaviour

- Results are cached per filter in the project's Cache collection (type `vtk_filter`)
- Cache key = pipeline JSON + filter name + simulation parameters
- Incremental: only computes timesteps not already in cache
- Changing the pipeline definition invalidates the cache (different JSON = different key)
- `overwrite=True` forces full recomputation

The cached results can be inspected and cleared from the toolkit:

```python
of.getVTKPipelineCacheTable(workflowName="wind_study_0001")          # DataFrame of cached filter results
of.getVTKPipelineCacheDocuments(filterName="ground_slice")          # the cache documents
of.clearVTKPipelineCache(workflowName="wind_study_0001", filterName="ground_slice")
```

!!! bug "Known issue"
    `registeredVTKPipeLine.clearCache()` queries the Simulations collection instead of the Cache collection (`hera/simulations/openFoam/postProcess/VTKPipeline.py:268`), so it does not find the cached results. Use `of.clearVTKPipelineCache(...)`.

---

## Exporting results

### ParaView CSV (one file per timestep)

`to_paraview_CSV` writes particle positions for ParaView. It expects a pandas/dask DataFrame with `globalX`, `globalY`, `globalZ` and `time` columns (the legacy particle format):

```python
of.presentation.to_paraview_CSV(
    data=particle_df,
    outputdirectory="/data/output/csv",
    filename="wind_study",
    timeFactor=1,          # multiply the time so it becomes a round number
)
# Creates: wind_study_100.csv, wind_study_200.csv, ...
```

!!! note
    The current Lagrangian readers (`getCaseResults`, `loadLagrangianDataParallel`) return `x`, `y`, `z` and `datetime` columns, not `globalX`/`time`, so their output has to be renamed before calling `to_paraview_CSV`.

### Unstructured VTK (for particle data)

```python
of.presentation.toUnstructuredVTK(
    data=particle_data,
    outputdirectory="/data/output/vtk",
    filename="particles",
    xcoord="x", ycoord="y", zcoord="z",   # defaults are globalX/globalY/globalZ
    timecoord="datetime",                 # pandas input is grouped by this column (default "time")
    fieldList=None,                       # None = all non-coordinate columns
    overwrite=False,
)
# Writes particles_0.vtu, particles_1.vtu, ... (one per timestep)
```

For a dask DataFrame the timesteps are taken from the index.

### Structured VTK (for regular grid data)

!!! bug "Known issue"
    `of.presentation.toStructuredVTK(data, outputdirectory, filename, extents, dxdydz, timeNameOutput=True)` currently fails: it calls `self.analysis.calcConcentrationTimeStepFullMesh`, which does not exist on the toolkit analysis layer (it lives on `of.stochasticLagrangian.analysis`) — `hera/simulations/openFoam/toolkit.py:1235`.

### Loading Lagrangian data for export

```python
# Load Lagrangian particle data as a Dask DataFrame (legacy loader, marked "to be debugged" in the code;
# prefer of.stochasticLagrangian.getCaseResults)
particles = of.presentation.loadLagrangianDataParallel(
    casePath="/data/simulations/dispersion_001",
    parallelCase=True,
    withMass=True,
    cloudName="kinematicCloud",
)
# Then export:
of.presentation.toUnstructuredVTK(data=particles, outputdirectory="/data/output/vtk",
                                  filename="particles", xcoord="x", ycoord="y", zcoord="z")
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
        "time": {
            "temporalType": "steadyState",    # or "dynamic" (of.TIME_STEADYSTATE / of.TIME_DYNAMIC)
        },
        "timeStep": 500,                      # flow timestep to use (nearest is taken).
                                              # Default: last (steadyState) / first (dynamic)
        "linkMeshSymbolically": True,         # symlink mesh (saves disk space, default True)
        "linkDataSymbolically": True,         # default True
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
    OriginalFlowField="wind_study_0001",   # DB name, directory or resource of the base flow
    dispersionDuration=3600,   # seconds
    # flowType=of.FLOWTYPE_INCOMPRESSIBLE, overwrite=False, useDBSupport=True
)
# The new flow is named <group>_DFF_<flowName> (here wind_study_DFF_dispersion_001),
# created under the toolkit files directory; `result` is its directory
# (None when useDBSupport=False).
# This:
#   1. Finds the original flow (DB or filesystem)
#   2. Detects parallel structure and available timesteps
#   3. Maps flow timesteps to dispersion time (steady-state: freeze; dynamic: shift)
#   4. Copies/links constant, system, and time directories
#   5. Writes ustar and Hmix fields at each dispersion timestep
#   6. Registers in the database
```

The same operation is available from the CLI as `hera-openFoam stochasticLagrangian dispersionFlow create` (see [CLI commands](#cli-commands)).

### Defining particle sources

```python
# Point source (single release point) — writes <case>/constant/kinematicCloudPositions
of.stochasticLagrangian.writeParticlePositionFile(
    x=500, y=250, z=10,
    nParticles=10000,
    dispersionName="/data/simulations/dispersion_case",   # the case directory
    type="Point",
)

# Cylinder source (e.g. smokestack plume); extra keyword arguments go to makeSource_<type>
of.stochasticLagrangian.writeParticlePositionFile(
    x=500, y=250, z=10,
    nParticles=10000,
    dispersionName="/data/simulations/dispersion_case",
    type="Cylinder",
    radius=5, height=20,
    fileName="kinematicCloudPositions",   # default
)

# The makeSource_* methods only return the particle positions as a DataFrame (x, y, z):
positions = of.stochasticLagrangian.makeSource_Cylinder(x=500, y=250, z=10, nParticles=10000,
                                                         radius=5, height=20)

# Available shapes (of.stochasticLagrangian.sourcesTypeList):
# makeSource_Point(x, y, z, nParticles)
# makeSource_Circle(x, y, z, nParticles, radius, distribution="uniform")
# makeSource_Sphere(x, y, z, nParticles, radius, distribution="uniform")
# makeSource_Cylinder(x, y, z, nParticles, radius, height, horizontalDistribution="uniform", verticalDistribution="uniform")
# makeSource_Rectangle(x, y, z, nParticles, lengthX, lengthY, rotateAngle=0)
# makeSource_Cube(x, y, z, nParticles, lengthX, lengthY, lengthZ, rotateAngle=0)
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
    # timeList=None, forceSingleProcessor=False, overwrite=False
)
# Returns Dask DataFrame with columns: x, y, z, id, procId, globalID, datetime,
# U_x, U_y, U_z (withVelocity), releaseTime, age (withReleaseTimes), mass (withMass)

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
# Returns xarray Dataset with variable C and a dt attribute (cached in the Cache collection)

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
# Intended to return a DataFrame with columns: time, cloudName, name, action (release/escape/stick), parcels, mass
```

!!! bug "Known issue"
    `getMassFromLog` currently raises `AttributeError`: it reads `self._datalayer.cloudName`, which the solver extension does not define (`hera/simulations/openFoam/lagrangian/abstractLagrangianSolver.py:831`).

---

## Templates

Solver templates are Hermes workflow JSON files stored as **data sources** of the toolkit, whose description carries the `solver` they belong to (they are usually loaded from a repository, e.g. with `hera-project project updateRepositories`). There is no save/load-case API; templates are listed and retrieved as data sources:

```python
# List the templates of a solver (DataFrame indexed by template name)
of.listHermesSolverTemplates(solverName="simpleFoam")

# Get a flow template (data source with desc component="Flow")
flowJSON = of.getHermesFlowTemplate("simpleFoam_base")

# Node templates
of.listHermesNodesTemplates()
node = of.getHermesNodeTemplate("blockMesh")
```

To start a new case from a template, write it to disk with `hera-openFoam <solver> templates get` (below), edit it, and add it to the DB with `of.addWorkflowToGroup(...)` or `hera-workflows add`.

---

## CLI commands

`hera-openFoam` is organised by solver. The solver subcommands are `simpleFoam`, `pimpleFoam`, `oguraPhillipShallowConvectionFOAM`, `buoyantReactingFoam`, `scalarTransportFoam`, `homogenousWindLogProfile` and `stochasticLagrangian`. Every solver subcommand has the same core actions (shown here for `simpleFoam`). When `--projectName` is omitted, the project name is read from `caseConfiguration.json`.

```bash
# Case management
hera-openFoam simpleFoam case createEmpty /path/to/case --fields U p k epsilon \
    [--fieldsDescription extraFields.json] [--projectName MY_PROJECT]
hera-openFoam simpleFoam case writeFieldDescription extraFields.json --fields myField

# Templates
hera-openFoam simpleFoam templates list --projectName MY_PROJECT
# Write template <name> to <projectPath>/<groupName>_1.json (defaults: cwd, groupName = template name)
hera-openFoam simpleFoam templates get myTemplate \
    --projectName MY_PROJECT [--groupName wind_study] [--projectPath /path/to/dir]

# Add a workflow file to the DB and build + execute it (requires hermes)
hera-openFoam simpleFoam executeWorkflow wind_study_1.json --projectName MY_PROJECT \
    [--noDB] [--scheduler local|central] [--scheduler-host HOST] [--scheduler-port PORT] [--dispatch-id ID]

# List the groups/workflows of this solver, with a table of the parameters that differ
hera-openFoam simpleFoam list --projectName MY_PROJECT \
    [--longTable] [--transpose] [--format pandas|json|latex|csv] [--file out]

# Delete workflows (optionally exporting them to disk first)
hera-openFoam simpleFoam delete wind_study_0001 wind_study_0002 --projectName MY_PROJECT [--export] [--overwrite]
```

Solver-specific commands:

```bash
# buoyantReactingFoam: hydrostatic initial pressure (currently broken, see above)
hera-openFoam buoyantReactingFoam IC hydrostaticPressure /path/to/case [--startTime 0]
# (buoyantReactingFoam blockMesh setBoundsFromFile / setDomainHeight are registered but not implemented)

# stochasticLagrangian
hera-openFoam stochasticLagrangian dispersionFlow writeTemplate dispersionFlow.json
hera-openFoam stochasticLagrangian dispersionFlow create wind_study_0001 dispersion_001 \
    --parameters dispersionFlow.json --dispersionDuration 3600 [--projectName MY_PROJECT] [--overwrite]
hera-openFoam stochasticLagrangian dispersion create <dispersionName> <dispersionFlowField> [--overwrite]
hera-openFoam stochasticLagrangian injectors positionFile cylinder <dispersionName> \
    --center X Y Z --radius R --height H --particles N
hera-openFoam stochasticLagrangian injectors positionFile createReleaseRateFile
hera-openFoam stochasticLagrangian postProcess toParquet <dispersionName> [--outputDirectory VTK] [--cloudName kinematicCloud] [--overwrite]
hera-openFoam stochasticLagrangian postProcess toVTK <dispersionName> [--outputDirectory VTK] [--cloudName kinematicCloud] [--overwrite]
```

General commands:

```bash
# Build blockMesh vertices/boundaries from an OBJ file
hera-openFoam objects createVerticesAndBoundary terrain.obj --fields U p

# Add a workflow to the DB and build + execute it (same as <solver> executeWorkflow)
hera-openFoam buildExecute wind_study_1.json --projectName MY_PROJECT [--noDB] [--scheduler ...]
```

!!! note
    The `--workflowGroup`, `--overwrite`, `--assignName` and `--allowDuplicate` flags of `executeWorkflow` / `buildExecute` are parsed but currently ignored: the workflow is added with the same logic as `hera-workflows add` (group inferred from the file name `<group>_<id>.json`).

Workflow groups and comparisons across solvers are handled by [`hera-workflows`](workflows.md#cli).

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

# === 2. Configure mesh via workflow (and store the change in the DB) ===
wf = of.getHermesWorkflowFromDB("haifa_wind_0001")
wf.set_blockMesh_blockBoundaries(minx=0, maxx=2000, miny=0, maxy=1000,
                                 minz=0, maxz=300, dx=20, dy=20, dz=10)
of.updateDocumentWorkflow(of.getWorkflowDocumentByName("haifa_wind_0001"), wf)

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
of.runOFSimulation("haifa_wind_0001")

# === 5. Post-process: extract a ground-level slice ===
pipeline = of.analysis.getVTKPipeline()
pipeline.addFilter("ground", filterType="Slice", write=True, params=[
    ("SliceType", "Plane"),
    ("SliceType.Origin", [1000, 500, 2]),
    ("SliceType.Normal", [0, 0, 1]),
])

registered = pipeline.registerPipeline("haifa_wind_0001")
results = registered.getData(
    regularMesh=True, filterName="ground",
    fieldNames=["U", "p"],
)

# === 6. Inspect results ===
ground_data = results["ground"]
print(ground_data)
# xarray Dataset with U (vector), p (scalar) on a regular grid

# === 7. Save for external tools ===
ground_data.to_netcdf("/data/output/haifa_wind_ground.nc")
```

---

## Further reading

- **[Hermes Workflows](workflows.md)** — workflow lifecycle, naming conventions, Luigi integration
- **[Developer Guide > OpenFOAM](../../developer_guide/simulations/openfoam.md)** — architecture, solver hierarchy, VTK pipeline internals, refactored method documentation
- **[Developer Guide > Hermes](../../developer_guide/simulations/workflows.md)** — Hermes workflow toolkit internals
- **[API Reference](../../developer_guide/api/simulations.md)** — auto-generated API docs
- **[CLI Reference](../../cli/reference.md#hera-openfoam)** — full CLI command reference
