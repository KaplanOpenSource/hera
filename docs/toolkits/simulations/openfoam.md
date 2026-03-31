# OpenFOAM

**Toolkit name:** `OpenFOAM` | **Constant:** `toolkitHome.SIMULATIONS_OPENFOAM`

The OpenFOAM toolkit manages the full CFD simulation lifecycle: creating cases, configuring meshes and boundary conditions, running solvers, and extracting results via VTK post-processing pipelines.

## Getting started

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

## Creating an empty case

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

## Working with fields (OFObjectHome)

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

## Configuring mesh with workflows

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

## Running simulations

### Via workflow (recommended)

```python
# Run a full simulation defined by a workflow JSON
of.runOFSimulation("wind_study_simpleFoam")

# The workflow handles: mesh generation, IC/BC setup, solver execution
```

### Batch runs with Slurm

```python
# Generate Slurm scripts for parameter sweep
of.prepareSlurmWorkflowExecution(
    workflowName="simpleFoam_sweep",
    variations={"windSpeed": [3, 5, 8, 12]},
    jobName="wind_sweep"
)
```

## Reading mesh data

```python
# Get cell centers (runs foamJob postProcess internally)
mesh = of.getMesh(caseDirectory="/data/simulations/wind_study")

# Get mesh bounding box
extent = of.getMeshExtent(caseDirectory="/data/simulations/wind_study")

# Get list of computed time steps
times = of.getTimeList("wind_study_simpleFoam")
```

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

In the JSON, this produces a nested `downstream` structure:

```json
{
  "filters": {
    "boundary": {
      "filterType": "ExtractBlock", "write": false,
      "params": [["Selectors", ["/Root/boundary/inlet", "/Root/boundary/outlet"]]],
      "downstream": {
        "boundary_values": {
          "filterType": "CellCenters", "write": true,
          "params": [], "downstream": {}
        }
      }
    }
  }
}
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

## Exporting results

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

## Lagrangian dispersion (stochastic solver)

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
```

## Further reading

- **[Simulations API Reference](../../developer_guide/api/simulations.md)** — auto-generated API docs for all OpenFOAM classes
- **[Simulations Implementation](../../developer_guide/simulations.md)** — developer guide with architecture, solver hierarchy, VTK pipeline internals
- **[CLI Reference](../../cli/reference.md#hera-openfoam)** — full CLI command reference

For the full API, see the [API Reference](../../developer_guide/api/simulations.md).
