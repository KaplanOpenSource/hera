# OpenFOAM toolkit

## Overall design

The OpenFOAM toolkit uses **composition over inheritance** — `OFToolkit` contains solver extensions, preprocessing objects, and post-processing pipelines as pluggable attributes rather than subclasses.

```
OFToolkit (hermesWorkflowToolkit)
├── OFObjectHome          — field definitions and creation
├── Analysis              — VTK post-processing pipeline
├── Presentation          — output formats (CSV, VTK)
├── Solver extensions (composition):
│   ├── simpleFoam_toolkitExtension
│   ├── buoyantReactingFoam_toolkitExtension
│   └── StochasticLagrangianSolver_toolkitExtension
│       └── OFLSMToolkit  — LSM coupling
└── Workflows (Hermes-integrated):
    ├── workflow_Eulerian
    │   ├── workflow_simpleFoam
    │   ├── workflow_buoyantReactingFoam
    │   ├── workflow_scalarTransportFoam
    │   ├── workflow_indoorFOAMBoussinesq
    │   └── workflow_homogenousWindLogProfile
    └── workflow_Lagrangian
        └── workflow_StochasticLagrangianSolver
```

## OFToolkit (`toolkit.py`)

The main class inherits from `hermesWorkflowToolkit` for workflow support and composes all other components:

| Attribute | Type | Purpose |
|-----------|------|---------|
| `OFObjectHome` | `OFObjectHome` | Field definitions catalog and field creation |
| `analysis` | `Analysis` | VTK pipeline and data extraction |
| `presentation` | `Presentation` | Output to CSV, VTK, parquet |
| `stochasticLagrangian` | `StochasticLagrangianSolver_toolkitExtension` | Lagrangian particle tracking |
| `buoyantReactingFoam` | `buoyantReactingFoam_toolkitExtension` | Compressible reactive solver |

**Flow type constants:**

| Constant | Value | Used for |
|----------|-------|----------|
| `FLOWTYPE_INCOMPRESSIBLE` | — | simpleFoam, wind profile simulations |
| `FLOWTYPE_COMPRESSIBLE` | — | buoyantReactingFoam |
| `FLOWTYPE_DISPERSION` | — | Lagrangian particle dispersion |

**Key methods:**

| Method | Description |
|--------|-------------|
| `createEmptyCase(caseDirectory, fieldList, flowType)` | Initialize OpenFOAM case directory structure |
| `writeEmptyField(fieldName, flowType, caseDirectory, ...)` | Create empty field files with correct dimensions |
| `getMesh(caseDirectory)` | Read cell centers (runs `foamJob postProcess`) |
| `getMeshExtent(caseDirectory)` | Get bounding box from `points` file |
| `runOFSimulation(nameOrWorkflowFileOrJSONOrResource)` | Build and execute full simulation via workflow |
| `prepareSlurmWorkflowExecution(...)` | Generate Slurm batch submission scripts |
| `getTimeList(...)` | Extract computed time steps from case |
| `getVTKPipelineCacheTable(...)` | Query cached VTK filter results |

**Nested classes:**

| Class | Methods | Purpose |
|-------|---------|---------|
| `Analysis` | `getVTKPipeline()`, `executeVTKFiltersAndLoad()`, `getFiltersDocuments()` | Create and run VTK filter pipelines |
| `Presentation` | `to_paraview_CSV()`, `toUnstructuredVTK()`, `toStructuredVTK()`, `loadLagrangianDataParallel()` | Export results in various formats |

## Workflow system (`OFWorkflow.py`)

Workflows define the steps to build and run an OpenFOAM case. They integrate with the Hermes framework.

**Workflow hierarchy:**

| Class | Inherits from | Nodes |
|-------|--------------|-------|
| `abstractWorkflow` | `hermes.workflow` | `controlDict`, `fvSolution`, `fvSchemes`, `fileWriter`, `Parameters` |
| `workflow_Eulerian` | `abstractWorkflow` | + `blockMesh`, optional `snappyHexMesh` |
| `workflow_Lagrangian` | `abstractWorkflow` | Base for particle simulations |
| `workflow_simpleFoam` | `workflow_Eulerian` | Steady incompressible |
| `workflow_buoyantReactingFoam` | `workflow_Eulerian` | Compressible reacting |
| `workflow_scalarTransportFoam` | `workflow_Eulerian` | Scalar transport |
| `workflow_indoorFOAMBoussinesq` | `workflow_Eulerian` | Indoor Boussinesq approximation |
| `workflow_homogenousWindLogProfile` | `workflow_Eulerian` | Log-law wind profile |
| `workflow_StochasticLagrangianSolver` | `workflow_Lagrangian` | Stochastic Lagrangian dispersion |

**workflow_Eulerian key methods:**

| Method | Description |
|--------|-------------|
| `prepareICandBC(configurationFile)` | Set initial and boundary conditions from config |
| `set_blockMesh_blockBoundaries(...)` | Configure block mesh vertices and boundaries |
| `set_blockMesh_blockHeight(...)` | Adjust domain height |
| `ICandBCHandler_node(icnode)` | Process IC/BC node entries |

## Eulerian solver extensions (`eulerian/`)

Solver extensions are composed into `OFToolkit` as attributes. They handle solver-specific case setup.

| Class | Solver | Properties |
|-------|--------|-----------|
| `absractEulerianSolver_toolkitExtension` | (base) | `solverName`, `incompressible`, `flowType`, `analysis`, `presentation` |
| `simpleFoam_toolkitExtension` | `simpleFOAM` | `incompressible=True` |
| `buoyantReactingFoam_toolkitExtension` | `buoyantReactingFoam` | `incompressible=False` |

**Key methods on the abstract solver:**

| Method | Description |
|--------|-------------|
| `blockMesh_setBoundFromBounds(wf, minx, maxx, ...)` | Configure mesh boundaries from coordinates |
| `blockMesh_setBoundFromFile(wf, fileName, dx, dy, dz)` | Configure mesh from OBJ geometry file |
| `blockMesh_setDomainHeight(wf, Z, dz)` | Adjust vertical extent |
| `writeFieldInCase(fieldName, caseDirectory, ...)` | Write a field file to the case |
| `IC_getHydrostaticPressure(caseDirectory)` | Generate hydrostatic pressure initial condition |

## Lagrangian solver extensions (`lagrangian/`)

**Hierarchy:**

| Class | Purpose |
|-------|---------|
| `absractStochasticLagrangianSolver_toolkitExtension` | Base for Lagrangian solvers — dispersion flow fields, source generation |
| `StochasticLagrangianSolver_toolkitExtension` | Concrete implementation with default fields (ustar, distanceFromWalls) |
| `OFLSMToolkit` | OpenFOAM + LSM coupling — particle data loading, concentration calculation |

**Abstract Lagrangian solver key methods:**

| Method | Description |
|--------|-------------|
| `createDispersionFlowField(flowName, flowData, ...)` | Create flow field for dispersion (links to original flow, maps time steps, writes ustar/Hmix) |
| `createDispersionCaseDirectory(hermes_dispersionWorkflow, ...)` | Set up dispersion case with processor links and DB sync |
| `makeSource_Point/Circle/Sphere/Cylinder/Rectangle/Cube(...)` | Generate particle source geometries |
| `writeParticlePositionFile(x, y, z, nParticles, ...)` | Write particle initial positions |

**OFLSMToolkit (`lagrangian/LSM/toolkit.py`):**

Inherits from `abstractToolkit` (not from OFToolkit) and provides OpenFOAM + LSM coupling:

| Property | Description |
|----------|-------------|
| `casePath` | OpenFOAM case directory |
| `cloudName` | Lagrangian particle cloud name |
| `parallelCase` | Decomposed (parallel) case flag |
| `sourcesFactory` | `sourcesFactoryTool` for source geometries |
| `topography` | GIS topography toolkit reference |

| Method | Description |
|--------|-------------|
| `loadData(times, saveMode, withVelocity, ...)` | Extract particle positions/velocities/masses as dask DataFrame |
| `makeDistanceFromGround(times, ground, ...)` | Compute ground distance field |
| `makeUstar(times, ...)` | Compute friction velocity field |
| `to_paraview_CSV(data, ...)` | Export particles as timestep CSV |
| `toUnstructuredVTK(data, ...)` | Export as VTU particle files |
| `toStructuredVTK(data, ...)` | Export as structured VTK grid |

**Analysis inner class** on OFLSMToolkit:

| Method | Description |
|--------|-------------|
| `calcConcentrationPointWise(data, dxdydz, ...)` | Cell-binned concentration from particle data |
| `calcDocumentConcentrationPointWise(dataDocument, ...)` | Same with DB caching |
| `calcConcentrationTimeStepFullMesh(...)` | Full domain regular grid concentration |

## Preprocessing objects (`preprocessOFObjects/`)

The preprocessing system provides an abstraction for OpenFOAM field files:

| Class | Purpose |
|-------|---------|
| `OFObject` | Abstract base — field name, dimensions, regions (internalField/boundaryField), write to disk |
| `OFField(OFObject)` | Concrete field — initialize empty, set uniform values, add boundaries, read from case, convert to DataFrame |
| `OFObjectHome` | Field catalog — predefined fields (U, p, k, epsilon, T, ...) with correct dimensions for incompressible/compressible flows |

**OFObjectHome predefined fields:** U, p, p_rgh, epsilon, omega, alphat, nut, k, T, Tbackground, cellCenters, Hmix, ustar, distanceFromWalls

**Key OFObjectHome methods:**

| Method | Description |
|--------|-------------|
| `getEmptyField(fieldName, flowType, noOfProc)` | Create empty field with correct dimensions |
| `getEmptyFieldFromCase(fieldName, flowType, caseDirectory, ...)` | Create field pre-loaded with boundaries from an existing case |
| `readFieldFromCase(fieldName, flowType, caseDirectory, timeStep)` | Read complete field data |
| `readFieldAsDataFrame(fieldName, caseDirectory, times, ...)` | Read field across multiple time steps as DataFrame |
| `addFieldDefinitions(fieldName, dimensions, fieldType, ...)` | Register custom field definitions |

## VTK post-processing pipeline (`postProcess/`)

The VTK pipeline provides a ParaView-integrated post-processing system with DB caching.

**Pipeline architecture:**

| Class | Purpose |
|-------|---------|
| `VTKPipeLine` | Pipeline container — creates/manages filters, exports to JSON |
| `registeredVTKPipeLine` | Pipeline bound to a specific case — executes filters, caches results |
| `VTKFilter` | Base filter node — tree structure with `downstream` children |
| `VTKPipelineExecutionContext` | Reader and pipeline configuration for paraview execution |
| `paraviewOpenFOAM` (`pvOpenFOAMBase.py`) | ParaView backend — reads OpenFOAM cases, writes results |

**Available filters:**

| Filter class | Type constant | Purpose |
|-------------|--------------|---------|
| `vtkFilter_Slice` | `FILTER_SLICE` | Planar slice through the domain |
| `vtkFilter_PlotOverLine` | `FILTER_PLOTOVERLINE` | Sample along a line |
| `vtkFilter_CellCenters` | `FILTER_CELLCENTERS` | Extract cell center values |
| `vtkFilter_ExtractBlock` | `FILTER_EXTRACTBLOCK` | Extract specific mesh regions/patches |
| `vtkFilter_DescriptiveStatistics` | — | Compute summary statistics |
| `vtkFilter_IntegrateVariables` | `FILTER_INTEGRATEVARIABLES` | Integrate fields over domain |

**Data flow:**

```
VTKPipeLine.registerPipeline(case)
    → registeredVTKPipeLine
        → getData(filterName, timeList, ...)
            → checks DB cache
            → if miss: runs ParaView pipeline
            → saves to parquet (non-regular) or ZARR (regular mesh)
            → stores metadata document in DB
            → returns dict of filter results
```

**Filter tree:** Filters form a tree — each filter can have downstream children. The pipeline traverses the tree depth-first, executing each filter and optionally writing results to disk.

## Template system

Templates are saved OpenFOAM case directories stored as data sources. The toolkit provides:
- `saveTemplate(templateName, caseDirectory)` — save a case as a template
- `loadTemplate(templateName, toDirectory, caseName)` — instantiate a template
