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

### Dispersion flow field pipeline (`createDispersionFlowField`)

The most complex operation in the Lagrangian toolkit — creates a dispersion case directory from an existing flow simulation. The method is decomposed into a readable orchestrator calling private helpers:

| Step | Helper | Purpose |
|------|--------|---------|
| 1 | `_locateOriginalFlowCase()` | Find flow in DB or filesystem |
| 2 | `_detectParallelAndTimesteps()` | Detect parallel case + discover time directories |
| 3 | `_buildTimeMapping()` | Build (original→dispersion) time mapping for steady-state or dynamic flows |
| 4 | `_checkAndRegisterDispersionInDB()` | Query DB for existing record, handle overwrite |
| 5 | `_prepareDispersionDirectory()` | Copy/link constant/system/0 + per-timestep data + write dispersion fields |
| 6 | `_finalizeDispersionInDB()` | Register or update DB document |

**Time mapping logic:**

- **Steady-state**: freezes the flow at one timestep — maps to `t=0` and `t=dispersionDuration`
- **Dynamic**: shifts each flow timestep so dispersion starts at `t=0`. If the flow ends before the dispersion duration, the last timestep is repeated.

### Case results and caching

`getCaseResults` and `getCaseConcentrationsEulerian` share the same cache/load/compute pattern via 7 shared helpers:

| Helper | Purpose |
|--------|---------|
| `_resolveCaseDescriptorName()` | Extract name from str or MetadataFrame |
| `_checkCache()` | Lookup cache, handle overwrite, clean stale entries |
| `_locateCaseDirectory()` | DB or filesystem case path resolution |
| `_discoverTimesteps()` | Find numeric time dirs, filter system dirs |
| `_loadCaseDataViaDask()` | Parallel/serial detection + Dask distributed map |
| `_saveToCacheParquet()` | Parquet save + DB registration (Lagrangian) |
| `_saveToCacheNetCDF()` | xarray groupby+sum + netCDF save (Eulerian concentrations) |

### Dispersion case directory (`createDispersionCaseDirectory`)

Validates DB consistency, handles name/parameter conflicts, and creates the case directory:

| Helper | Purpose |
|--------|---------|
| `_resolveDispersionFlowField()` | DB lookup for flow field, name canonicalization |
| `_checkDBConsistency()` | Detect name conflicts + parameter duplicates |
| `_resolveDirectoryConflict()` | Handle existing directory with rewrite flag |
| `_syncDispersionToDB()` | Add new or update existing DB record |

### Concentration field computation (`calcConcentrationFieldFullMesh`)

Bins Lagrangian particles onto a Cartesian mesh, stores as partitioned netCDF:

| Helper | Purpose |
|--------|---------|
| `_removeConcentrationCache()` | Clean old cache files + DB record |
| `_createConcentrationCacheDoc()` | Counter-based path + DB registration |
| `_processConcentrationPartition()` | Per-Dask-partition timestep binning |
| `_padRemainingTimesteps()` | Zero-fill to full dispersion duration |
| `_writeConcentrationPartition()` | Concat + rename axes + write netCDF |

### Source generation

| Method | Description |
|--------|-------------|
| `makeSource_Point/Circle/Sphere/Cylinder/Rectangle/Cube(...)` | Generate particle source geometries |
| `writeParticlePositionFile(x, y, z, nParticles, ...)` | Write particle initial positions |

### OFLSMToolkit (`lagrangian/LSM/toolkit.py`)

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
| `makeUstar(times, ...)` | Compute friction velocity field (decomposed into 6 helpers) |
| `to_paraview_CSV(data, ...)` | Export particles as timestep CSV |
| `toUnstructuredVTK(data, ...)` | Export as VTU particle files |
| `toStructuredVTK(data, ...)` | Export as structured VTK grid |

**makeUstar helper decomposition:**

| Helper | Purpose |
|--------|---------|
| `_loadOrComputeCellData()` | Cached cell height data lookup/computation |
| `_parseOFFieldHeader()` | Find internalField start + cell count |
| `_parseBoundaryStart()` | Find boundaryField line number |
| `_extractVelocityField()` | Parse OF vector, compute magnitude |
| `_interpolateNearGroundValues()` | Regularize to 2D grid, interpolate with (x,y) cache |
| `_writeScalarField()` | Assemble header + values + boundary → write file |

**Analysis inner class** on OFLSMToolkit:

| Method | Description |
|--------|-------------|
| `calcConcentrationPointWise(data, dxdydz, ...)` | Cell-binned concentration from particle data |
| `calcDocumentConcentrationPointWise(dataDocument, ...)` | Same with DB caching |
| `calcConcentrationTimeStepFullMesh(...)` | Full domain regular grid concentration |
| `calcConcentrationFieldFullMesh(...)` | Full pipeline: load particles → bin → pad → write netCDF (decomposed into helpers) |

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

**`getData` decomposition** (was 181 lines, CC 15):

| Step | Helper | Purpose |
|------|--------|---------|
| 1 | `_resolveRequestedFilters()` | Coerce None/str/list → concrete filter name list |
| 2 | `_parseTimeList()` | Handle None (all times), `"start:end"` range, or explicit list |
| 3 | `_filterCachedTimesteps()` | Subtract already-cached timesteps for incremental computation |
| 4 | `_buildAndExecuteParaViewPipeline()` | Create reader → build filter tree → execute → write results |
| 5 | `_updateCacheDB()` | Merge new timesteps into existing DB records |

**`writeCase` decomposition** (was 127 lines, CC 20):

| Step | Helper | Purpose |
|------|--------|---------|
| 1 | `_removeOldOutputs()` | Clean existing files before overwrite |
| 2 | `_ensureOutputDirs()` | Create output directories |
| 3 | `_resolveTimeList()` | Determine timesteps from reader or caller |
| 4 | `_writeTimeStepBlocks()` | Stream timesteps in fixed-size blocks to temp files |
| 5 | `_collectTmpFiles()` | Glob temp block files per filter |
| 6 | `_mergeToFinalOutput()` | Dispatch to `_mergeZarr` or `_mergeParquet` |
| 7 | `_atomicReplace()` | `.final` staging → rename |
| 8 | `_cleanupTmpFiles()` | Remove temp blocks after merge |

The merge step supports both regular meshes (xarray → zarr with time-concatenation) and unstructured meshes (dask → parquet with repartitioning to ~100 MB chunks). The atomic `.final` → rename pattern ensures no partial writes are visible.

## Template system

Templates are saved OpenFOAM case directories stored as data sources. The toolkit provides:
- `saveTemplate(templateName, caseDirectory)` — save a case as a template
- `loadTemplate(templateName, toDirectory, caseName)` — instantiate a template
