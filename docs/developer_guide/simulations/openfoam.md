# OpenFOAM toolkit

## Overall design

The OpenFOAM toolkit uses **composition over inheritance** — `OFToolkit` contains solver extensions, preprocessing objects, and post-processing pipelines as pluggable attributes rather than subclasses.

```
OFToolkit (hermesWorkflowToolkit)
├── OFObjectHome          — field definitions and creation
├── Analysis              — VTK post-processing pipeline
├── Presentation          — output formats (CSV, VTK)
├── Solver extensions (composition):
│   ├── buoyantReactingFoam_toolkitExtension   (attribute buoyantReactingFoam)
│   └── StochasticLagrangianSolver_toolkitExtension (attribute stochasticLagrangian)
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

The main class inherits from `hermesWorkflowToolkit` for workflow support and composes all other components. It is registered as `toolkitHome.SIMULATIONS_OPENFOAM` (`"OpenFOAM"`) and passes `toolkitName="OFworkflowToolkit"` to the base class.

`simpleFoam_toolkitExtension` (`eulerian/simpleFoam.py`) exists but is **not** composed into `OFToolkit`; instantiate it directly with the toolkit if needed. `OFLSMToolkit` is a separate toolkit (registry name `OF_LSM`), not an attribute.

| Attribute | Type | Purpose |
|-----------|------|---------|
| `OFObjectHome` | `OFObjectHome` | Field definitions catalog and field creation |
| `analysis` | `Analysis` | VTK pipeline and data extraction |
| `presentation` | `Presentation` | Output to CSV, VTK; parallel Lagrangian loading |
| `stochasticLagrangian` | `StochasticLagrangianSolver_toolkitExtension` | Lagrangian particle tracking |
| `buoyantReactingFoam` | `buoyantReactingFoam_toolkitExtension` | Compressible reactive solver |

**Flow type constants:**

| Constant | Value | Used for |
|----------|-------|----------|
| `FLOWTYPE_INCOMPRESSIBLE` | `"incompressible"` | simpleFoam, wind profile simulations |
| `FLOWTYPE_COMPRESSIBLE` | `"compressible"` | buoyantReactingFoam |
| `FLOWTYPE_DISPERSION` | `"dispersion"` | Lagrangian particle dispersion |

Other class constants (defined in `hera/simulations/openFoam/__init__.py` and re-exported on `OFToolkit`): `FIELDTYPE_SCALAR/VECTOR/TENSOR`, `CASETYPE_DECOMPOSED` (`"Decomposed Case"`), `CASETYPE_RECONSTRUCTED` (`"Reconstructed Case"`), `TIME_STEADYSTATE` / `TIME_DYNAMIC`, `DOCTYPE_OF_FLOWDISPERSION = "flowDispersion"`.

**Key methods:**

| Method | Description |
|--------|-------------|
| `createEmptyCase(caseDirectory, fieldList, flowType, additionalFieldsDescription=dict())` | Initialize OpenFOAM case directory structure |
| `writeEmptyField(fieldName, flowType, caseDirectory, timeOrLocation=0, readBoundaryFromCase=False, writeProcBoundary=False)` | Create empty field files with correct dimensions |
| `getMesh(caseDirectory, readParallel=True, time=0)` | Read cell centers (runs `foamJob ... postProcess -func writeCellCentres`) |
| `getMeshFromName(nameOrWorkflowFileOrJSONOrResource, ...)` | Same, resolving the case from a workflow name/document |
| `getMeshExtent(caseDirectory)` | Get bounding box from `constant/polyMesh/points` |
| `runOFSimulation(nameOrWorkflowFileOrJSONOrResource, scheduler="local", schedulerHost=None, schedulerPort=None, dispatch_id=None)` | Build and execute full simulation via workflow |
| `prepareSlurmWorkflowExecution(baseConfiguration, jsonVariations, slurmExecutionFileName="submit_all.sh", caseListFileName="cases.txt", ...)` | Generate Slurm batch submission scripts |
| `getTimeList(nameOrWorkflowFileOrJSONOrResourceorDirectory, singleProcessor=False, returnFirst=True)` | Extract computed time steps from case |
| `datasetToCaseFields(...)`, `netcdfToCaseFields(...)`, `datasetToSetFieldsDict(...)`, `netcdfToSetFieldsDict(...)`, `xarrayToSetFieldsDictDomain(...)` | Map xarray/netCDF data onto case fields or a `setFieldsDict` |
| `getVTKPipelineCacheDocuments(...)` / `getVTKPipelineCacheTable(...)` / `clearVTKPipelineCache(...)` | Query or clear cached VTK filter results (filter by `regularMesh`, `filterName`, `workflowName`, `groupName`) |

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
| `abstractWorkflow` | `hermes.workflow` | Node accessors `controlDict` (`ControlDict` node), `fvSolution`, `fvScheme`, `fileWriter`, `parameters` (`Parameters` node), `buildAllRun`, `defineNewBoundaryConditions`; `workflowGroup` / `workflowGroupID` |
| `workflow_Eulerian` | `abstractWorkflow` | + `blockMesh`, optional `snappyHexMesh` |
| `workflow_Lagrangian` | `abstractWorkflow` | Base for particle simulations |
| `workflow_simpleFoam` | `workflow_Eulerian` | Steady incompressible |
| `workflow_buoyantReactingFoam` | `workflow_Eulerian` | Compressible reacting |
| `workflow_scalarTransportFoam` | `workflow_Eulerian` | Scalar transport |
| `workflow_indoorFOAMBoussinesq` | `workflow_Eulerian` | Indoor Boussinesq approximation |
| `workflow_homogenousWindLogProfile` | `workflow_Eulerian` | Log-law wind profile |
| `workflow_StochasticLagrangianSolver` | `workflow_Lagrangian` | Stochastic Lagrangian dispersion; properties `dispersionName`, `originalFlowFieldName`, `dispersionFlowFieldName`, `dispersionDuration` |

**workflow_Eulerian key methods:**

| Method | Description |
|--------|-------------|
| `prepareICandBC(configurationFile)` | Set initial and boundary conditions from config |
| `set_blockMesh_blockBoundaries(minx, maxx, miny, maxy, minz, maxz, dx, dy, dz)` | Configure block mesh vertices and boundaries |
| `set_blockMesh_blockHeight(Z, dz)` | Adjust domain height |
| `ICandBCHandler_node(icnode)` | Process IC/BC node entries |
| `foam_snappyhexmesh_addobject(dataOrFile, objectFile)` | Add a geometry object to snappyHexMesh |

## Eulerian solver extensions (`eulerian/`)

Solver extensions handle solver-specific case setup. Only `buoyantReactingFoam_toolkitExtension` is composed into `OFToolkit` (as `buoyantReactingFoam`).

| Class | Solver | Properties |
|-------|--------|-----------|
| `absractEulerianSolver_toolkitExtension` | (base) | `solverName`, `incompressible`, `flowType`, `analysis`, `presentation` |
| `simpleFoam_toolkitExtension` | `simpleFOAM` | `incompressible=True`; adds `buildSpecializedField(workflow, turbulenceType)` |
| `buoyantReactingFoam_toolkitExtension` | `buoyantReactingFoam` | `incompressible=False` |

**Key methods on the abstract solver:**

| Method | Description |
|--------|-------------|
| `blockMesh_setBoundFromBounds(eulerianWF, minx, maxx, miny, maxy, minz, maxz, dx, dy, dz)` | Configure mesh boundaries from coordinates |
| `blockMesh_setBoundFromFile(eulerianWorkFlow, fileName, dx, dy, dz)` | Configure mesh from OBJ geometry file |
| `blockMesh_setDomainHeight(eulerianWorkFlow, Z, dz)` | Adjust vertical extent |
| `writeFieldInCase(fieldName, caseDirectory, components, xarrayData, boundaryConditions=None, data=None)` | Write a field file to the case |
| `IC_getHydrostaticPressure(caseDirectory, fieldName="p", groundPressure=101000)` | Generate hydrostatic pressure initial condition |

## Lagrangian solver extensions (`lagrangian/`)

**Hierarchy:**

| Class | Purpose |
|-------|---------|
| `absractStochasticLagrangianSolver_toolkitExtension` (`abstractLagrangianSolver.py`) | Base for Lagrangian solvers — dispersion flow fields, source generation, case results; holds an `analysis` object (concentration methods) |
| `StochasticLagrangianSolver_toolkitExtension` (`StochasticLagrangianSolver.py`) | Concrete implementation; overrides `createDispersionFlowField` (default dispersion fields such as `ustar`, `distanceFromWalls`) |
| `OFLSMToolkit` (`LSM/toolkit.py`) | Separate OpenFOAM + LSM coupling toolkit — particle data loading, concentration calculation |

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
| `makeSource_Point/Circle/Sphere/Cylinder/Rectangle/Cube(x, y, z, nParticles, ...)` | Generate particle source geometries (`sourcesTypeList` lists them) |
| `writeParticlePositionFile(x, y, z, nParticles, dispersionName, type="Point", ...)` | Write particle initial positions |

### OFLSMToolkit (`lagrangian/LSM/toolkit.py`)

Inherits from `abstractToolkit` (not from OFToolkit) and provides OpenFOAM + LSM coupling. Registry name `OF_LSM` (no `toolkitHome` constant); constructor `OFLSMToolkit(projectName, casePath=None, cloudName="kinematicCloud", filesDirectory=None, parallelCase=False, connectionName=None)`.

!!! warning "Known issue in master"
    `OFLSMToolkit.__init__` calls `toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_TOPOGRAPHY, ...)` (`lagrangian/LSM/toolkit.py:97`), but `ToolkitHome` defines no `GIS_TOPOGRAPHY` constant (only `GIS_RASTER_TOPOGRAPHY` / `GIS_VECTOR_TOPOGRAPHY`), so instantiation raises `AttributeError`.

| Property | Description |
|----------|-------------|
| `casePath` | OpenFOAM case directory |
| `cloudName` | Lagrangian particle cloud name |
| `parallelCase` | Decomposed (parallel) case flag |
| `sourcesFactory` | `sourcesFactoryTool` for source geometries |
| `topography` | GIS topography toolkit reference |

| Method | Description |
|--------|-------------|
| `loadData(times=None, saveMode=TOOLKIT_SAVEMODE_NOSAVE, withVelocity=False, withReleaseTimes=False, withMass=False, releaseTime=0, cloudName="kinematicCloud", casePath=None, **kwargs)` | Extract particle positions/velocities/masses as dask DataFrame |
| `makeSource(x, y, z, nParticles, type="Point", fileName="kinematicCloudPositions", **kwargs)` | Write a particle source file via `sourcesFactory` |
| `makeDistanceFromGround(times, ground="ground", fileName="cellHeights", resolution=10, ...)` | Compute ground distance field |
| `makeUstar(times, fileName="ustar", ground="ground", heightLimits=[1, 2], ...)` | Compute friction velocity field (decomposed into 6 helpers) |
| `createRootCaseMeshLink(rootCase)` | Link the root case mesh into the case |
| `to_paraview_CSV(data, outputdirectory, filename, timeFactor=1)` | Export particles as timestep CSV |
| `toUnstructuredVTK(data, outputdirectory, filename, timeNameOutput=True)` | Export as VTU particle files |
| `toStructuredVTK(data, outputdirectory, filename, extents, dxdydz, timeNameOutput=True)` | Export as structured VTK grid |

**makeUstar helper decomposition:**

| Helper | Purpose |
|--------|---------|
| `_loadOrComputeCellData()` | Cached cell height data lookup/computation |
| `_parseOFFieldHeader()` | Find internalField start + cell count |
| `_parseBoundaryStart()` | Find boundaryField line number |
| `_extractVelocityField()` | Parse OF vector, compute magnitude |
| `_interpolateNearGroundValues()` | Regularize to 2D grid, interpolate with (x,y) cache |
| `_writeScalarField()` | Assemble header + values + boundary → write file |

**Analysis class** on OFLSMToolkit (`lagrangian/LSM/toolkit.py`; the Lagrangian solver extension's `analysis` class in `abstractLagrangianSolver.py` has methods of the same names, defaulting to `x/y/z` instead of `globalX/globalY/globalZ`):

| Method | Description |
|--------|-------------|
| `calcConcentrationPointWise(data, dxdydz, xfield="globalX", ...)` | Cell-binned concentration from particle data |
| `calcDocumentConcentrationPointWise(dataDocument, dxdydz, ...)` | Same with DB caching |
| `calcConcentrationTimeStepFullMesh(timeData, extents, dxdydz, ...)` | Full domain regular grid concentration |
| `calcConcentrationFieldFullMesh(dataDocument, extents, dxdydz, ...)` | Full pipeline: load particles → bin → write netCDF (the helper decomposition above applies to the solver-extension version, not this one) |
| `getConcentrationField(dataDocument, returnFirst=True, **metadata)` | Load a cached concentration field |
| `getMassFromLog(logFile, solver="StochasticLagrangianSolver")` | Parse released/escaped mass from a solver log |

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
| `getEmptyField(fieldName, flowType, noOfProc=None, addParallelProc=False)` | Create empty field with correct dimensions |
| `getEmptyFieldFromCase(fieldName, flowType, caseDirectory, internalValue=None, readParallel=True)` | Create field pre-loaded with boundaries from an existing case |
| `readFieldFromCase(fieldName, flowType, caseDirectory, timeStep=0, readParallel=True)` | Read complete field data |
| `readFieldAsDataFrame(fieldName, caseDirectory, times=0, readParallel=True, filterInternalPatches=False)` | Read field across multiple time steps as DataFrame |
| `addFieldDefinitions(fieldName, dimensions, fieldType, fieldComputation=FIELDCOMPUTATION_EULERIAN, compressible_dimensions=None, overwrite=False)` | Register custom field definitions |

## VTK post-processing pipeline (`postProcess/`)

The VTK pipeline provides a ParaView-integrated post-processing system with DB caching.

**Pipeline architecture:**

| Class | Purpose |
|-------|---------|
| `VTKPipeLine` | Pipeline container — creates/manages filters (`addFilter(name, filterType, write=True, params=None)`), exports to JSON; obtained via `OFToolkit.analysis.getVTKPipeline()` |
| `registeredVTKPipeLine` | Pipeline bound to a specific case — executes filters, caches results |
| `VTKFilter` | Base filter node — tree structure with `downstream` children |
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
VTKPipeLine.registerPipeline(nameOrWorkflowFileOrJSONOrResource, serverName=None, caseType=CASETYPE_DECOMPOSED)
    → registeredVTKPipeLine
        → getData(regularMesh, filterName=None, timeList=None, latestTime=False, fieldNames=None, overwrite=False)
          (or getRegularData(...) / getNonRegularData(...))
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

**`writeCase` decomposition** (`paraviewOpenFOAM.writeCase` in `pvOpenFOAMBase.py`; was 127 lines, CC 20):

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

There is no `saveTemplate` / `loadTemplate` API in master. Templates are Hermes flow/node templates stored as toolkit datasources (`desc.component` = `"Flow"` or `"Node"`, flows also tagged with `solver`) and read through `hermesWorkflowToolkit`:

- `listHermesSolverTemplates(solverName)` / `getHermesFlowTemplate(hermesFlowName)` — list and fetch flow templates
- `listHermesNodesTemplates()` / `getHermesNodeTemplate(hermesNodeName)` — list and fetch node templates
- CLI: `hera-openFoam <solver> templates list` and `hera-openFoam <solver> templates get`

`OFToolkit.template_add(name, objFile, workflowObj=None)` exists but is an unimplemented stub (`pass`, `toolkit.py:502-525`).

## CLI

The `hera-openFoam` CLI (`hera/bin/hera-openFoam`, handlers in `openFoam/CLI.py`) exposes per-solver subcommands (`simpleFoam`, `pimpleFoam`, `oguraPhillipShallowConvectionFOAM`, `buoyantReactingFoam`, `scalarTransportFoam`, `homogenousWindLogProfile`, `stochasticLagrangian`) with actions such as `templates`, `executeWorkflow`, `list`, `delete`, `blockMesh`, `IC`, `case`, `postProcess`. See the [CLI reference](../../cli/reference.md#hera-openfoam).
