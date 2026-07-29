# Dead Code / Unused Code Report — hera

Consolidated report from a comprehensive scan of the entire `hera/` package (datalayer, measurements/GIS, measurements/meteorology, measurements/experiment, simulations, riskassessment, utils, toolkit.py, bin/CLI). Every finding was verified with grep across the whole repo (including tests, docs, JSON) before being flagged as "dead" — to rule out dynamic usage (pydoc.locate, getattr, factories, `__all__`).

**General recommendation:** before actually deleting anything — run the existing test suite, and confirm with the code owners (especially in simulations/riskassessment, which expose a lot of public API) that there's no external usage (notebooks/scripts outside the repo).

---

## 🔴 Urgent — broken code (not just dead, will crash if run)

These are the most severe findings: code that looks "in use" (real call paths reach it) but contains a bug that will cause a crash. Some of these effectively block entire toolkits.

1. **`hera/simulations/gaussian/gasCloud.py`** — the entire gaussian toolkit is broken in practice:
   - Line 393-397: `instantaneousReleaseGasCloud` uses default arguments `dz=1*m, dt=1*min` — `m`/`min` are undefined names (should be `ureg.m`/`ureg.min`). Default argument expressions are evaluated **at import time**, so the `NameError` happens as soon as the module is loaded. Verified by actually executing (`ast.parse` + running with stubs).
   - As a result, `hera/simulations/gaussian/toolkit.py` (which does `from hera.simulations.gaussian.gasCloud import abstractGasCloud`) **cannot be imported** — meaning the `gaussianToolkit` registered in `toolkitHome` is currently broken. Likely masked because `pydoc.locate` silently swallows import errors.
   - Lines 325-390: an older, duplicate definition of `instantaneousReleaseGasCloud` — immediately shadowed by the second definition (fully dead, unreachable).
   - **Fix:** replace `m`→`ureg.m`, `min`→`ureg.min` (line 397), delete the first duplicate definition (325-390).

2. **`hera/simulations/gaussian/source.py:12`** — an actual syntax error (`self.Q = #need to see...` with no value) — `SyntaxError`, the module cannot be imported at all. Nothing currently imports it, so it doesn't break anything live, but it's 100% invalid code.

3. **`hera/simulations/utils/interpolation/` (whole directory)** — `interpolations.py:137` contains a `TabError` (mixed tabs/spaces) — cannot be imported. This is a broken duplicate of `hera/simulations/utils/interpolations.py` (which IS actually used by `windProfile/toolkit.py`). **Full deletion of the duplicate directory is recommended** (also flagged by a prior report, `CODE_REVIEW_REPORT.md:188`).

4. **`hera/riskassessment/agents/effects/thresholdGeoDataFrame.py:78,90`** — `isinstance(x, collections.Iterable)` — `collections.Iterable` was removed as of Python 3.10 (moved to `collections.abc.Iterable`). Any call to `project()` with a list of angles will crash. **Easy fix**: `collections.abc.Iterable`.

5. **`hera/riskassessment/protectionpolicy/ProtectionPolicy.py:128`** — `os.JSONpath.exists(...)` — a typo (`os` has no `JSONpath`), should be `os.path.exists`. Currently not triggered because `addActions` is always called with a `dict`, but the documented API (string/JSON file path) is completely broken.

6. **`hera/riskassessment/riskToolkit.py:274` (`analysis.getRiskAreas`)** — calls `calculateRaw(...)`, which doesn't exist anywhere on `Injury`/`InjuryLevel`. Will always crash if called.

7. **`hera/riskassessment/analysis/riskAreas.py:21-24` (`getRiskAreaAlgorithm`)** — uses `pydoc.locate("pyriskassessment...")` — an old package name from before the rename to `hera` — will always return `None` and raise `ValueError`.

8. **`hera/simulations/openFoam/eulerian/abstractEulerianSolver.py`** — a whole cluster of methods that will hit `NameError`/`AttributeError` the moment they're called (`flowType` line 30 references the nonexistent `SIMULATIONTYPE_COMPRESSIBLE`; `blockMesh_setBoundFromFile`/`blockMesh_setDomainHeight` lines 88/111/116 reference nonexistent parameters; they call `set_blockMesh_boundaries`, which doesn't exist anywhere in `OFWorkflow.py`). "Public" API that has apparently never actually run.

9. **`hera/simulations/openFoam/OFWorkflow.py:74-83`** — `@workflowGroup.setter` (should be `@workflowGroupID.setter`) — accidentally overwrites the original `workflowGroupID` property (rendering it unreachable).

10. **`hera/simulations/openFoam/postProcess/pvOpenFOAMBase.py`** — `write_parquet` is defined twice, identically (lines 631-634 and 637-640) — the first definition is dead-on-arrival. Both bodies call a `writeNonRegularCase` function that doesn't exist (should be `self.writeCase`) — guaranteed crash. Also `write_netcdf`, `to_pandas`, `to_xarray`, `to_dataFrame`, `to_dataArray` (all `@deprecated`) call nonexistent methods/parameters — completely broken.

11. **`hera/simulations/machineLearningDeepLearning/torch/modelContainer.py:211-222, 266-277` (`fit`/`only_validate`)** — `ckpt_path = None`, then `elif os.path.exists(ckpt_path):` is checked while `ckpt_path` is still `None` → `TypeError`. This is **live** code (not just "dead") with a real bug — worth fixing, not deleting.

---

## 1️⃣ Unused functions/classes/methods/variables (Unused Definitions)

### datalayer / utils / toolkit.py / bin
| file:line | what | why dead |
|---|---|---|
| `datalayer/collection.py:171` | `AbstractCollection.addDocumentFromJSON` | grep — single hit (the definition) |
| `utils/query.py:1` | `andClause(excludeFields=[], **kwargs)` | 0 call sites; also a bug: mutable default argument |
| `utils/zipUtils.py:9` | `add_directory_to_zip` | 0 call sites; even `zip_items` in the same file doesn't use it — duplicated logic |
| `utils/unitHandler.py:304` | `extractUnumUnitsFromPint` | 0 call sites; fully duplicates the body of `pintToUnum` |
| `utils/unitHandler.py:357` | `unumToBaseUnits` | 0 call sites, not in `__all__` |
| `toolkit.py:128` | `abstractToolkit.classLoggerName` (property) | 0 call sites |
| `bin/hera-OF-postProcess.old` | whole file | intentionally filtered by `setup.py` (`.old` suffix) — dead code sitting in the tree |
| `bin/jupyter-lab-server` | CLI file | **broken registration**: the glob pattern in `setup.py:19` (`hera-*`) doesn't match it since it lacks the `hera-` prefix → never installed, despite being documented |

### GIS (raster + vector)
- `raster/hill2stl.py` — **entire file** is dead: executes `generate_solid_stl(...)` at import time (dangerous side-effect — writes `test1.stl` to disk if anyone imports it).
- `raster/landcover.py:561-585` — `roughnesslength2sandgrainroughness` is defined **twice** (also at 717) — the first definition is unreachable.
- `raster/landcover.py:488-522` — `_handleType1` — 0 call sites at all.
- `raster/tiles.py` — `tileScaleAtLatLonZoom`, `listImages`, `setDefaultTileServer` — 0 call sites, no `test_tiles.py`.
- `vector/buildings/toolkit.py:202-238,241-292` — `get_buildings_height`, `filter_buildings_in_area` (static) — 0 call sites.
- `vector/toolkit.py:9,37-57` — `TOOLKIT_VECTOR_REGIONNAME`, `geopandasToGeoJson` — 0 call sites.
- `vector/topography.py:99-113,158-217` — `geoPandasToSTL`, `toDEM` — 0 call sites.
- `vector/demography.py:113-119` — `projectPolygonOnPopulation` — deprecated shim (`DeprecationWarning`), 0 call sites — its own docstring says to remove it.

### meteorology (highfreqdata)
- `toolkit.py:26-27` — `DOCTYPE_STATIONS`, `DOCTYPE_MEASUREMENTS` — 0 call sites.
- `analysis/meandatacalculator.py:197-198` — `_UV_to_SpdDir` — 0 call sites.
- `analysis/turbulencestatistics.py:1352-1425` — entire class `SinglePointStatisticsSpark` — 0 imports/usage anywhere.
- `analysis/turbulencestatistics.py:1428-1549` — `InMemoryRawData.append/read_hdf/to_hdf` — 0 call sites (the subclass `InMemoryAvgData` is used, but not through these methods).
- `parsers/CampbellBinary.py` — `instrument`, `station`, `firstTime`, `lastTime` (properties), `getRecordByTime` — all 0 call sites; live code reads `headers[0].split(...)` directly instead of going through the properties.
- `analysis/abstractcalculator.py:62-65` — `JoinMethod` property — 0 call sites.

### experiment
- `dataEngine.py` — `pandasDataEngineDB`, `daskDataEngineDB` (entire classes) — only ever constructed through a factory that always receives the default `PARQUETHERA`; both are also internally broken (`self.experimentObj` is never set, `dask_mongo` is never imported).
- `analysis.py:29-52,232-254` — `getDeviceLocations`, `getDeviceTypePlannedMessageCount` (calls the nonexistent `getOptimalFrequencyHz`) — no real test coverage.
- `presentation.py:219-248,449-497,648-700` — `plotMap`, `plotDevices`, `generateLatexTable` — use nonexistent parameters/attributes (`self.trialSet` should be `self.datalayer.trialSet`, etc.) — broken and also never called.
- `parsers.py:356-372` — `Parser_TOA5` — 0 call sites, empty body (`pass`).
- `experiment.py:221-225,382-399` — `experimentDataType` (always `None`), `_initAnalysisAndPresentation`, `trialsOfDefaultTrialSet` — 0 call sites.

### riskassessment
- Three factory-eligible classes that are never dispatched via JSON/config: **`CalculatorHaber`**, **`InjuryExponential`**, **`InjuryLevelExponential`** (verified: 0 occurrences of the strings `"Haber"`/`"Exponential"` in any JSON/config in the repo).
- `Agents.py:60-80` — `Agent.fullDescription`, `Agent.effectproperties` — 0 call sites.
- `Injury.py:190-195,285-357` — `_postCalculatePointWise` (body is `pass`), `calculatePointWiseFractionInjured` — 0 call sites.
- `ProtectionPolicy.py:205-213` — `hdfkey` (at the policy level, not the action level) — 0 call sites.
- `riskToolkit.py:127-152` — `listAgentsNames`, `loadAgent` — 0 call sites, no CLI wiring.
- `analysis/riskAreas.py` — the entire module (`getRiskAreaAlgorithm`, `riskAreaAlgorithm_Sweep`) — 0 call sites and broken (see Urgent #7).

### simulations (the largest list — see also "Urgent" for the broken cases)
- `analysis/errorCalculation.py` — entire class `errorCalculation` (7 methods) — 0 usage anywhere in the repo.
- `gaussian/DropletCloud.py`, `gaussian/MeshUtils.py` — entire files, 0 imports anywhere.
- `hydrodynamics/nearWallFlow.py` — entire file/package, not registered in `toolkitHome` at all.
- `gaussian/gasCloud.py:668-716` — class `Continuous` ("Yehuda's Code") — 0 usage.
- `gaussian/Meteorology.py:317-322` — `MeteorologyProfile` — empty stub (`pass`), not registered in the factory.
- `windProfile/toolkit.py:115` — `_getStationsInRegion` — 0 call sites, also has a relative-path bug.
- `CLI.py:480,525` — `workflowNodes_list`, `workflowNodes_listParameters` — not wired into the `hera-workflows` CLI.
- `WRF/wrfDatalayer.py:16` — `wrfDatalayer` — re-exported from `__init__.py` but not registered in `toolkit.py`, never constructed anywhere. **Needs manual verification.**
- `LSM/hermesWorkflowToolkit.py` — **entire file (783 lines)**: `workflowToolkit` is an old, stale fork of `hera/simulations/hermesWorkflowToolkit.py` that nothing inherits from (`LSMToolkit` inherits directly from `abstractToolkit`).
- `deposition/models.py`, `evaporation/models.py`, `evaporation/monaghan.py` — classes `depositionModels`, `evaporationModels`, `MonaghanConstantConditions` — 0 external usage; `monaghan.py` also imports `pyriskassessment`, which isn't installed at all.
- `machineLearningDeepLearning/dataanalisys/` — **entire subpackage** (`stat()`, class `ml`) — 0 imports anywhere.
- `machineLearningDeepLearning/toolkit.py:323-415` — class `analysis` (`sensitivityAnalysis_morris`) — never instantiated in `__init__`.
- `LSM/toolkit.py:33-34`, `LSM/template.py:26-28` — constants `TRUE`/`FALSE`, `STABILITY_*` — 0 call sites (live code uses string literals directly).
- `openFoam/toolkit.py` — a long list of methods with no callers: `template_add` (497-520, body is `pass`), `clearVTKPipelineCache`, `getVTKPipelineCacheTable`, `getHermesWorkflow_Flow`, `getMeshFromName`, `getMeshExtentFromName`, `xarrayToSetFieldsDictDomain` (docstring literally says "Not debugged").
- `openFoam/postProcess/VTKpipelineExecutionContext.py` — **entire class** unused (0 references besides its own definition); redundant duplicate of `VTKPipeline.py`.
- `openFoam/postProcess/VTKPipeline.py:201-619` — class `registeredVTKPipeLine` (~420 lines) and `registerPipeline` — 0 call sites; `toolkit.py:814` only constructs the base `VTKPipeLine` without ever calling `registerPipeline`.
- `openFoam/preprocessOFObjects/OFList.py` — **entire class**, imported in `__init__.py` but never instantiated; also internally broken (`self.columnNames`/`getHeader()` don't exist).
- `openFoam/preprocessOFObjects/OFObject.py:22-23,59-81` — `REGION_INTERNSALFIELD` (also a typo in the name), `internalField()`, `.processors`, `.processorItems`, `.dimensionsStr` — 0 usage.
- `openFoam/lagrangian/abstractLagrangianSolver.py` — `getOriginalFlowFieldMesh`, `getDispersionDocument`, `getDispersionFlowDocument` — 0 call sites.
- `openFoam/lagrangian/LSM/toolkit.py:699-730` — `createRootCaseMeshLink` — 0 call sites.
- `eulerian/simpleFoam.py:4` — `simpleFoam_toolkitExtension` — not defined/wired in `OFToolkit.__init__` (unlike `buoyantReactingFoam_toolkitExtension`, which is wired).
- `openFoam/toberewritten/netcdf2of.py`, `xarrayDataset2OF.py` — **actual syntax errors** (`SyntaxError`/`TabError`), cannot be imported at all. `wrf2of.py` — not syntactically broken but has hardcoded absolute paths (`/data5/NOBACKUP/...`) and nothing imports it. (`toberewritten/utils.py`, however, **is** used by `LSM/toolkit.py` — do not delete.)

---

## 2️⃣ Commented-out code — no longer needed

Significant blocks (more than a few lines) recommended for deletion:

- `hera/simulations/openFoam/postProcess/pvOpenFOAMBase.py:643-718` — an entire `writeRegularCase` method commented out (~76 lines).
- `hera/simulations/openFoam/preprocessOFObjects/utils.py:52-270` — a huge block (~219 lines): `emptyParallelField`, `extractFile`, `extractBoundaryField` — also syntactically incomplete if uncommented.
- `hera/simulations/openFoam/OFWorkflow.py:445-935` — a huge block (~490 lines) of old geometry/blockMesh handler functions.
- `hera/simulations/openFoam/CLI.py:772-981` — an entire old copy (~210 lines) of `stochasticLagrangian_dispersion_create`, duplicating the live version above it.
- `hera/simulations/openFoam/lagrangian/abstractLagrangianSolver.py:1879-2056` — an old implementation (~178 lines) of `_extractFile`/`_readRecord`, superseded by newer functions.
- `hera/simulations/openFoam/lagrangian/LSM/toolkit.py:1190-1324` — includes **leftover unresolved git merge-conflict markers** (`# >> >> >> > b93a26...`) — evidence of a botched merge resolution that was never cleaned up.
- `hera/simulations/CLI.py:434-446` — a whole commented-out branch that contradicts the live logic below it.
- `hera/simulations/machineLearningDeepLearning/dataanalisys/ml.py` — several large blocks (lines ~38-51, 161-171, 454-479, 487-521) of disabled experimental code.
- `hera/measurements/GIS/utils.py:456-587` — **five** entire functions commented out: `ITMtolatlong`, `PolygonDataFrameIntersection`, `getBoundaries`, `makePolygonFromEndPoints`, `ConvexPolygons` (the last one an old duplicate of the live version in `buildings/analysis.py:33`).
- `hera/measurements/experiment/presentation.py:703-784` — two entire methods commented out: `plotNDIRFrequencyDistribution`, `plotMessageFrequencyDistribution`.
- `hera/measurements/meteorology/highfreqdata/analysis/meandatacalculator.py:452-472` — an entire `zoL` method commented out, replaced by `zOverL`.
- `hera/riskassessment/agents/effects/InjuryLevel.py:9-12,276-279` — old `unum` imports (the project migrated to `pint`) + an alternate implementation with a documented "bug".
- `hera/datalayer/collection.py:247-249`, `hera/datalayer/document/__init__.py:22` — small leftover scraps of old code.

---

## 3️⃣ Conditional logic that can never occur (Unreachable Code)

- `hera/simulations/gaussian/gasCloud.py:325-390` — duplicate class definition, immediately shadowed (see Urgent #1).
- `hera/simulations/openFoam/CLI.py` — three functions defined **twice**; Python only keeps the last definition: `foam_mesh_blockMesh` (259 vs. 640), `foam_mesh_setDomainHeight` (263 vs. 654), `IC_hydrostaticPressure` (267 vs. 708) — the first definitions are completely dead.
- `hera/measurements/GIS/vector/buildings/analysis.py:201-202` — `_LambdaF = 0`/`_LambdaP = 0` as class attributes, shadowed by later `def` definitions of the same names (385, 494) — the initial values can never be observed.
- `hera/measurements/GIS/vector/topography.py:50-95` — `cutRegionFromSource` calls `super()...` with keyword arguments that don't match the parent's signature (will raise `TypeError` on every call); the other branch calls `getRegionData`, which doesn't exist anywhere in the repo.
- `hera/simulations/hydrodynamics/nearWallFlow.py:174,298` — `Lambda = 2*numpy.log()` — called with no arguments, always `TypeError`.
- `hera/simulations/machineLearningDeepLearning/dataanalisys/ml.py:117,154,252` — three large blocks behind constant conditions like `if 4==5:`/`if 5==6:` — will never run.
- `hera/simulations/LSM/template.py:449,458` — `getSimulationByID`/`getSimulationByName` call methods that don't exist on `LSMTemplate` (`getDocumentByID`/`getSimulationsDocuments` — should be `self.Toolkit....`) — guaranteed `AttributeError`.
- `hera/simulations/openFoam/eulerian/abstractEulerianSolver.py:30` — the `flowType` property references a nonexistent constant.
- See also the "Urgent" section above for the most severe cases in this category.

---

## 4️⃣ Test files testing components that were deleted/are unused

**Good news:** no genuinely broken imports were found in the test files (`hera/tests/`) — every imported symbol checked has a matching class/function in the source.

One note: `hera/tests/test_datalayer.py:959` — `@pytest.mark.xfail(reason="Project.load references _iter_pickled_docs which is not implemented")` — the reason is **no longer accurate**: `_iter_pickled_docs` IS implemented, at `project.py:342`. The xfail itself isn't "broken," but its documented reason is stale — worth updating/removing the xfail if the test now passes.

**Significant coverage gap (not "broken," but "never tested at all"):** the entire `simulations/openFoam/{postProcess,preprocessOFObjects,toberewritten}` directories — **have no test file at all** in the repo. This explains why so many of the bugs found there (methods calling nonexistent functions) were never caught.

**Additional finding (not a "broken test" but a latent source bug):** in `hera/datalayer/autocache.py` — if the module is imported directly before the lazy loading in `hera/__init__.py` (`_load_deferred()`) has run, the line `from hera import Project` inside `autocache.py` triggers a re-entrant call to `_load_deferred()`, which tries to re-import `autocache` while it's still `partially initialized` → `ImportError: cannot import name 'cacheFunction'...`. In practice this doesn't currently blow up because the three relevant tests (`test_datalayer.py:898,917,931`) are marked `@requires_mongo` and are skipped in an environment without live MongoDB — but it would crash with a real MongoDB. Needs a fix in `hera/datalayer/autocache.py` (e.g., a local/lazy import instead of a module-level import).

---

## 5️⃣ Environment variables / config settings never read anywhere

| Variable | Read at | Set anywhere? |
|---|---|---|
| `HERA_REPO_ROOT` | `toolkit.py:771`, `experiment.py:141` | ✅ Yes (`set_hera_environment.sh`, `activate_hera.sh`) |
| `TEST_HERA` | `conftest.py:223` | ✅ Yes (`.github/workflows/ci.yml:22`) |
| `HERA_FULL_LOGGING_TESTS` | `logging/test_helpers.py:16` | ⚠️ Documented only, not set in CI — tests gated on it never run in CI |
| `RESULT_SET` | `conftest.py:245` | ⚠️ Documented only, not set in CI |
| `GDF_TOL_AREA` | `conftest.py:417` | ⚠️ Documented only, not set in CI |
| **`PYARGOS_PATH`** | `tests/dynamic_loading_tests_pack/test_experiment_cli_shortcuts.py:58` | ❌ **Never set anywhere** in the repo (not in CI, not in a Dockerfile, not in docs, not in .env) — the "deadest" of them all |

Unused constants:
- `hera/measurements/meteorology/highfreqdata/toolkit.py:26-27` — `DOCTYPE_STATIONS`, `DOCTYPE_MEASUREMENTS`.
- `hera/measurements/GIS/vector/toolkit.py:9` — `TOOLKIT_VECTOR_REGIONNAME`.
- `hera/simulations/LSM/toolkit.py:33-34`, `hera/simulations/LSM/template.py:26-28` — `TRUE`/`FALSE`/`STABILITY_*` constants.
- `hera/simulations/openFoam/preprocessOFObjects/OFObject.py:22-23` — `REGION_INTERNSALFIELD`/`REGION_BOUNDARYFIELD` (also a typo in the first name).
- `hera/simulations/openFoam/postProcess/pvOpenFOAMBase.py:12,31` — unused imports: `CASETYPE_RECONSTRUCTED`, `hera_logging`.

---

## Recommendations by priority

1. **Fix immediately** (live code that's broken): `gasCloud.py` (m/min), `thresholdGeoDataFrame.py` (collections.Iterable), `ProtectionPolicy.py` (os.JSONpath), `modelContainer.py` (ckpt_path), `abstractEulerianSolver.py`.
2. **Delete whole duplicate/dead files or folders**: `simulations/utils/interpolation/`, `simulations/LSM/hermesWorkflowToolkit.py`, `openFoam/postProcess/VTKpipelineExecutionContext.py`, `openFoam/preprocessOFObjects/OFList.py`, `openFoam/toberewritten/{netcdf2of,xarrayDataset2OF}.py`, `GIS/raster/hill2stl.py`.
3. **Delete stale commented-out blocks** (section 2 above) — low-risk, easy cleanup.
4. **Verify manually before deleting** (public API that may be used by external code): `wrfDatalayer`, the untested presentation/analysis methods, `PYARGOS_PATH`.
