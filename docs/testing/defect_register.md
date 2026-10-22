# Defect Register

!!! warning "Not yet in master"
    The register and the tests that pin each defect live on `fix/groupAD` and the
    `tests/batchNN-*` branches. The fixes listed at the end live on separate `fix/*` branches.
    **None of it is merged into `master`.** Everything described here is still present in `master`.
    Pending merge through open PRs [#1067](https://github.com/KaplanOpenSource/hera/pull/1067) (batches 1–10), [#1083](https://github.com/KaplanOpenSource/hera/pull/1083) (11–20), [#1084](https://github.com/KaplanOpenSource/hera/pull/1084) (21–30) and [#1085](https://github.com/KaplanOpenSource/hera/pull/1085) (31–33).

While building the [hermetic unit layer](unit_tests.md), every public function was called directly,
many of them for the first time. Defects found this way were **recorded, not fixed**. This page
summarises that record.

## How defects are recorded

| Item | Convention |
|---|---|
| ID | `B<n>`, numbered in order of discovery across all 33 batches (B1–B332). Phase 0 infrastructure and process findings use `P1`–`P8` |
| Register | `docs/superpowers/findings/2026-08-24-test-expansion-findings.md` on `fix/groupAD`, a running log with one section per batch. **Written in Hebrew** |
| Entry | Heading `### B<n>. <summary>`, the file and line (`**קובץ:**`), a code excerpt, and the observed effect |
| Pin | Every defect that can be imported is pinned by a `@pytest.mark.xfail(strict=True, reason="B<n>: ...")` test that asserts the correct behaviour, usually paired with a "Characterisation of B<n>" test that asserts today's behaviour. See [Pinning known defects](unit_tests.md#pinning-known-defects-xfailstricttrue) |
| Verification | The register header states that each finding was confirmed by running code. Later batches (32–33) add "verified by reading source" where a stub prevents running it. A few, such as B286 in `hera/bin/hera-LSM`, cannot be pinned because the script has no `.py` extension and cannot be imported |
| Retraction | A finding that does not reproduce under the pinned `requirements.txt`, or that was fixed upstream in parallel, is marked in the register and its xfail turned into a passing assertion |

Related files on `fix/groupAD`:

- `docs/superpowers/findings/CONSOLIDATED-ISSUE-DRAFT.md`: a draft single GitHub issue covering B1–B78. Marked "draft, not opened".
- `docs/superpowers/findings/2026-08-25-how-the-new-test-layer-works.md`: a plain-language explainer.
- `docs/superpowers/plans/2026-08-24-batch1-utils.md`: the only source for B1–B4.

### Gaps in the numbering

| IDs | Where they are recorded |
|---|---|
| B1–B4 | Only in the batch 1 plan and in test `reason=` strings, not as register headings |
| B17 | Not found anywhere (looks like an unused number) |
| B82 | Commit `9be926ea` and `coverage_floor.txt` only. The register mentions it inside B89 (its twin) |
| B287–B295 | Only in test-file docstrings (`test_openfoam_eulerian_solver_more.py`, `test_openfoam_objecthome_more.py`, `test_openfoam_stochasticlagrangiansolver.py`, and others) |

The register itself has **317 numbered entries** (B5–B332, minus B17, B82 and B287–B295), plus P1–P8.

## Counts

### By batch

| Batches | IDs | Count |
|---|---|---|
| 1–9 | B1–B45 (no B17) | 44 IDs (`coverage_floor.txt` says "48 defects"; the difference is not explained) |
| 10–19 | B46–B78 | 33 |
| 20–29 | B79–B109 | 31 |
| 30 | B110–B141 | 32 |
| 31 | B142–B211 | 70 |
| 32 | B212–B309 | 98 |
| 33 | B310–B332 | 23 |

### By area

Classified by the source file named in each of the 317 register entries. An entry that names no file is
assigned to the file of the entry before it in the same section. Treat the counts as approximate.

| Area | Entries |
|---|---|
| OpenFOAM (`simulations/openFoam/**`, `hermesWorkflowToolkit.py`) | 75 |
| GIS (`measurements/GIS/**`) | 38 |
| Other simulations (deposition, evaporation, hydrodynamics, windProfile, WRF, `simulations/utils`, `simulations/CLI.py`) | 38 |
| Experiment (`measurements/experiment/**`) | 33 |
| Risk assessment (`riskassessment/**`) | 28 |
| Meteorology (`measurements/meteorology/**`) | 26 |
| `hera/utils` (incl. `utils/data`) | 26 |
| LSM (`simulations/LSM/**`, `hera/bin/hera-LSM`) | 20 |
| Gaussian (`simulations/gaussian/**`) | 16 |
| Core (`toolkit.py`, `datalayer/**`, `hera/bin/hera-project`) | 13 |
| Presentation (`presentation/basicplots.py`) | 3 |
| Test infrastructure (`test_no_invalid_escapes.py`) | 1 |

### Retracted or resolved elsewhere

| IDs | Status |
|---|---|
| B53, B60, B65, B84, B85, B86 | Real when recorded, then fixed independently upstream (commits `14508043`, `193249e6`, `db405cec`, `276d9a93`). Tests now assert the working behaviour |
| B58, B94 | **Retracted.** They reproduced only with a locally drifted pandas 3.0.2, not the pinned 2.2.3. B94 would come back with a pandas upgrade that makes copy-on-write the default |
| B42, B61 | Resolved by deletion: `hill2stl.py` and `OFList.py` were removed upstream as dead code |

## Most consequential defects

These were chosen for their effect on results or because they disable a whole feature. All of them
are in `master` today. "Fixed on" names a branch that is **not merged**.

| ID | Module | Effect | Fixed on |
|---|---|---|---|
| B12 | `utils/unitHandler.py` | `tounit` builds the Quantity in pint's default registry, not Hera's `ureg`, so Hera-only units such as `dunam` raise. This also breaks `InjuryLevelThreshold.getPercent` for numeric loads (B39) | — |
| B19 | `toolkit.py` | Add-document wrappers `setdefault` into the caller's dict. With a shared metadata dict, a second toolkit's document is saved under the first toolkit's name | — |
| B22 | `toolkit.py` | `setDefaultRepository` always raises `RuntimeError`, so the feature never worked | — |
| B27 | `gaussian/MeshUtils.py` | `sigmaYName` points at the X column, so the Gaussian cloud is forced isotropic | `fix/groupM-gaussian` |
| B31 | `evaporation`, `deposition` | `from ..utils import` resolves to an empty package, so both physics packages cannot be imported | `fix/groupAD` (commit `8eaa510a`) |
| B32 | `analysis/errorCalculation.py` | `calculateR` mixes N and N−1 normalisations. Every correlation is scaled by (N−1)/N, so a perfect model scores 0.667 at N=3 | — |
| B43 | `GIS/raster/landcover.py` | Two contradictory IGBP roughness tables in one file, differing ×100 for water and ×1300 for snow/ice. z0 feeds the wind profile and from there all dispersion | — |
| B59 | `protectionpolicy/ProtectionPolicy.py` | The indoor-dilution model writes into a copy returned by xarray indexing, so it always computes zero | — |
| B63 | `agents/effects/Injury.py` | No way to pass `inUnits`, so any pandas input to `calculateToxicLoads` raises `KeyError` | — |
| B66 | `riskassessment/Agents.py` | `getVolatility` returns g²/(mol·cm³) instead of g/cm³ | — |
| B70 | `highfreqdata/.../meandatacalculator.py` | The constructor compares `type(x)` to the *modules* `pandas`/`dask`, so building it from a DataFrame, as documented, always raises | `fix/groupZ` |
| B78 | `gaussian/gasCloud.py` | `continuousReleaseGasCloud` has the wrong base class. All four public methods raise, so continuous-release Gaussian dispersion is unusable | `fix/groupM-gaussian` |
| B82 | `highfreqdata/parsers/CampbellBinary.py` | TOB1 parser's column indices are off by two, dropping the first data column and crashing `getPandasFromFile`/`parse` (twin: B89) | — |
| B96 | `gaussian/FallingNonEvaporatingDroplets.py` | The constructor calls a `MeteorologyFactory` method that does not exist, so the class cannot be built | — |
| B110 | `windProfile/toolkit.py` | `getWindProfile` recomputes U* at the same z it evaluates, so the profile has no shear and z0 has no effect | — |
| B115–B117 | Campbell TOB1 reader (both copies) | Timestamp precision loss, an off-by-one in field conversion, and `key / 256` instead of `key // 256` | `fix/groupJ-parsers` |
| B126 | `GIS/CLI.py` | `hera-GIS` is entirely dead: `_setup` imports `WSG84`/`ITM` from the wrong module | `fix/groupE-data-handling` |
| B140 | `hermesWorkflowToolkit.py` | Every DB-backed workflow comparison passes `resource=` where hermes expects `Resource_path=` | `fix/groupO-mixed` |
| B153 | `utils/data/toolkit.py` | A `Cache` section in a repository JSON can never load: a parameter-name mismatch hidden by a bare `except` | — |
| B166 | `openFoam/CLI.py` | Handlers defined twice. The surviving bodies are `pass`, so mesh commands exit 0 and do nothing | `fix/groupF-openfoam-cli` |
| B183 | `GIS/raster/landcover.py` | Urban roughness writes an IGBP class code into z0: 12–14 m for rural cells instead of ~0.15 m | — |
| B187 | `GIS/raster/landcover.py` | `np.vectorize` without `otypes` infers int from cell (0,0), truncating every sub-metre z0 to 0 | `fix/groupAD` |
| B191 | `GIS/raster/topography.py` | `create_xarray`'s `minx`/`miny` mean (lat, lon) for WSG84 but (x, y) otherwise, silently returning a grid in the wrong place | — |
| B201 | `openFoam/lagrangian/LSM/toolkit.py` | `_readRecord` uses `self.logger`, which no base class defines, so the fallback OpenFOAM parser is unreachable | — |
| B232 | `GIS/vector/topography.py` | `addHeight` reads `TOOLKIT_SAVEMODE_*` from the wrong `toolkit` module, so every save mode raises and computed columns are discarded | `fix/groupAD` |
| B257 | `deposition/models.py` | Boltzmann constant written as 1.83e-23 instead of 1.380649e-23, making Brownian diffusivity 32.6% too high | `fix/groupE-data-handling` |
| B259, B260 | `toolkit.py` (ToolkitHome) | `auto_register_and_get` and `import_toolkits_from_json` call `registerToolkit` with the wrong signature | — |
| B261 | `datalayer/collection.py` | `addDocumentFromJSON` silently drops any document with `_id`. Since `exportProject` writes `_id`, re-importing an export imports nothing | `fix/groupP-more` |
| B263 | `analysis/riskAreas.py` | `riskAreaAlgorithm_Sweep` calls `.datalayer()` where the method is `project`, so the only risk-area algorithm never completes | — |
| B272 | `LSM/template.py` | `_getSimulationsList` cross-joins descriptions and parameters, attaching one run's id to another run's params | `fix/groupM-gaussian` |
| B280 | `LSM/toolkit.py` | `prepareSlurmLSMExecution` passes `desc=` instead of `template_desc=`, after already writing scripts | `fix/groupAA` |
| B286 | `hera/bin/hera-LSM` | `parser.parse_args()()` calls the Namespace, so the LSM CLI cannot dispatch any subcommand | — |
| B309 | `presentation/basicplots.py` | `seaborn.set()` at import rebinds matplotlib colour shorthands process-wide, depending on import order | — |
| B314 | `GIS/raster/tiles.py` | Mosaic sized `lr-ul` instead of `lr-ul+1`, so the last tile row/column is never fetched | `fix/groupN-gis2` |
| B320 | `openFoam/postProcess/VTKPipeline.py` | `clearCache` deletes from the Simulations collection, but the cache lives in Cache, so it clears nothing | — |
| B330 | `openFoam/lagrangian/abstractLagrangianSolver.py` | Four call sites use `getCaseListDocumentFromDB`, which is defined nowhere, so `createDispersionCaseDirectory` cannot run | `fix/groupX` |
| B332 | `openFoam/lagrangian/abstractLagrangianSolver.py` | `from hera.utils.unitHandler import *` replaces the builtin `min` with a unum unit, breaking time-step mapping | `fix/groupY` |

The recurring pattern, as the consolidated draft puts it: the most serious defects are not wrong
arithmetic but **code that never ran**: undefined names, wrong base classes, type checks that can
never be true. Line coverage alone would not find them.

## Fixes since the register was written

Fixes come as small `fix:` commits that cite B-numbers in the subject, one branch per group, each
built on top of the batch branches. Each fix removes the matching `xfail` marker and deletes the
characterisation test. **No fix branch is merged into `master`.**

Across all branches, **124 distinct B-numbers** are cited in `fix:` commit subjects, not counting retractions.

| Branch | B-numbers |
|---|---|
| `fix/groupAD` (this branch) | B187, B232. Also commit `8eaa510a`, which cites no IDs but fixes B31, B26, B41 and the import-time side effect B42 |
| `fix/groupA-unrunnable` | B69, B71, B72, B103, B127, B147 ¹, B248, B151, B152, B181 |
| `fix/groupB-error-paths` | B6, B133, B134, B327, B132 |
| `fix/groupC-plotting` | B143, B144, B224, B225, B322 |
| `fix/groupD-utils` | B5, B7, B10, B11, B15 |
| `fix/groupE-data-handling` | B126, B228, B254, B257, B323 |
| `fix/groupF-openfoam-cli` | B166, B167, B168, B171, B172, B173 |
| `fix/groupG-missing-names` | B18, B146, B158, B199, B233 |
| `fix/groupH-entangled` | B174, B176, B177, B295 (+ two unnumbered). Built on groupA |
| `fix/groupI-experiment` | B75, B76, B113, B162, B212, B213, B214 |
| `fix/groupJ-parsers` | B115, B116, B117, B118 |
| `fix/groupL-risk` | B24, B123, B124, B265 |
| `fix/groupM-gaussian` | B27, B28, B78, B272, B274 |
| `fix/groupN-gis2` | B190, B229, B314, B323 |
| `fix/groupO-mixed` | B88, B137, B140, B278 |
| `fix/groupP-more` | B87, B99, B261, B278, B294 |
| `fix/groupQ-workflow` | B102, B136, B138, B267 |
| `fix/groupR-more` | B268, B269, B270, B276, B277 |
| `fix/groupS-more` | B20, B139, B141, B262 |
| `fix/groupT-more` | B130, B219, B220, B273 |
| `fix/groupU-more` | B236, B238, B242 |
| `fix/groupV-more` | B217, B218, B223 |
| `fix/groupW-more` | B45, B74, B81 |
| `fix/groupX` | B296, B303, B304, B308, B325, B328, B330, B331 |
| `fix/groupY` | B302, B332 |
| `fix/groupZ` | B70, B80, B129 |
| `fix/groupAA` | B279, B280, B301, B306, B307 |
| `fix/groupAC` | B98 |
| `fix/bugs-structural` | B143 |
| `fix/structural-utils` | B5, B6, B7, B9, B10, B11 |
| `fix/structural-datalayer-cli` | B20, B87, B88, B262 |

¹ B147 is only partly fixed: the `AttributeError` is gone, but the test stays `xfail` with a new reason,
because the rest of the problem is inside `hermes.workflow`.

Some fix commits also record defects they attempted and then reverted: B9 in `fix/groupB` (later fixed on
`fix/structural-utils`), B33 deliberately not touched in `fix/groupA`, and B241 in `fix/groupAD`. The `fix/structural-*`
and `fix/bugs-structural` commits are dated 2026-10-20 and are still in progress.

Defects that are still open include several with the largest scientific effect: **B12, B19, B32, B43, B59, B63, B110, B183**.
