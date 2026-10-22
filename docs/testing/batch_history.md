# Unit-Test Batch History

!!! warning "Not yet in master"
    The batches below live on the stacked branches `tests/batch1-utils` … `tests/batch33-coverage90`
    and are collected on `fix/groupAD`. None of them is merged into `master` yet.
    Pending merge through open PRs [#1067](https://github.com/KaplanOpenSource/hera/pull/1067) (batches 1–10), [#1083](https://github.com/KaplanOpenSource/hera/pull/1083) (11–20), [#1084](https://github.com/KaplanOpenSource/hera/pull/1084) (21–30) and [#1085](https://github.com/KaplanOpenSource/hera/pull/1085) (31–33).

The [hermetic unit layer](unit_tests.md) was built in a Phase 0 (infrastructure) followed by 33
**batches**. Each batch targets one area of the codebase, adds tests, and pins every defect it finds
with `xfail(strict=True)` (see the [Defect Register](defect_register.md)). Each batch branch is
created from the previous one.

**Sources:** the measurement log below line 1 of `coverage_floor.txt` (batches 1–20), the commit messages on
`origin/master..fix/groupAD`, and the per-branch commit lists (`git log tests/batch(N-1)..tests/batchN`).

## Phases

| Phase | Batches | Dates (commits) | Outcome |
|---|---|---|---|
| Phase 0 – infrastructure | – | 2026-08-24 | mongomock seam, `sys.modules` stubs, isolation guards, coverage config, CI split. Baseline: 22% combined, floor 20 |
| Phase 1 | 1–9 | 2026-08-24 → 08-25 | 1,000 hermetic tests, floor 20 → 30, combined 22% → 32%, "48 defects" (per `coverage_floor.txt`; B1–B45 without B17 is 44 IDs) |
| Phase 2 ("phase 2 sweep") | 10–19 | 2026-08-25 | 1,370 tests, floor 30 → 37, combined 32% → 38%, function coverage 27% → 46.5%, defects through B78 |
| Function-by-function pass | 20–29 | 2026-08-30 → 09-02 | Local-only and not re-measured. Defects B79–B109. Several earlier findings retracted after checking against the pinned `requirements.txt` |
| Coverage push | 30–33 | 2026-09-03 → 09-07 | 5,600 passing tests, "unit coverage" 55.0% → ~90%, 223 new defects (B110–B332) |

The names "Phase 0", "phase 2 sweep" and the 1–9 / 10–19 boundaries come from `coverage_floor.txt`.
The labels for batches 20–29 and 30–33 are this page's own grouping.

## Batches

Column notes:

- **Coverage after** is the *combined* line+branch figure from `coverage_floor.txt` unless marked.
  For batches 30–33 it is the "unit coverage" from the commit subject. That metric is not defined in
  the commits and is **not** comparable with the combined figure (see the [note on the unit-tests page](unit_tests.md#current-numbers)).
- **Tests added** is given only where a source states it. "Σ" marks a total across several batches.
- "—" means the source does not give a value.

| # | Branch | Scope | Tests added | Coverage after | Floor after | Defects |
|---|---|---|---|---|---|---|
| 0 | (on `tests/batch1-utils`) | Infrastructure | — | 22% | 20 | P1–P8 |
| 1 | `batch1-utils` | `hera/utils`: angles, query, Filter, `_LazyModule`, slurm, latex/bibtex, zip, unitHandler, JSON config, logging | Σ 1–9: 1,000 | 28% | 26 | B1–B16 (B1–B4 in the batch plan only) |
| 2 | `batch2-datalayer` | `toolkit.py` + datalayer: abstractToolkit datasources/documents, datahandler registry, Project config, ToolkitHome | — | 28% | 26 | B18–B23 |
| 3 | `batch3-gaussian` | `simulations/gaussian`: Sigma (Briggs 1973), Meteorology, MeshUtils | — | 30% | 28 | B24–B30 |
| 4 | `batch4-physics-utils` | evaporation, deposition, hydrodynamics, windProfile, `simulations/utils` | — | 31% | 29 | B31–B35 |
| 5 | `batch5-riskassessment` | riskassessment + presentation: toxic-load calculators, injury levels | — | 31% | 29 | B36–B39 |
| 6 | `batch6-meteorology` | `measurements/meteorology`: seasonal classification, import surface | — | 31% | 29 | B40–B41 |
| 7 | `batch7-gis` | `measurements/GIS`: CRS conversion, land-cover roughness tables | — | 31% | 29 | B42–B43 |
| 8 | `batch8-experiment` | `measurements/experiment`: pydoc factory, orphaned TOA5 parser | — | — | 29 | B44 |
| 9 | `batch9-simulations` | openFoam, LSM, mlDL, WRF: LSM parameter preparation, stub layer end to end | (Σ 1,000) | 32% | 30 | B45 |
| 10 | `batch10-pure-physics` | evaporation/deposition physics, opened up by the import fix `8eaa510a` | Σ 10–12: 82 | — | — | B46–B51 |
| 11 | `batch11-openfoam` | OpenFOAM preprocessing objects | (Σ 82) | — | 30 | B52 |
| 12 | `batch12-risk-pure` | riskassessment: protectionpolicy, riskAreas, riskToolkit | (Σ 82) | 33% | 31 | B53–B60 |
| 13 | `batch13-openfoam-deep` | OFObject, OFList, `preprocessOFObjects/utils.py` | 28 | 34% | 32 | B61 |
| 14 | `batch14-risk-agents` | `riskassessment/agents`: Agents, Calculator, Injury, InjuryLevel, thresholdGeoDataFrame | Σ 14–19: 260 ¹ | 34% | 33 | B62–B66, B97 ² |
| 15 | `batch15-openfoam-field` | `OFField.py`, the remaining five methods | (Σ) | 35% (fn 40.1%) | 34 | B67–B68 |
| 16 | `batch16-meteorology` | highfreqdata: abstractcalculator, meandatacalculator, turbulencestatistics, parsers | (Σ) | 37% (fn 43.1%) | 36 | B69–B73 |
| 17 | `batch17-gis` | GIS, partial: stlFactory, TilesToolkit tile math, LandCoverToolkit (hill2stl tests later dropped) | (Σ) | 37% (fn 44.0%) | 36 | B74 |
| 18 | `batch18-experiment` | experimentHome, experimentAnalysis, dataEngine | (Σ) | measured with 19 | 36 | B75–B77 |
| 19 | `batch19-gaussian-gascloud` | `gasCloud.py` (deferred from batch 3) | (Σ) | 38% (fn 46.5%) | 37 | B78 |
| 20 | `batch20-gis-buildings` | GIS buildings/vector/demography morphology, Campbell TOB1 parser | 105 | not re-measured | 37 | B79–B82 |
| 21 | `batch21-openfoam-toolkit` | OFToolkit filesystem/parsing helpers | Σ 20–29: ~507 ³ | — | 37 | B83 |
| 22 | `batch22-openfoam-cli` | `Foam_parser_FieldDescription`, eulerian solver extension | (Σ) | — | 37 | B84–B86 (retracted) |
| 23 | `batch23-lsm-lagrangian` | Lagrangian particle-source generators, point-wise concentration | (Σ) | — | 37 | none |
| 24 | `batch24-datalayer` | remaining `DataHandler_*` format handlers | (Σ) | — | 37 | B87 |
| 25 | `batch25-datalayer-project` | remaining `datalayer/project.py` query/save helpers | (Σ) | — | 37 | B88 |
| 26 | `batch26-experiment-parsers` | `experiment/parsers.py`, `utils/data` repository export, lowfreq presentationLayer | (Σ) | — | 37 | B89–B91 |
| 27 | `batch27-lsm` | LSMToolkit construction, templates, `getTemplatesTable` | (Σ) | — | 37 | B92 |
| 28 | `batch28-simulations-utils` | `interpolations.py`, FallingNonEvaporatingDroplets, `utils/data/CLI.py` helpers | (Σ) | — | 37 | B93–B96 (B94 retracted) |
| 29 | `batch29-datahandler-more` | 30 commits: datahandler tif/npz/Class, LSM singleSimulation/template/sources, autocache, DropletCloud, OFToolkit VTK cache, hermesWorkflowToolkit names, OFLSMToolkit, highfreq parsers & HighFreqToolKit, risk CLI | (Σ) | "55.0%" ⁴ | 37 | B98–B109 |
| 30 | `batch30-coverage60` | windProfile toolkit, gaussian & experiment presentation, experiment parsers, Injury/InjuryFactory, droplet physics, hermesWorkflowToolkit DB lifecycle, GIS/experiment/simulations CLIs | 486 (1,877 → 2,363) | 62.4% ⁴ | 37 | B110–B141 (32) |
| 31 | `batch31-coverage72` | `utils/data` CLI + dataToolkit, openFoam CLI/toolkit/OFWorkflow, Lagrangian LSM toolkit, GIS raster landcover/topography, `experiment.py` with an argos v3.0.0 zip, buildings & basicplots presentation | 1,559 (1,863 → 3,422) ⁵ | 75.4% ⁴ | 37 | B142–B211 (70) |
| 32 | `batch32-coverage82` | 50 new files, 170 reachable functions: experiment analysis/presentation, GIS vector/raster, physics/utils tail, ToolkitHome, datahandler, risk tail, LSM, gaussian, meteorology | 1,605 (3,422 → 5,027) | 83.4% ⁴ | 37 | B212–B309 (98) |
| 33 | `batch33-coverage90` | `pvOpenFOAMBase.py`, TilesToolkit, VTKPipeline, `abstractLagrangianSolver.py` | 573 (5,027 → 5,600) | ~90% ⁴ | 37 | B310–B332 (23) |

**Footnotes**

1. Derived: 1,370 tests after batch 19 (`coverage_floor.txt`), minus 1,000 after batch 9, minus 82 and 28.
2. B97 was found on 2026-09-01 (commit `636f7d1a`) while retracting B65 in the batch 14 code.
3. Derived: 1,877 (start of batch 30, per its commit) minus 1,370. Approximate, because tests were also
   *removed* in this period (hill2stl and OFList coverage was dropped after those files were deleted upstream).
4. "Unit coverage" as stated in the batch 30–33 commit subjects (batch 30's subject gives 55.0% as the
   starting value, which is used here for batch 29). The metric is not defined in the commits.
5. **Sources disagree:** batch 30 ends at 2,363 tests, but batch 31 starts from 1,863. The 31 → 33 chain
   (3,422 → 5,027 → 5,600) is internally consistent and is the one used for current totals.

## Follow-up commits on `fix/groupAD`

| Commit | What |
|---|---|
| `8eaa510a` | First production change: fixes a circular import in `autocache.py`, wrong-level `from ..utils` imports (B31), import-time side effects (hill2stl demo, a hard-coded path in `highfreqdata/__main__.py`), and `test_no_invalid_escapes.py` never failing (B26). The commit reports "704 passed, 0 failed" for the integration suite |
| `ea25d331`, `636f7d1a`, `3672df49`, `bcd299e0` | Retractions after checking against the pinned `requirements.txt`: B53/B58/B60, B65, B84–B86, B94 |
| `5263ce30`, `41d5ae2b` | Real placeholder classes for `argos.experimentSetup.dataObjects` and `hermes.workflow` in `_stubs.py` |
| `4c267d67` | CI: install `requirements-test.txt` in the integration job too |
| `9bb0170e` | `fix:` B187 and B232. The commit reports 5,579 passed, 376 xfailed, 39 failed |

The other `fix/group*` branches, each fixing a few B-numbers on top of the batches, are listed in the
[Defect Register](defect_register.md#fixes-since-the-register-was-written).
