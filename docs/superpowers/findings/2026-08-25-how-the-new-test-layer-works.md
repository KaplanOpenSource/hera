# How the new test layer works, and why we built it

A plain-language explainer for the new hermetic unit-test layer added under `hera/tests/unit/`. Not formal docs, not part of the mkdocs site — just "what is this and how do I use it."

---

## What existed before

Before this work, hera had **713 integration tests** (under `hera/tests/`, not `hera/tests/unit/`). They worked well, but:

- All of them require a **real MongoDB** running (`make mongo-up`).
- Some require **real S3 test data** (a few hundred MB, fetched via `make test` / `bootstrap_unittest_data.sh`).
- They test toolkits "from the top" — call `Project`/`Toolkit` and compare the final result to a baseline. That's great for checking "is the overall result correct," but it **doesn't exercise every function individually**.
- Result: **only 27% of the codebase's public functions (366 of 1,334) had ever been called by any test.** 968 functions never ran in any test — including functions with serious bugs (like referencing an undefined variable, or a class that breaks inheritance) that nobody could have caught without actually calling them.

## What we built

A **second, parallel layer — not a replacement** — at `hera/tests/unit/`, built on two ideas:

### 1. mongomock instead of a real MongoDB
`hera/tests/unit/_seam.py` replaces the function that connects to Mongo (`connectToDatabase`) with one that points at `mongomock` — an in-memory database that speaks the same protocol. **This is not a mock of hera's own code** — datalayer, Project, Toolkit — all of hera's own layers run in full, just against in-memory storage instead of a real server. That's what makes it possible to test things like "was the document saved under the right name" without Mongo at all.

### 2. Stubs for external packages that aren't installed
`hera/tests/unit/_stubs.py` registers placeholders for 7 packages (PyFoam, paraview, FreeCAD, hermes, argos, evtk) so that code importing them at module level (`import PyFoam...`) doesn't crash with `ModuleNotFoundError`. This makes it possible to test the **pure logic** inside openFoam/LSM/etc. even without the real software installed — at the cost that things which genuinely need PyFoam to do work (e.g. actually writing a real OpenFOAM file) stay only smoke-tested, not deeply verified.

**Result:** 1,370 new tests that run in **about 6-10 seconds**, with no Mongo, no S3, no network. Function coverage rose from 27% to **46.5%** (665 of 1,429 — the total count also grew slightly as small blocker fixes were made along the way). Along the way, **78 real bugs** were found — not theoretical issues, but code that crashes when actually run. Every one is pinned with a `@pytest.mark.xfail(strict=True)` test, so if someone fixes the bug later, the test will "fail" (because it expected a failure and got a success) and force the marker to be removed. That's how "this is a known bug" gets documented without fixing it right now, and without it being forgotten.

---

## How to run it

| Goal | Command | Needs Mongo/S3? | Time |
|---|---|---|---|
| Just the new layer (fast) | `make test-unit` | No | ~6-10s |
| New layer + its own coverage report | `make coverage-unit` | No | ~10s |
| **Everything** — both layers + combined coverage report + HTML | `make coverage` | **Yes** (fetches Mongo and S3 data if missing) | a few minutes |
| Just the old integration layer (as before) | `make test` | Yes | a few minutes |

For day-to-day work, `make test-unit` is enough for fast feedback. Before a push/PR, `make coverage` is worth running for the full picture.

**Important:** the two layers **cannot run in the same pytest process** — `hera/tests/unit/conftest.py` moves `HOME` and installs the seam at collection time, which affects the whole process. There's a guard that throws an explicit error if both are run together by mistake (`--ignore=hera/tests/unit` is what separates them above — not `-m "not unit"`, because 147 pre-existing tests are already marked `@pytest.mark.unit` and need to run in the integration layer).

---

## What runs automatically, and when

**The trigger is `pull_request` targeting `master` — and only that.** A plain `git push` to a branch triggers nothing. CI (`.github/workflows/ci.yml`) only starts when a PR targeting `master` is opened or updated.

When that happens, two jobs run **in sequence**:

### Job 1 — `unit` (no services)
```
pytest hera/tests/unit -m unit -q --cov=hera
```
No MongoDB, no S3. If this layer accidentally needs something external, it fails — **on purpose**: that's the proof it's actually hermetic.

### Job 2 — `test` (with a real MongoDB, needs: unit)
Only runs if Job 1 passes. Spins up a MongoDB container, fetches S3 test data (cached by version), then:
```
pytest hera/tests/ -m "not notebook" --ignore=hera/tests/unit --cov=hera
```
This is the entire old integration layer (713+ tests) running as before.

### At the end of Job 2 — **the coverage gate**
```
coverage combine   # merges the coverage from both layers (Job 1 + Job 2) into one report
coverage report
coverage report --fail-under=<the first number in coverage_floor.txt>
```
This is the step that can **fail the PR**. It checks neither the new layer alone nor the old layer alone — it checks the **combined** coverage, on purpose: gating on the new layer alone would reward mocking everything instead of testing real behavior.

**The current floor: 37%.** It's always set a few points below what's actually measured (38%) — a deliberate margin, since CI runs Python 3.11 with the versions pinned in `requirements.txt`, while this was measured on 3.12, so there's room for drift between environments. **The floor only ever goes up, never down** — that's the intent: it's meant to prevent regressions (a PR that lowers coverage fails), not to require every PR to raise it.

---

## Why this was necessary at all

Three reasons, briefly:

1. **699 functions ran for the first time in this project's history, and 78 of them crashed.** Not "returned the wrong answer" — crashed with `TypeError`/`AttributeError`/`NameError`. No line-coverage metric would have found these bugs, because they require someone to **actually call** the function. The most extreme example: `continuousReleaseGasCloud` — a whole class in the Gaussian dispersion code where all four of its public methods crash, because it inherits from the wrong parent class. Nobody knew.
2. **The coverage gate prevents silent regression.** Without a gate, new code can merge with no test at all, and coverage can quietly drift backward without anyone noticing. Now: if a PR drops the combined coverage below the floor, it fails automatically, before anyone has to remember to check.
3. **The unit layer is fast enough to run on every commit, not just before a release.** 6-10 seconds with no external infrastructure means a feedback loop you can run locally before every push, instead of relying only on CI that runs once per PR and needs MongoDB+S3.

---

## What's still not covered, and why that's okay (for now)

Of the 1,429 public functions, 764 are still untested. Most of it concentrates in areas that need infrastructure that wasn't available in the hermetic environment: a real OpenFOAM solver, VTK/paraview for rendering, a fully namespaced torch stub, real hermes/luigi for workflow orchestration, and real shapefiles/DEM data for GIS. Each one is documented in the findings file (`docs/superpowers/findings/2026-08-24-test-expansion-findings.md`) with a concrete reason it was deferred — not forgotten, just out of scope for this round.
