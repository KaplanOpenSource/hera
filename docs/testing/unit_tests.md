# Hermetic Unit Test Layer

!!! warning "Not yet in master"
    This page describes work that lives on the branch `fix/groupAD` and on the stacked
    `tests/batchNN-*` branches (batches 1–33). None of it is merged into `master` yet.
    File paths, commands and numbers below refer to `fix/groupAD` unless stated otherwise.
    Pending merge through open PRs [#1067](https://github.com/KaplanOpenSource/hera/pull/1067) (batches 1–10), [#1083](https://github.com/KaplanOpenSource/hera/pull/1083) (11–20), [#1084](https://github.com/KaplanOpenSource/hera/pull/1084) (21–30) and [#1085](https://github.com/KaplanOpenSource/hera/pull/1085) (31–33).

## Overview

The original Hera suite under `hera/tests/` is an **integration** suite: every test needs a live
MongoDB, most need the S3 test data set, and toolkits are exercised top-down against baselines
(see [Getting Started](getting_started.md) and [Test Flow](flow.md)).

The unit layer is a **second, parallel layer** under `hera/tests/unit/`. It does not replace the
integration suite. Its goal is to call individual functions directly, with:

- no MongoDB server,
- no S3 test data,
- no network.

When the work started, only 27% of Hera's public functions (366 of 1,334) had ever been called by
any test, according to the design notes in `docs/superpowers/findings/2026-08-25-how-the-new-test-layer-works.md`.
Calling them one by one is what found the defects listed in the [Defect Register](defect_register.md).

!!! tip "Key principle"
    Only the **database connection** and the **heavy external packages** are replaced. Hera's own
    code (`Project`, the datalayer, `abstractToolkit`, the toolkits themselves) runs for real.

---

## How it stays hermetic

### 1. mongomock instead of MongoDB (`_seam.py`)

`abstractToolkit` inherits from `Project`, and `Project.__init__` reads from and writes to the
database. So no toolkit can be built without one. Mocking every call site would mean the datalayer
was never really tested. The seam replaces **one function** instead:
`hera.datalayer.document.connectToDatabase`.

```python
# hera/tests/unit/_seam.py (abridged)
def install():
    import mongomock
    from mongoengine import connect, disconnect
    import hera.datalayer.document as hdoc

    def _mongomock_connect(mongoConfig, alias=None):
        ...
        return connect(alias=alias, db=mongoConfig["dbName"],
                       mongo_client_class=mongomock.MongoClient,
                       uuidRepresentation="standard")

    hdoc.connectToDatabase = _mongomock_connect
    hdoc.createDBConnection(connectionName=getpass.getuser(), mongoConfig=...)
    _rebind_module_level_collections()   # Measurements / Simulations / Cache / All
```

`createDBConnection` looks up `connectToDatabase` by its module-level name, so replacing that
attribute is enough. `hera/datalayer/__init__.py` creates the `Measurements`, `Simulations`, `Cache`
and `All` collection objects when it is imported, before the seam exists, so the seam creates them again.
Query building, desc filtering, version resolution and document tagging then run unchanged against an
in-memory store.

`requirements-test.txt` pins `mongomock==4.3.0` and `pytest-cov==5.0.0`. The seam uses
`mongo_client_class=`, not the `mongomock://` URI, because mongoengine removed that scheme.

### 2. `sys.modules` stubs for packages CI cannot install (`_stubs.py`)

The module docstring says nineteen Hera modules import PyFoam, paraview, FreeCAD, hermes, argos or evtk
at module level. `_stubs.install()` registers placeholder modules before those imports run, so
the pure-Python logic around them can be tested:

| Group | Modules | How it is stubbed |
|---|---|---|
| Namespace packages (always) | `PyFoam`, `PyFoam.RunDictionary`, `PyFoam.Basics`, `paraview`, `paraview.simple` | `types.ModuleType` with `__path__ = []`, so `import a.b` works |
| Leaf mocks (always) | `PyFoam.RunDictionary.ParsedParameterFile`, `PyFoam.RunDictionary.BoundaryDict`, `PyFoam.Basics.DataStructures`, `evtk`, `evtk.hl`, `dask.distributed` | `unittest.mock.MagicMock()` |
| Only if absent | `argos`, `argos.experimentSetup`, `argos.experimentSetup.dataObjects` | namespace packages, plus **real empty classes** `ExperimentZipFile`, `TrialSet`, `Trial`, `EntityType`, `Entity` |
| Only if absent | `hermes` | module whose `workflow` attribute is a real empty class |
| Only if absent | `FreeCAD` | `MagicMock()` |

Real classes are needed where Hera **subclasses** the stubbed name. `class Foo(mod.Bar)` needs
`Bar` to be a type. Subclassing a `MagicMock` attribute does not raise an error. It quietly produces
another MagicMock, and later `isinstance()` calls fail with a confusing `TypeError` (fixed in commit `41d5ae2b`).

!!! note "What is *not* stubbed"
    - **`torch`**: a single MagicMock is not enough, because `modelContainer.py` imports submodules
      (`'torch' is not a package`). It would need namespace-package treatment for `torch`,
      `torch.nn`, `torch.utils` and so on. That is deliberately out of scope, since `requirements.txt` declares the real package.
    - **GDAL / `osgeo`**: not stubbed. CI installs `libgdal-dev gdal-bin`. Tests that need `osgeo`
      or `sklearn` use `pytest.importorskip` (for example, the `GFS.py` import test).

    Because PyFoam/paraview are MagicMocks, code that really *writes* OpenFOAM files can only be
    smoke-tested. Batch 33's commit message says this directly: much of that batch checks argument
    forwarding and branch selection, not physics.

### 3. A fixed bootstrap order (`conftest.py`)

`hera/tests/unit/conftest.py` runs five statements at collection time. Their order matters:

1. Point `HOME` at a fresh `tempfile.mkdtemp(prefix="hera_unit_home_")`. Only the standard library may run before this.
2. Import `_seam` and `_stubs`. Both import only the standard library.
3. `_seam.write_pyhera_config()` writes a placeholder `~/.pyhera/config.json`. It has to exist before
   `hera.datalayer` is imported.
4. `_stubs.install()`.
5. `_seam.install()`. This is the first statement that imports `hera.datalayer`.

### 4. Guards that enforce isolation

All of these are autouse fixtures or hooks in `hera/tests/unit/conftest.py`:

| Guard | What it does |
|---|---|
| `_block_network` | Replaces `socket.socket` with a function that raises immediately. A test that tries to reach MongoDB fails at once instead of waiting out pymongo's 30 s timeout |
| `_no_trace_guard` | Overrides the parent conftest's session fixture of the same name, which purges projects in the real MongoDB. Before this override, every run took 62.95 s |
| `_no_home_writes` | Fails a test that creates `.hera` in the isolated home. Pass `filesDirectory` instead |
| `_reset_unit_database` | Drops the in-memory database after **every** test. Test order must never matter |
| `_matplotlib_agg` | Forces the `Agg` backend and calls `plt.close("all")` after each test |
| `pytest_collection_modifyitems` | Raises `pytest.UsageError` if anything outside `hera/tests/unit` is collected in the same process |

!!! danger "Two layers, two processes"
    The bootstrap affects the whole process. Running the unit and integration layers in a single
    pytest invocation made integration subprocesses authenticate with the placeholder credentials.
    `-m "not unit"` does **not** prevent this, because marker filtering happens after collection and
    the damage is done at import time. Always run the layers separately:

    ```bash
    pytest hera/tests/unit -m unit
    pytest hera/tests --ignore=hera/tests/unit -m "not notebook"
    ```

    The integration run must not also add `-m "not unit"`. The CI comments note that 147
    pre-existing integration tests carry `@pytest.mark.unit` and would then run in neither job.

---

## Directory layout

```
hera/tests/unit/
├── __init__.py
├── conftest.py          # bootstrap (HOME, config, stubs, seam), fixtures, guards
├── _seam.py             # mongomock seam: write_pyhera_config(), install(), reset()
├── _stubs.py            # sys.modules stubs: install()
├── _factories.py        # timeseries(), elevation_grid(), points_geodataframe()
├── test_infra_*.py      # 5 files that test the layer itself (seam, stubs, guards, isolation, require-data)
└── test_<area>_<module>[_more|_extra|...].py   # 185 behavioural test files
```

There are 190 `test_*.py` files on `fix/groupAD`. By filename prefix: `openfoam` 30, `utils` 23, `gis` 22,
`risk` 18, `meteorology` 18, `simulations` 13, `gaussian` 13, `experiment` 11, `datalayer` 11,
`toolkit` 5, `lsm` 5, `infra` 5, with smaller groups (`hermes`, `deposition`, `evaporation`, `windprofile`,
`riskassessment`, `wrf`, `ml`, `presentation`, `risktoolkit`) making up the rest.

---

## Running it

### Make targets

| Goal | Command | Needs MongoDB / S3? |
|---|---|---|
| Unit layer only | `make test-unit` | No |
| Unit layer + its own coverage report | `make coverage-unit` | No |
| Both layers + combined report + HTML (`cache/coverage_html`) | `make coverage` | Yes |
| Integration layer only (unchanged behaviour) | `make test` | Yes |

What the targets run (from the `Makefile` on `fix/groupAD`):

```bash
# make test-unit
PYTHONPATH=. pytest hera/tests/unit -m unit -q

# make coverage-unit
COVERAGE_FILE=.coverage.unit pytest hera/tests/unit -m unit -q --cov=hera --cov-report=
COVERAGE_FILE=.coverage.unit coverage report

# make coverage   (unit half, then integration half, then combine)
COVERAGE_FILE=.coverage.unit        pytest hera/tests/unit -m unit -q --cov=hera --cov-report=
COVERAGE_FILE=.coverage.integration pytest hera/tests -m "not notebook" --ignore=hera/tests/unit -q --cov=hera --cov-report=
coverage combine && coverage report && coverage html -d cache/coverage_html
```

Each layer writes its own `COVERAGE_FILE`. pytest-cov already merges its own parallel fragments, and
`coverage combine` exits 1 when there is nothing to merge, so only `make coverage` runs combine.

Install the test-only dependencies first. They are deliberately kept out of `requirements.txt`:

```bash
pip install -r requirements-test.txt
```

### Markers and options

`pytest.ini` declares four markers:

| Marker | Meaning |
|---|---|
| `unit` | No external dependencies. Every class or function in `hera/tests/unit` carries it |
| `integration` | Requires MongoDB |
| `slow` | Deselect with `-m "not slow"` |
| `notebook` | Notebook execution. Run with `make test-notebooks`, excluded from `make test` |

`addopts = -v --tb=short`. The integration layer's `hera/tests/conftest.py` also adds
`--require-data`, or the `HERA_REQUIRE_TEST_DATA=1` environment variable. Without it, missing
`TEST_HERA` data **skips** tests (the old local behaviour). With it, missing data **fails** them. CI arms it so a partial S3 download
cannot turn green.

### What CI runs

`.github/workflows/ci.yml` runs on `pull_request` targeting `master` only. Both jobs use
**Python 3.11** on `ubuntu-latest`:

1. **`unit`** (no services): `pytest hera/tests/unit -m unit -q --cov=hera --cov-report=` with
   `COVERAGE_FILE=.coverage.unit`. The data file is uploaded as an artifact.
2. **`test`** (`needs: unit`, MongoDB service container, S3 data cached by version, hermes and pyargos
   checked out onto `PYTHONPATH`):
   `pytest hera/tests/ -v -m "not notebook" --ignore=hera/tests/unit --cov=hera --cov-report=`.
   It then downloads the unit coverage, runs `coverage combine`, and enforces the floor (next section).

!!! note "Changes pending on other branches"
    - `tests/batch10-pure-physics` (PR #1067) has later review commits that `fix/groupAD` lacks.
      Commit `de2998de` deleted `coverage_floor.txt` and made **every** stub conditional, so a real
      install always wins. Commit `2554e5d6` moved the floor inline into CI (`coverage report --fail-under=30`).
      Commit `ab6264da` gives the `unit` job real hermes/pyargos checkouts. According to that commit, 106 of
      111 failures and 278 of 286 errors on `tests/batch33-coverage90` in CI came from the hermes/argos
      placeholder stubs.
    - `fix/bugs-structural` (issue #1125) adds a CI matrix of Python 3.11 and 3.14 against MongoDB 9.0.

---

## The coverage gate

- **Configuration:** `.coveragerc` sets `branch = True`, `source = hera`, `parallel = True`, and omits
  `hera/tests/*`, `hera/doc/*`, `hera/bin/*` and four paths that do not parse as Python.
  The gated figure is therefore **combined line + branch** coverage.
- **Floor:** the first line of `coverage_floor.txt` is the floor. CI reads it with
  `grep -oE '^[0-9]+' coverage_floor.txt | head -1` and runs `coverage report --fail-under="${FLOOR}"`.
- **Combined, not per layer:** the gate checks the merge of both layers. The Makefile and CI comments
  explain that gating on the unit layer alone would reward mocking every boundary.
- **Never lower it.** The header of `coverage_floor.txt` reads: *"Raise it in the same PR that raises
  coverage; never lower it."* The floor is deliberately set about two points below the measured value.
  Measurements were taken on Python 3.12, while CI runs 3.11 with the pinned `requirements.txt`, so the
  gate should fail on regressions and not on environment drift.
- **History:** every raise is recorded below the first line of `coverage_floor.txt`, so the file is also
  the measurement log. See [Batch History](batch_history.md).

---

## Pinning known defects: `xfail(strict=True)`

These batches changed **no production code**. A defect found by a test is not fixed in place. It is
recorded with a test that asserts the **correct** behaviour and is marked
`@pytest.mark.xfail(strict=True)`. Because the marker is strict, fixing the defect turns the test into an
XPASS, which fails the run. Whoever fixes the bug therefore has to remove the marker.

Most pins come with a **characterisation** test that asserts what the code does *today*. The fixer
deletes it together with the marker. A real example from
`hera/tests/unit/test_gaussian_gascloud.py` (defect B78):

```python
@pytest.mark.unit
class TestContinuousReleaseGasCloudIsUnusable:
    """B78: continuousReleaseGasCloud inherits directly from abstractGasCloud
    -- NOT from instantaneousReleaseGasCloud -- but all four of its own
    public methods call self.getDosage_inst_noQ / self.getDosage_inst_NoERF_noQ,
    which are only ever defined on instantaneousReleaseGasCloud. ..."""

    @pytest.mark.xfail(
        strict=True,
        reason="B78: getConcentration_cont calls self.getDosage_inst_noQ, "
               "defined only on instantaneousReleaseGasCloud. "
               "See the consolidated findings issue.",
    )
    def test_concentration_cont_should_scale_the_erf_dosage_by_the_release_rate(self, cloud):
        cloud.getConcentration_cont()

    def test_all_four_public_methods_currently_raise(self, cloud):
        """Characterisation of B78."""
        for method_name in ("getConcentration_cont", "getConcentration_cont_NoERF",
                            "getDosage_cont_NoERF", "getDosage_cont_doubleNoERF"):
            with pytest.raises(AttributeError, match="getDosage_inst"):
                getattr(cloud, method_name)()
```

Conventions:

- The `reason=` string always starts with the defect ID (`"B78: ..."`), so you can grep for it.
- The ID refers to an entry in the [Defect Register](defect_register.md).
- A defect that turns out not to exist in the pinned environment is **retracted**: the marker is
  removed and the test asserts the working behaviour (examples: B58, B94).

---

## Adding tests or a new batch

### Writing a test

1. **File:** `hera/tests/unit/test_<area>_<module>.py`, for example `test_gaussian_gascloud.py`
   or `test_meteorology_meandatacalculator.py`. Later passes over the same module add a suffix:
   `_more`, `_extra`, `_remaining`.
2. **Class:** `Test<Subject>`, decorated with `@pytest.mark.unit`. Without the marker the test is
   not selected by `-m unit`.
3. **Fixtures** from `hera/tests/unit/conftest.py`:

    | Fixture | Gives you |
    |---|---|
    | `unit_files_directory` | a per-test `tmp_path/heraFiles` directory |
    | `unit_project` | a real `Project("UNIT_TEST_PROJECT")` on mongomock |
    | `unit_toolkit_factory` | `factory(toolkitHome.<CONSTANT>)` builds any real toolkit through `toolkitHome.getToolkit` |

    ```python
    @pytest.mark.unit
    def test_real_toolkit_is_usable(unit_toolkit_factory):
        from hera import toolkitHome
        toolkit = unit_toolkit_factory(toolkitHome.GIS_RASTER_TOPOGRAPHY)
        toolkit.addDataSource("SRC", "/rel/path.tif", "GeoTIFF", version=[0, 0, 1])
        names = [d.desc.get("datasourceName") for d in toolkit.getDataSourceDocuments("SRC")]
        assert names == ["SRC"]
    ```

4. **Synthetic inputs:** use `hera/tests/unit/_factories.py` (`timeseries`, `elevation_grid`,
   `points_geodataframe`), or build small files in `tmp_path`. The commit messages for batches 30–31 note that most
   modules once listed as "needs real data" needed only small synthetic files: TOB1 binaries, GeoTIFFs, an argos v3.0.0 zip.
5. **Hygiene rules** from the design spec: `sorted()` around every glob, `pytest.approx` for floats,
   explicit edge cases (empty, NaN, single row, constant column), and `tmp_path` for every disk write.
6. **Monkeypatch classes, not the `toolkitHome` instance.** Batch 30 found that patching `getToolkit`
   on the instance left a permanent instance attribute after teardown. That silently broke 64 tests in another file,
   depending on collection order.
7. **Found a bug?** Do not fix production code in a test batch. Assign the next free B-number,
   add it to the findings log, and pin it with `xfail(strict=True)` plus a characterisation test.

### Closing a batch

The pattern used by batches 1–19:

1. Run `make coverage` and record the combined figure, and any per-area figures, as a new
   `# After batch NN (...)` block in `coverage_floor.txt`.
2. Raise the number on line 1 only if the gain leaves headroom (batch 2 kept the floor at 26 after a 0.4-point gain).
   Never lower it.
3. Record new defects in `docs/superpowers/findings/2026-08-24-test-expansion-findings.md`.
4. Branch the next batch from the previous one: `tests/batchNN-<scope>`.

---

## Current numbers

| Measure | Value | Source |
|---|---|---|
| Coverage floor (line 1) | **37** | `coverage_floor.txt` on `fix/groupAD` |
| Last recorded combined line+branch coverage | **38%** (18,592 stmts, 10,873 missed) | `coverage_floor.txt`, after batch 19 |
| Last recorded function coverage | **46.5%** (665 of 1,429 public functions) | `coverage_floor.txt`, after batch 19 |
| Unit coverage reported after batch 33 | **~90%** (from 83.4%) | batch 33 commit `d68d3ac1` |
| Passing unit tests after batch 33 | **5,600** | batch 33 commit `d68d3ac1` |
| Latest run on `fix/groupAD` | 5,579 passed, 376 xfailed, 39 failed | fix commit `9bb0170e` |
| Test functions (static count) | 5,502 `def test_` in 190 files, 349 `strict=True` markers | `git grep` on `fix/groupAD` (parametrisation expands the collected count) |
| Defects recorded | B1–B332 | [Defect Register](defect_register.md) |

!!! warning "Two different coverage measures"
    `coverage_floor.txt` was last updated after **batch 19**. Batches 20–29 were not re-measured, because the
    integration layer could not run locally. Batches 30–33 report a **"unit coverage"** percentage
    (55.0 → 62.4 → 75.4 → 83.4 → ~90) in their commit subjects, but the commit messages do not
    define it. It is not the combined line+branch figure the gate enforces, since the floor stayed at 37
    through batch 33. The figures continue the function-coverage series (46.5% after batch 19), but that
    is an inference. Re-measure with `make coverage` before raising the floor.

The 39 failures in the latest run come from the environment, not from the defects (the related
fix commit `72073d3f` calls them "the same 39 environmental failures"). Commit `21b766db` on
`fix/groupH-entangled` brings them down to 25.
