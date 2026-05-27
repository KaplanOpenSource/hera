# CLAUDE.md — Hera Development Guide

This file provides Claude Code with the context, conventions, and rules needed to work effectively on the Hera codebase.

---

## What is Hera?

Hera is a Python scientific data management platform (v2.16.1). It provides a unified MongoDB-backed data layer and a set of domain-specific **Toolkits** for GIS, meteorology, atmospheric dispersion (CFD + Lagrangian), and risk assessment.

**Repo:** `github.com/KaplanOpenSource/hera`  
**Language:** Python 3  
**Database:** MongoDB (3 collections: Measurements, Simulations, Cache)  
**Key dependencies:** pandas, dask, geopandas, xarray, pint, mongoengine, luigi, OpenFOAM (external)

---

## Architecture Overview

The system has three core abstractions — always keep these in mind:

```
ToolkitHome (singleton registry)
    └── abstractToolkit (base class, inherits from Project)
            ├── Data layer    — datasource management (inherited from Project)
            ├── Analysis      — domain-specific processing (toolkit.analysis)
            └── Presentation  — domain-specific plots (toolkit.presentation)

Project
    ├── Measurements_Collection  (observational data, toolkit datasources)
    ├── Simulations_Collection   (simulation outputs)
    └── Cache_Collection         (intermediate results, configs)
```

**Never** instantiate toolkits directly. Always use `toolkitHome.getToolkit(...)`.

---

## Package Layout

```
hera/
├── toolkit.py                         # ToolkitHome + abstractToolkit
├── datalayer/                         # Project, collections, datatypes
├── measurements/
│   ├── GIS/
│   │   ├── raster/                    # TopographyToolkit, LandCoverToolkit, TilesToolkit
│   │   └── vector/                    # VectorToolkit, BuildingsToolkit, DemographyToolkit
│   └── meteorology/
│       ├── lowfreqdata/               # lowFreqToolKit
│       └── highfreqdata/              # HighFreqToolKit
│   └── experiment/                    # experimentHome
├── simulations/
│   ├── openFoam/                      # OFToolkit (composition pattern)
│   ├── LSM/                           # LSMToolkit
│   ├── gaussian/                      # gaussianToolkit
│   ├── windProfile/                   # WindProfileToolkit
│   └── hermesWorkflowToolkit.py       # workflow base class
├── riskassessment/                    # RiskToolkit, agents, effects, policies
├── utils/                             # units (pint), angles, JSON, logging, Slurm
└── bin/                               # CLI entry points (hera-project, hera-toolkit, ...)
```

---

## Coding Rules

### General

- All new code must be Python 3. No Python 2 compatibility.
- Follow existing module structure — put new toolkits under the appropriate domain folder (`measurements/`, `simulations/`, etc.).
- Inherit from `abstractToolkit` for new toolkits. Use `hermesWorkflowToolkit` if the toolkit needs Luigi-based workflow support.
- Never modify `hera/toolkit.py`'s hardcoded `_toolkits` dict to add external toolkits — use dynamic registration instead.
- Use `pydoc.locate`-compatible class names when naming subclasses of `Calculator`, `InjuryLevel`, `Injury`, and `Action` in the risk assessment domain — naming convention is enforced by the factory.

### Toolkit Structure

Every new toolkit must follow this pattern:

```python
from hera.toolkit import abstractToolkit

class MyToolkit(abstractToolkit):
    toolkitName = "MyToolkit"

    def __init__(self, projectName, **kwargs):
        super().__init__(projectName=projectName, toolkitName=self.toolkitName, **kwargs)
        self._analysis = MyAnalysis(self)
        self._presentation = MyPresentation(self)

    @property
    def analysis(self):
        return self._analysis

    @property
    def presentation(self):
        return self._presentation
```

Analysis and presentation layers are separate classes that hold a back-reference to the parent toolkit. They must not be subclasses of the toolkit.

### Data Layer

- Use the correct collection for each data type:
  - `addMeasurementsDocument` → observational data, toolkit datasources
  - `addSimulationsDocument` → simulation outputs
  - `addCacheDocument` → intermediate results, configs
- Never mix collection types. Simulation results do not go in Measurements.
- Always use `isRelativePath: True` in repository JSON files. Never hardcode absolute paths.
- Always call `setDataSourceDefaultVersion(name, version)` after registering a datasource that may have multiple versions.

### Data Formats

| Use case | Format |
|----------|--------|
| Tabular data | `parquet` |
| Spatial vector data | `geopandas` |
| Multi-dimensional arrays / time series | `netcdf_xarray` |
| Small configs / metadata | `JSON_dict` |
| Directory paths | `string` |

Prefer `parquet` over CSV. Prefer lazy (Dask) evaluation for large datasets — never call `.compute()` before filtering.

### Naming Conventions

| Item | Convention | Example |
|------|-----------|---------|
| Project names | UPPERCASE with underscores | `WIND_ANALYSIS_2024` |
| Datasource names | UPPERCASE or lowercase_underscore | `YAVNEEL`, `tel_aviv_station` |
| Repository names | lowercase with underscores | `meteo_data_v1` |
| Toolkit class names | PascalCase + `Toolkit` suffix | `MyDomainToolkit` |
| Analysis class names | PascalCase + `Analysis` suffix | `MyDomainAnalysis` |
| Presentation class names | PascalCase + `Presentation` suffix | `MyDomainPresentation` |
| CLI scripts | `hera-<domain>` | `hera-myDomain` |
| Risk subclasses | `Calculator<Name>`, `InjuryLevel<Name>`, `Action<Name>` | `CalculatorThermal` |

### Units

Always use `pint` via Hera's `ureg` for physical quantities. Never use bare floats for dimensional values when units matter.

```python
from hera.utils import ureg
mass = 1 * ureg.kg
concentration = 5 * ureg.mg / ureg.m**3
```

Use `toMeteorologicalAngle` / `toMathematicalAngle` from `hera.utils` when converting wind directions. Never do angle arithmetic manually.

### Angles

- **Meteorological angles**: 0° = North, clockwise. Used for wind direction input.
- **Mathematical angles**: 0° = East, counter-clockwise. Used internally for geometry.
- Always convert explicitly — never assume which convention is in use.

```python
from hera.utils import toMeteorologicalAngle, toMathematicalAngle
```

### Coordinate Systems

Always use the named constants — never raw EPSG integers inline:

```python
from hera.measurements.GIS import WSG84, ITM, convertCRS
# WSG84 = 4326, ITM = 2039
```

### Caching

Use the Cache collection for expensive intermediate results. Check before computing:

```python
cached = proj.getCacheDocuments(type="MyResult", desc={"param": value})
if cached:
    result = cached[0].getData()
else:
    result = expensive_computation()
    proj.addCacheDocument(resource=result, dataFormat="JSON_dict",
                          type="MyResult", desc={"param": value})
```

---

## Testing Rules

- All tests live in `hera/tests/`.
- Use `pytest` with session-scoped fixtures. Load test data via a `test_repository.json`.
- Test class names: `Test<ToolkitName>`.
- Compare results with expected outputs using the `compare_outputs()` helper — do not assert equality of DataFrames manually.
- Never connect to production MongoDB in tests — use the test project name defined in `conftest.py`.
- Add expected output files to `hera/tests/expected/` when adding new tests.

```python
# conftest.py pattern
@pytest.fixture(scope="session")
def my_toolkit(hera_test_project):
    return MyToolkit(projectName=PYTEST_PROJECT_NAME)
```

---

## CLI Rules

- All CLI entry points go in `hera/bin/` and are registered in `setup.py`.
- CLI scripts use `argparse` or `click`, consistent with existing CLIs in the same domain.
- Command structure: `hera-<domain> <subcommand> <action> [options]`
- Always require `--projectName` (or `-p`) as an argument for any command that touches data.

---

## Common Pitfalls to Avoid

1. **Do not call `.compute()` on Dask DataFrames prematurely** — filter first, then compute.
2. **Do not use absolute file paths in repository JSON** — always use `isRelativePath: True`.
3. **Do not guess toolkit name strings** — use `toolkitHome.<CONSTANT>` constants.
4. **Do not store simulation results in the Measurements collection**.
5. **Do not add new built-in toolkits by editing the `_toolkits` dict directly** — the roadmap is to move to a registry JSON; prefer dynamic registration.
6. **Do not write angle conversion math inline** — use `hera.utils` helpers.
7. **Do not skip versioning** — every datasource must have a `[major, minor, patch]` version tuple.
8. **Do not subclass the toolkit for analysis/presentation logic** — use separate composition classes.

---

## Repository JSON Format (reference)

```json
{
    "MeteoLowFreq": {
        "Config": { "stationType": "IMS" },
        "Datasource": {
            "YAVNEEL": {
                "isRelativePath": "True",
                "item": {
                    "resource": "measurements/meteorology/YAVNEEL.parquet",
                    "dataFormat": "parquet"
                }
            }
        }
    }
}
```

---

## Key Entry Points (quick reference)

```python
from hera import Project, toolkitHome

# Open a project
proj = Project(projectName="MY_PROJECT")

# Get a toolkit
topo = toolkitHome.getToolkit(toolkitHome.GIS_RASTER_TOPOGRAPHY, projectName="MY_PROJECT")
meteo = toolkitHome.getToolkit(toolkitHome.METEOROLOGY_LOWFREQ, projectName="MY_PROJECT")
lsm = toolkitHome.getToolkit(toolkitHome.LSM, projectName="MY_PROJECT")
risk = toolkitHome.getToolkit(toolkitHome.RISKASSESSMENT, projectName="MY_PROJECT")

# Register a datasource
toolkit.addDataSource("NAME", "/path/to/file", "parquet", version=[0, 0, 1])
toolkit.setDataSourceDefaultVersion("NAME", [0, 0, 1])

# Load data
df = toolkit.getDataSourceData("NAME")

# Project config
proj.setConfig(key="value")
config = proj.getConfig()
```

---

## Workflow for Adding a New Toolkit

1. Create `hera/<domain>/<name>/toolkit.py` with a class inheriting `abstractToolkit`.
2. Create `hera/<domain>/<name>/analysis.py` with the analysis layer class.
3. Create `hera/<domain>/<name>/presentation.py` with the presentation layer class.
4. Add a CLI in `hera/<domain>/<name>/CLI.py` and register in `setup.py`.
5. Write tests in `hera/tests/test_<name>.py` following the session-fixture pattern.
6. Add expected output files in `hera/tests/expected/`.
7. Register the toolkit dynamically (or add to `toolkits_registry.json` per roadmap).
8. Document the toolkit following the existing MkDocs structure under `docs/toolkits/`.

---

## UI Development

- The UI client lives in `ui/client/` (React + TypeScript + Vite).
- **Before reporting any UI work as done, follow the full validation checklist in `ui/client/TEST_UI.md`.** Do not run tsc, tests, or build as individual ad-hoc commands — always follow the checklist in order.
