
# Dynamic Loading Test Suite — README

## Overview

This test suite validates Hera's dynamic loading system for two major components:

1. **Dynamic Experiment Loading** — ability to load experiment repositories, retrieve data, and interact with experiments via CLI.
2. **Dynamic Toolkit Loading** — ability to load toolkit Python classes dynamically at runtime using `classpath`.

These tests ensure that both mechanisms work in isolation (unit tests) and through the CLI & data layer (integration tests).

---

## Test Files and Their Purpose

### ✔ test_toolkit_dynamic_loading_unit.py
**Type:** Unit Test  
**Purpose:**
- Validates the _core mechanism_ of dynamic class loading.
- Creates a temporary Python package.
- Uses Hera’s `getData(..., desc={'classpath': ...})`.
- Ensures:
  - Class import works.
  - Instantiation succeeds.
  - No CLI / DB is required.

---

### ✔ test_experiment_cli_shortcuts.py
**Type:** Integration Test  
**Purpose:**
Validates CLI shortcuts:

- `hera-experiment list`
- `hera-experiment table`
- `hera-experiment data`

The test loads a dummy experiment and checks:

- `list` prints the experiment name.
- `table` prints datasource fields (dataFormat, toolkit, etc.)
- `data` prints actual parquet content and device metadata.

This ensures the CLI works end-to-end.

---

### ✔ test_toolkit_dynamic_loading_integration.py
**Type:** Integration Test  
**Purpose:**
- Tests full dynamic toolkit loading.
- Loads `DummyToolkit`.
- Verifies:
  - ToolkitDataSource is registered correctly.
  - `DataHandler.getData()` dynamically instantiates the toolkit.
  - Toolkit responds correctly (via `ping()`).

This confirms dynamic toolkit loading is working at runtime.

---

### ✔ Skipped Tests

One helper fixture (`load_dummy_experiment_to_project`) prepares a project & loads the dummy experiment.  
It is not executed alone.  
It is used by the integration tests.

Skipped tests occur because:
- They are fixtures or helpers, not stand-alone tests.
- Or marked to run only under specific markers (`unit`, `integration`).

---

## How to Run the Tests

### Run ALL tests:
```
pytest -q hera/tests/dynamic_loading_tests_pack
```

### Run only UNIT tests:
```
pytest -q hera/tests/dynamic_loading_tests_pack -m unit
```

### Run only INTEGRATION tests:
```
pytest -q hera/tests/dynamic_loading_tests_pack -m integration
```

---

## What We Achieved

- Built a reproducible DummyExperiment with:
  - repository file
  - runtime configuration
  - Sonic parquet
- Implemented CLI validation.
- Implemented dynamic toolkit validation.
- Ensured correct classpath loading.
- Ensured fully automated testing with no manual DB or CLI preparation.

This delivers a comprehensive baseline for anyone extending:
- experiment loaders
- toolkit loaders
- repository mapping logic

---


