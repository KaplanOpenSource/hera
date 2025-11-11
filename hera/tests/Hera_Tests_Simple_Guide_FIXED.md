# Hera JSON Tests — Clean Setup & Run (English)

This guide reflects the **new flow**: the shell script reads a **local `tests/env.template`**, no user-specific exports inside the script, and **named result sets** are mandatory.

---

## 0) Folder layout you should expect

```
<repo-root>/
  tests/
    run_all_json_tests.sh
    env.template              ← local, gitignored
    json_definitions/         ← your JSON test definitions
    expected/
      BASELINE/               ← example result set (created by prepare)
        ...
  run_all_definitions.py      ← (or tests/run_all_definitions.py)
```

> Either `run_all_definitions.py` lives at the repo root **or** under `tests/` — the script auto-detects both.

---

## 1) Create `tests/env.template` (edit to your machine)

Create the file `tests/env.template` with **English comments only**:

```bash
# =======================
# Local developer settings
# =======================

# ---- Python ----
# Path to your Python executable.
PYTHON_BIN=python3

# Add the repo root so its modules are importable.
# IMPORTANT: default expansion so set -u won't fail if PYTHONPATH is unset.
PYTHONPATH=/home/ilay/hera:${PYTHONPATH:-}

# ---- Data roots ----
# Root folder that contains the unit-test datasets (HGT, GeoTIFF, SHP, Parquet, etc.).
HERA_DATA_PATH=/home/ilay/hera_unittest_data

# Some tests read from this var as well. Keep equal to HERA_DATA_PATH unless you know otherwise.
HERA_UNITTEST_DATA=/home/ilay/hera_unittest_data

# ---- Test runner defaults (optional) ----
# Named result set for expected outputs. Can be overridden by --result-set.
RESULT_SET=BASELINE
```

---

## 2) One-time permission

```bash
chmod +x tests/run_all_json_tests.sh
```

---

## 3) Prepare mode (creates/updates expected results)

```bash
./tests/run_all_json_tests.sh prepare --result-set BASELINE
```

---

## 4) Run mode (compares current outputs vs expected)

```bash
./tests/run_all_json_tests.sh run --result-set BASELINE
```

---

## 5) Sanity checks

```bash
source tests/env.template
echo "$PYTHON_BIN"
echo "$PYTHONPATH"
echo "$HERA_DATA_PATH"
```

---

## 6) Common pitfalls  
“Missing tests/env.template” → create it (see step 1) and set your paths.

“Expected results set 'NAME' not found … run prepare first” → you tried run before prepare.

Module import errors → confirm PYTHONPATH in env.template points to your repo.

Data not found → fix HERA_DATA_PATH / HERA_UNITTEST_DATA to your test-data root.
---

## 7) FAQ  
Q: Where are expected files stored?
A: Under tests/expected/<RESULT_SET>/…. You can keep multiple named sets.

Q: Can I place run_all_definitions.py under tests/?
A: Yes. The shell script looks in both the repo root and tests/.
---

## 8) How this solves the review notes  
“Remove exports from run_all_json.sh”
The script contains no user-specific exports. It only sources tests/env.template.

“Script should only source a template file (env.template) which is gitignored”
We load tests/env.template, and you should add it to .gitignore.

“Preparation must support more than one result set and be mandatory parameter”
We added --result-set <NAME> (or RESULT_SET in env.template).
You can keep several sets (e.g., BASELINE, REGRESSION_*, etc.).

“Compare without result set name or if set does not exist must error and request to create it first”
run checks if tests/expected/<RESULT_SET>/ exists. If not, it prints a clear error telling you to run prepare first.

“Load/save consistency under subfolders (expected_outputs/…)”
The loader now resolves all relative paths into the active result-set directory, so files saved during prepare are found during run.