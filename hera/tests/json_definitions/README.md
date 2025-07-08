# 📦 Hera Toolkit – JSON Test Runner

This module provides a flexible framework for running functional tests on Hera toolkits using JSON definition files.

---

## 🚀 How It Works

Each test is defined in a `.json` file and includes:

- `class_path`: Full import path of the toolkit class  
  _Example_: `"hera.measurements.GIS.raster.landcover.LandCoverToolkit"`

- `method_name`: Name of method (or attribute) to test  
  _Example_: `"getLandCoverAtPoint"`

- `parameters`: Dictionary of parameters to pass.  
  These can be:
  - Raw values (`int`, `str`, etc.)
  - Paths based on environment variables (`fromEnv`)
  - References to previous outputs (`fromFunction`)

- `output_type`: Type of expected output (`dataframe`, `geojson`, `ndarray`, etc.)

- `output_filename`: Path to expected output file, relative to `expected_outputs/`

- `postprocess` (optional): Additional steps to apply on the result (e.g. `getData`, `result`, `rawData`)

---

## ▶️ Running the Tests

Make sure you're inside the `hera/tests/` directory and working in the virtual environment (`heraenv`).

### 🧪 To compare current outputs to saved expected outputs:

```bash
export HERA_UNITTEST_DATA=/home/ilay/hera_unittest_data
export HERA_DATA_PATH=$HERA_UNITTEST_DATA
export PREPARE_EXPECTED_OUTPUT=0

./run_all_json_tests.sh
```

**What it does:**
- Loads each `.json` test definition
- Resolves parameters
- Executes the method on the toolkit
- Applies any post-processing
- Loads the saved output file (e.g., `.json`, `.geojson`, `.nc`)
- Compares the result to the saved output
- Prints a colored test summary

---

### 💾 To generate and save expected outputs (first time only):

```bash
export HERA_UNITTEST_DATA=/home/ilay/hera_unittest_data
export HERA_DATA_PATH=$HERA_UNITTEST_DATA
export PREPARE_EXPECTED_OUTPUT=1

./run_all_json_tests.sh
```

**What it does:**
- Runs all test methods just like in comparison mode
- Saves the result of each test into the appropriate file in `expected_outputs/`
- File format is based on `output_type`: `.json`, `.parquet`, `.geojson`, `.npz`, etc.

---

## 📁 Folder Structure

```
hera/
├── tests/
│   ├── run_all_definitions.py        # Python test runner
│   ├── run_all_json_tests.sh         # Shell script wrapper
│   └── json_definitions/             # JSON test definition files
│       ├── test_definitions_topography.json
│       ├── test_definitions_lowfreq.json
│       └── ...
├── expected_outputs/                 # Saved expected results
```

---

## 🔍 Supported Toolkits

The framework supports dynamic testing of any Hera toolkit, including:

- 🗺️ `TopographyToolkit`
- 🌍 `LandCoverToolkit`
- 🌦️ `LowFreqToolKit`
- 🌡️ `HighFreqToolKit`
- 👨‍👩‍👧 `DemographyToolkit`

Each toolkit has a `.json` file in `json_definitions/` with method tests.

---

## ⚙️ Advanced Features

- ✅ **Environment paths**: Load files using `fromEnv` + `relative` structure.
- 🔁 **Post-processing**: Automatically run `.getData()`, `.result()`, etc.
- 🔄 **Cross-test reuse**: You can reuse output of previous functions using `"fromFunction": "previous_method_name"`.
- 🧠 **Smart comparison**: Uses tolerant comparison for DataFrames, arrays, floats, and nested types.
- 💥 **Exception logging**: Failing tests show type mismatches, structure differences, and sample values.

---

## 💡 Example: Test Definition

```json
{
  "class_path": "hera.measurements.meteorology.lowfreq.toolkit",
  "method_name": "getDailyStatistics",
  "init_parameters": {
    "dataSourceName": {
      "fromEnv": "HERA_UNITTEST_DATA",
      "relative": "measurements/meteorology/lowfreqdata/YAVNEEL.parquet"
    }
  },
  "parameters": {
    "stat": "mean",
    "var": "TG"
  },
  "output_type": "dataframe",
  "output_filename": "lowfreq/getDailyStatistics_mean_TG.parquet"
}
```

---

## ✅ Summary Output

After running the tests, you'll see a color-coded summary like:

```
✅ getDailyStatistics → matched expected output
❌ getElevation → mismatch in result type or content
```

```
Total: 10 | Passed: 8 | Failed: 2
```

---

📬 For any questions or suggestions, contact the Hera dev team.
