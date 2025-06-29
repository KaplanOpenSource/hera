# Hera Toolkit - JSON Test Runner

This package includes a test framework for running dynamic integration tests on Hera toolkits using JSON definition files.

## 🚀 How It Works

Each JSON file defines a test case:
- Class path (e.g., `hera.measurements.GIS.raster.landcover.LandCoverToolkit`)
- Method to call (e.g., `getLandCoverAtPoint`)
- Parameters (raw values, env-based paths, or references to previous outputs)
- Expected output type and path

The runner `run_all_definitions.py` dynamically imports classes and runs the tests.

## ▶️ Running the Tests

```bash
cd hera/tests
chmod +x run_all_json_tests.sh
./run_all_json_tests.sh
```

Make sure to export these before running:

```bash
export HERA_UNITTEST_DATA=/home/ilay/hera_unittest_data
export HERA_DATA_PATH=$HERA_UNITTEST_DATA
export PREPARE_EXPECTED_OUTPUT=1  # First run: generate outputs
```

Set `PREPARE_EXPECTED_OUTPUT=0` to validate against saved outputs.

## 📁 Folder Structure

```
hera/
├── tests/
│   ├── run_all_definitions.py
│   ├── run_all_json_tests.sh
│   └── json_definitions/
│       ├── *.json
├── expected_outputs/
```

## 💾 Output Files

- `.json` – for basic types (float, int, dict)
- `.geojson` – for GeoDataFrame
- `.nc` – for xarray

# Hera Toolkit JSON Test Runner

This project runs a full suite of functional tests for Hera toolkits using JSON definition files.

## 🧪 What does this do?

This script executes all defined method tests from the following toolkits:

- 🗺️ TopographyToolkit
- 🌦️ LowFreqToolKit
- 👨‍👩‍👧 DemographyToolkit
- 🌡️ HighFreqToolKit
- 🌍 LandCoverToolkit

Each JSON file contains a series of test definitions, including method name, parameters, and expected output.

## 🚀 How to run the tests

Make sure you are in your virtual environment (`heraenv`) and inside the `hera/tests` directory.

Then run:

```bash
chmod +x run_all_json_tests.sh
./run_all_json_tests.sh
