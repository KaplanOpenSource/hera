
# 🛠 Hera Scripts

This folder contains utility scripts used in the Hera project.

## 📜 Script: `run_all_tests.sh`

This script runs **all unit tests** across the Hera system, including:

- Low-frequency meteorology toolkit
- High-frequency meteorology toolkit
- GIS raster toolkits (landcover, topography)
- GIS vector toolkit (demography)

### ▶️ How to Run

Navigate to the root project directory:

```bash
cd /home/ilay/hera/hera
make test
```

This will run the `run_all_tests.sh` script which executes all test modules.

> 💡 You can also run it manually with:
> `bash scripts/run_all_tests.sh`

---

## 📂 Environment Variable: `HERA_UNITTEST_DATA`

To allow unit tests to locate required input files (such as `.parquet` or `.hgt`), the system uses an environment variable named:

```
HERA_UNITTEST_DATA
```

This variable should point to the base path where your data directory is stored. For example:

```bash
export HERA_UNITTEST_DATA=/home/ilay/hera_unittest_data
```

> 🧠 This allows the tests and documentation to dynamically resolve file paths like:
>
> `$HERA_UNITTEST_DATA/measurements/meteorology/lowfreqdata/YAVNEEL.parquet`

Make sure that:
- You have the required data files under the expected folder structure.
- The environment variable is set in your shell (or in `.bashrc`, `.zshrc`, etc.) before running the tests.

---

## 🔁 Future Scripts

This folder is intended for future utility scripts as well — such as:
- Automatic data validation
- Report generation
- CI testing hooks
