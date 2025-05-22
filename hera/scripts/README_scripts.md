
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

## 🔁 Future Scripts

This folder is intended for future utility scripts as well — such as:
- Automatic data validation
- Report generation
- CI testing hooks
