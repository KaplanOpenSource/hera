# Test Data for GIS Raster Topography Toolkit

This folder contains test data for running unit tests on the `TopographyToolkit`.

## Files
- `N33E035.hgt`: NASA SRTM elevation data for a tile covering latitude 33°N and longitude 35°E.

## Why is this needed?
The unit tests (`test_unit_test_gis_raster_topography.py`) need a real HGT file to check that the elevation calculations work correctly.

Without this file, tests will fail because they cannot open the raster dataset.

## How to use
- Make sure this file `N33E035.hgt` exists inside the `hera/tests/test_data/` folder.
- If missing, you can download it manually from NASA's SRTM data sources.
- Alternatively, you can place any `.hgt` file here but **update the test code accordingly**.

## Important
- The tests automatically set the environment variable `HERA_TEST_DATAPATH` to point to this folder.
- No need to manually configure anything if the file exists in the right place.
