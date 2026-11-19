# Test fixtures

- `dots.fgb`: three WGS84 points with the title `dots max=2`, written by GDAL as `ui/server/geo_shapes.py` writes dots. The `flatgeobuf` package cannot write a title.
- `ZipSample.zip`: a project export from hera's `Project.export` (without progress bar). Project `ZipSample`, files folder `/projects/ZipSample`, one settings record and one data record.
- `ZipSampleBatched.zip`: the same project exported with the progress bar, so each zip entry holds a list of records.
