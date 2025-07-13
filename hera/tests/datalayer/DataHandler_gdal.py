from osgeo import gdal
import numpy as np
from hera.toolkits.gis.topography import TopographyToolkit


file_path = "/home/ilay/hera/hera/tests/UNIT_TEST_GIS_RASTER_TOPOGRAPHY/N33E035.hgt"

ds = gdal.Open(file_path)
if ds is None:
    raise FileNotFoundError(f"GDAL failed to open file: {file_path}")

band = ds.GetRasterBand(1)
elevation_array = band.ReadAsArray()

print("Shape:", elevation_array.shape)
print("Min elevation:", np.min(elevation_array))
print("Max elevation:", np.max(elevation_array))
print(toolkit.getDataSourceList())
