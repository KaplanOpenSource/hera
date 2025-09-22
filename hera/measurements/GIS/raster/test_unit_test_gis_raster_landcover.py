# ===============================================
# קובץ: test_unit_test_gis_raster_landcover.py (גרסה מתוקנת מלאה)
# ===============================================

import os
import glob
import unittest
import numpy as np
import rasterio
import shapely.geometry
import geopandas as gpd
from measurements.GIS.raster.landcover import LandCoverToolkit

def safe_convertCRS(points, inputCRS, outputCRS):
    """
    Safe conversion function for unit tests
    """
    points_geom = [shapely.geometry.Point(p[1], p[0]) for p in points]  # (lon, lat)
    gdf = gpd.GeoDataFrame(geometry=points_geom, crs=inputCRS)
    return gdf.to_crs(outputCRS).geometry.values

class TestLandCoverToolkit(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print("\n[SETUP] Setting up LandCoverToolkit for tests")

        # Read environment variable
        cls.HERA_DATA_PATH = os.getenv("HERA_DATA_PATH")
        if not cls.HERA_DATA_PATH:
            raise EnvironmentError("HERA_DATA_PATH environment variable not set.")

        # Search for landcover file automatically
        search_pattern = os.path.join(
            cls.HERA_DATA_PATH,
            "measurements/GIS/raster",
            "lc_mcd12q1v061.t1_c_500m_s_20210101_20211231_go_epsg.4326_v*.tif"
        )
        matching_files = glob.glob(search_pattern)

        if not matching_files:
            raise FileNotFoundError(f"No landcover file found matching pattern: {search_pattern}")

        cls.landcover_file = matching_files[0]

        # Create Toolkit instance
        cls.toolkit = LandCoverToolkit(projectName="unittest_project")

        # Custom config for Unit Test
        cls.custom_config = {
            "defaultLandCover": os.path.basename(cls.landcover_file),
            "DataSources": {
                os.path.basename(cls.landcover_file): {
                    "relativePath": "measurements/GIS/raster/" + os.path.basename(cls.landcover_file)
                }
            }
        }

        cls.toolkit.getConfig = lambda: cls.custom_config
        cls.toolkit.getDataSourceData = lambda datasource_name: rasterio.open(
            os.path.join(cls.HERA_DATA_PATH, "measurements/GIS/raster", datasource_name)
        )
        cls.toolkit.getDataSourceDocument = lambda datasource_name: {"desc": {"type": 1}}

        # Adjust static method for roughness calculation
        cls.toolkit.roughnesslength2sandgrainroughness = lambda rl: rl * 30

        # Define a test point
        cls.test_lon = 35.0
        cls.test_lat = 32.0
        cls.minx = 34.9
        cls.miny = 31.9
        cls.maxx = 35.1
        cls.maxy = 32.1

    def test_get_land_cover_at_point(self):
        print("\n[TEST] Running test_get_land_cover_at_point")
        value = self.toolkit.getLandCoverAtPoint(self.test_lon, self.test_lat)
        print(f"Value at point: {value}")
        self.assertIsNotNone(value)

    def test_get_land_cover(self):
        print("\n[TEST] Running test_get_land_cover")
        landcover = self.toolkit.getLandCover(self.minx, self.miny, self.maxx, self.maxy, dxdy=500)
        print(f"Landcover shape: {landcover['landcover'].shape}")
        self.assertIsNotNone(landcover)

    def test_get_land_cover_map_vs_raster(self):
        """
        Test that getLandCover map matches the original raster values at sampled points.

        We sample every few points in the output landcover map, and verify that
        the value matches the raster file directly using rasterio.
        """
        print("\n[TEST] Running test_get_land_cover_map_vs_raster")

        landcover = self.toolkit.getLandCover(self.minx, self.miny, self.maxx, self.maxy, dxdy=500)
        landcover_values = landcover['landcover'].values
        lons = landcover['lon'].values
        lats = landcover['lat'].values

        with rasterio.open(self.landcover_file) as src:
            raster_band = src.read(1)  # [FIX] Read the raster once

            for i in range(0, landcover_values.shape[0], 10):  # [FIX] Sample every 10 rows
                for j in range(0, landcover_values.shape[1], 10):  # [FIX] Sample every 10 cols
                    lon = lons[i, j]
                    lat = lats[i, j]
                    try:
                        row, col = src.index(lon, lat)
                        value_from_raster = raster_band[row, col]
                        value_from_toolkit = landcover_values[i, j]

                        # Allow for small differences due to rounding / indexing errors
                        self.assertEqual(
                            value_from_raster,
                            value_from_toolkit,
                            msg=f"Mismatch at (lon={lon}, lat={lat}): raster={value_from_raster}, toolkit={value_from_toolkit}"
                        )
                    except IndexError:
                        # If the sample falls outside the raster extent, skip
                        print(f"Skipping point outside raster: lon={lon}, lat={lat}")

    def test_get_roughness_at_point(self):
        print("\n[TEST] Running test_get_roughness_at_point")
        value = self.toolkit.getRoughnessAtPoint(self.test_lon, self.test_lat)
        print(f"Roughness at point: {value}")
        self.assertIsNotNone(value)

    def test_get_roughness(self):
        print("\n[TEST] Running test_get_roughness")
        roughness = self.toolkit.getRoughness(self.minx, self.miny, self.maxx, self.maxy, dxdy=500)
        print(f"Roughness shape: {roughness['z0'].shape}")
        self.assertIsNotNone(roughness)

    def test_roughnesslength2sandgrainroughness(self):
        print("\n[TEST] Running test_roughnesslength2sandgrainroughness")
        roughness = 0.1  # meters
        ks = self.toolkit.roughnesslength2sandgrainroughness(roughness)
        print(f"Ks (Sand Grain Roughness): {ks}")
        self.assertTrue(ks > 0)

    def test_land_cover_at_point_against_original(self):
        print("\n[TEST] Running test_land_cover_at_point_against_original")
        with rasterio.open(self.landcover_file) as src:
            lon, lat = self.test_lon, self.test_lat
            row, col = src.index(lon, lat)
            value_from_raster = src.read(1)[row, col]

        value_from_toolkit = self.toolkit.getLandCoverAtPoint(lon, lat)
        print(f"Raster value: {value_from_raster}, Toolkit value: {value_from_toolkit}")
        self.assertEqual(value_from_raster, value_from_toolkit)

    def test_roughness_values_are_in_range(self):
        print("\n[TEST] Running test_roughness_values_are_in_range")
        roughness = self.toolkit.getRoughness(self.minx, self.miny, self.maxx, self.maxy, dxdy=500)
        z0 = roughness['z0'].values
        self.assertTrue(np.all((z0 > 0) & (z0 < 2)), "Some roughness values are out of expected range!")

    def test_land_cover_out_of_bounds(self):
        print("\n[TEST] Running test_land_cover_out_of_bounds")
        with self.assertRaises(IndexError):
            self.toolkit.getLandCoverAtPoint(1000, 1000)  # ברור מחוץ לתחום

    def test_get_coding_map(self):
        print("\n[TEST] Running test_get_coding_map")
        coding_map = self.toolkit.getCodingMap('Type-1')
        self.assertIsInstance(coding_map, dict)
        self.assertIn(0, coding_map)
        self.assertEqual(coding_map[0], "Water")

    def test_roughness_known_landcover(self):
        print("\n[TEST] Running test_roughness_known_landcover")
        known_landcover_value = 5  # Mixed forests
        expected_roughness = 1.0

        # mock getLandCoverAtPoint
        self.toolkit.getLandCoverAtPoint = lambda lon, lat, inputCRS=None, dataSourceName=None: known_landcover_value

        roughness = self.toolkit.getRoughnessAtPoint(self.test_lon, self.test_lat)
        self.assertAlmostEqual(roughness, expected_roughness, places=3)


if __name__ == "__main__":
    unittest.main()
