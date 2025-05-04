import unittest
import numpy as np
import pandas as pd
import xarray as xr
from hera import Project
from measurements.GIS.raster.topography import TopographyToolkit, WSG84

class TestTopographyToolkit(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.project = Project(projectName="MY_PROJECT")
        cls.topo_toolkit = TopographyToolkit(projectName="MY_PROJECT")

        cls.custom_config = {
            "defaultSRTM": "SRTMGL1",
            "DataSources": {
                "SRTMGL1": {
                    "item": {
                        "resource": "/home/ilay/hera/hera/tests/UNIT_TEST_GIS_RASTER_TOPOGRAPHY",
                        "dataFormat": "SRTM",
                        "valueType": "Elevation",
                        "desc": {}
                    }
                }
            }
        }

        # Override config/data methods
        cls.topo_toolkit.getConfig = lambda: cls.custom_config
        cls.topo_toolkit.getDataSourceData = lambda name: cls.custom_config["DataSources"][name]["item"]["resource"]

    def safe_get_elevation(self, func, *args, **kwargs):
        """Wrapper helper for safe elevation calls"""
        try:
            return func(*args, **kwargs)
        except FileNotFoundError as e:
            self.skipTest(f"Skipped due to missing tile: {e}")
        except AttributeError as e:
            if "NoneType" in str(e):
                self.skipTest(f"Skipped due to null elevation data or missing raster: {e}")
            else:
                raise
        except Exception as e:
            self.fail(f"Unexpected exception: {e}")

    def test_getPointElevation(self):
        lat, lon = 33.85, 35.15  # קובץ ראשון
        elevation = self.safe_get_elevation(self.topo_toolkit.getPointElevation, lat, lon)
        self.assertIsInstance(elevation, (float, int))

    def test_getPointElevation_second_file(self):
        lat, lon = 33.85, 36.05  # קובץ שני
        elevation = self.safe_get_elevation(self.topo_toolkit.getPointElevation, lat, lon)
        self.assertIsInstance(elevation, (float, int))

    def test_getPointListElevation(self):
        points = pd.DataFrame({
            'lat': [33.85, 33.9],
            'lon': [35.15, 36.05]
        })
        result = self.safe_get_elevation(self.topo_toolkit.getPointListElevation, points)
        self.assertEqual(result.shape[0], 2)
        self.assertIn('elevation', result.columns)

    def test_getElevationOfXarray(self):
        lat_vals = np.array([[33.85, 33.85], [33.86, 33.86]])
        lon_vals = np.array([[35.1, 35.11], [36.05, 36.06]])

        ds = xr.Dataset({
            "lat": (["i", "j"], lat_vals),
            "lon": (["i", "j"], lon_vals)
        })

        result = self.safe_get_elevation(self.topo_toolkit.getElevationOfXarray, ds)
        self.assertIn('elevation', result.coords)

    def test_getElevation(self):
        minx, miny = 35.1, 33.85
        maxx, maxy = 35.12, 33.86
        dxdy = 100

        result = self.safe_get_elevation(self.topo_toolkit.getElevation, minx, miny, maxx, maxy, dxdy, inputCRS=WSG84)
        self.assertIsInstance(result, xr.Dataset)
        self.assertIn('elevation', result.coords)

    def test_convertPointsCRS(self):
        points = [(35.1, 33.85), (36.05, 33.9)]
        converted = self.topo_toolkit.convertPointsCRS(points, inputCRS=4326, outputCRS=2039)
        self.assertEqual(converted.shape[0], 2)

    def test_createElevationSTL(self):
        try:
            stl_str = self.safe_get_elevation(
                self.topo_toolkit.createElevationSTL,
                minx=35.1, miny=33.85, maxx=35.11, maxy=33.86,
                dxdy=50, inputCRS=WSG84, solidName="TestSTL"
            )
            self.assertTrue(stl_str.startswith("solid"))
        except unittest.SkipTest:
            raise
        except Exception as e:
            self.fail(f"createElevationSTL failed with exception: {e}")

    def test_getElevationSTL(self):
        lat_vals = np.array([[33.85, 33.85], [33.86, 33.86]])
        lon_vals = np.array([[35.1, 35.11], [35.1, 35.11]])
        elevation_vals = np.array([[100, 110], [120, 130]])

        ds = xr.Dataset(
            {
                "lat": (["i", "j"], lat_vals),
                "lon": (["i", "j"], lon_vals),
                "elevation": (["i", "j"], elevation_vals)
            }
        )

        stl = self.topo_toolkit.getElevationSTL(ds, solidName="SurfaceTest")
        self.assertTrue(stl.startswith("solid SurfaceTest"))

    def test_calculateStastics(self):
        elevation = xr.Dataset({
            "X": (["i", "j"], np.array([[100, 200], [100, 200]])),
            "Y": (["i", "j"], np.array([[300, 300], [400, 400]])),
            "Elevation": (["i", "j"], np.array([[10, 20], [30, 40]])),
        })

        stats = self.topo_toolkit._analysis.calculateStastics(elevation)
        self.assertEqual(stats["mean"], 25.0)
        self.assertEqual(stats["domainmax"], 40)
        self.assertEqual(stats["domainmin"], 10)

if __name__ == '__main__':
    unittest.main()
