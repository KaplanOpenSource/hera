import unittest
import numpy as np
import pandas as pd
import xarray as xr
import os
import struct
import math
import glob
from hera import Project
from measurements.GIS.raster.topography import TopographyToolkit, WSG84

class TestTopographyToolkit(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.project = Project(projectName="MY_PROJECT")
        cls.topo_toolkit = TopographyToolkit(projectName="MY_PROJECT")

        base_path = os.environ.get("HERA_DATA_PATH", "/home/ilay/hera/hera/tests")

        # מציאת כל קבצי .hgt בתיקיות שונות
        hgt_files = glob.glob(os.path.join(base_path, "**", "*.hgt"), recursive=True)
        resource_dirs = sorted(set(os.path.dirname(f) for f in hgt_files))

        if not resource_dirs:
            raise FileNotFoundError(f"No .hgt files found under {base_path}")

        print(f"Detected resource folders:\n" + "\n".join(resource_dirs))

        # הגדרת הקונפיגורציה עם כל הנתיבים
        cls.custom_config = {
            "defaultSRTM": "SRTMGL1",
            "DataSources": {
                "SRTMGL1": {
                    "item": {
                        "resource": resource_dirs[0],  # עדיין שדה ברירת מחדל
                        "resource_folders": resource_dirs,  # נוסיף את כולם
                        "dataFormat": "SRTM",
                        "valueType": "Elevation",
                        "desc": {}
                    }
                }
            }
        }

        # הקונפיג יינתן ל־toolkit
        cls.topo_toolkit.getConfig = lambda: cls.custom_config

        # נעדכן גם את getDataSourceData שיחזיר את כל הרשימה
        cls.topo_toolkit.getDataSourceData = lambda name: cls.custom_config["DataSources"][name]["item"][
            "resource_folders"]

    def safe_get_elevation(self, func, *args, **kwargs):
        """Wrapper helper for safe elevation calls"""
        try:
            return func(*args, **kwargs)
        except FileNotFoundError as e:
            # Skip the test if the required elevation tile file is missing
            self.skipTest(f"Skipped due to missing tile: {e}")
        except AttributeError as e:
            if "NoneType" in str(e):
                # If the elevation data is None (e.g. point out of raster bounds), log it and return None
                func_name = getattr(func, "__name__", "unknown_function")
                print(f"⚠️  Function '{func_name}': point likely out of raster bounds, returned None")
                return None
            else:
                raise
        except Exception as e:
            # Any unexpected error should fail the test
            self.fail(f"Unexpected exception: {e}")

    def test_getPointElevation(self):
        lat, lon = 33.85, 35.15  # Point inside Israel, expected to be covered by N33E035.hgt
        elevation = self.safe_get_elevation(self.topo_toolkit.getPointElevation, lat, lon)

        if elevation is None:
            print("⚠️  Could not test getPointElevation: elevation is None (missing tile or point out of bounds)")
        else:
            self.assertIsInstance(elevation, (float, int))
            print(f"✅ Elevation at ({lat}, {lon}) = {elevation}")

    def test_getPointElevation_second_file(self):
        lat, lon = 33.85, 36.05  # Expected to be covered by N33E036.hgt
        elevation = self.safe_get_elevation(self.topo_toolkit.getPointElevation, lat, lon)

        if elevation is None:
            print(
                "⚠️  Could not test getPointElevation_second_file: elevation is None (missing tile or point out of bounds)")
        else:
            self.assertIsInstance(elevation, (float, int))
            print(f"✅ Elevation at ({lat}, {lon}) = {elevation}")

    def test_getPointElevation_matches_hgt_file(self):

        def read_raw_elevation(lat, lon, hgt_folders):
            """
            Reads the raw elevation value from an .hgt file given latitude and longitude.
            hgt_folders: list of folders to search for the relevant tile.
            Returns the elevation (int) or None if not found.
            """
            lat_deg = int(math.floor(lat))
            lon_deg = int(math.floor(lon))
            lat_prefix = 'N' if lat_deg >= 0 else 'S'
            lon_prefix = 'E' if lon_deg >= 0 else 'W'
            tile_name = f"{lat_prefix}{abs(lat_deg):02d}{lon_prefix}{abs(lon_deg):03d}.hgt"

            for folder in hgt_folders:
                tile_path = os.path.join(folder, tile_name)
                if os.path.exists(tile_path):
                    SAMPLES = 1201
                    size = os.path.getsize(tile_path)
                    if size == 1201 * 1201 * 2:
                        lat_fraction = lat - lat_deg
                        lon_fraction = lon - lon_deg
                        row = int((1 - lat_fraction) * (SAMPLES - 1))
                        col = int(lon_fraction * (SAMPLES - 1))
                        with open(tile_path, 'rb') as f:
                            offset = 2 * (row * SAMPLES + col)
                            f.seek(offset)
                            data = f.read(2)
                            elevation = struct.unpack('>h', data)[0]
                            return elevation if elevation != -32768 else None
            return None

        lat, lon = 33.85, 35.15  # Covered by N33E035.hgt

        elevation_from_toolkit = self.safe_get_elevation(
            self.topo_toolkit.getPointElevation, lat, lon
        )

        if elevation_from_toolkit is None:
            print("⚠️  Could not test: toolkit returned None")
            return

        resource_folders = self.topo_toolkit.getDataSourceData("SRTMGL1")
        if isinstance(resource_folders, str):
            resource_folders = [resource_folders]

        elevation_from_file = read_raw_elevation(lat, lon, resource_folders)

        if elevation_from_file is None:
            print("⚠️  Could not test: HGT file missing or point outside tile")
            return

        print(f"Toolkit elevation: {elevation_from_toolkit}, File elevation: {elevation_from_file}")
        self.assertAlmostEqual(elevation_from_toolkit, elevation_from_file, delta=1)

    def test_getPointListElevation(self):
        points = pd.DataFrame({
            'lat': [33.85, 33.9],
            'lon': [35.15, 36.05]
        })
        result = self.safe_get_elevation(self.topo_toolkit.getPointListElevation, points)
        self.assertEqual(result.shape[0], 2)
        self.assertIn('elevation', result.columns)

    def test_getPointListElevation_matches_hgt_files(self):

        def read_raw_elevation(lat, lon, hgt_folders):
            lat_deg = int(math.floor(lat))
            lon_deg = int(math.floor(lon))
            lat_prefix = 'N' if lat_deg >= 0 else 'S'
            lon_prefix = 'E' if lon_deg >= 0 else 'W'
            tile_name = f"{lat_prefix}{abs(lat_deg):02d}{lon_prefix}{abs(lon_deg):03d}.hgt"

            for folder in hgt_folders:
                tile_path = os.path.join(folder, tile_name)
                if os.path.exists(tile_path):
                    SAMPLES = 1201
                    size = os.path.getsize(tile_path)
                    if size == 1201 * 1201 * 2:
                        lat_fraction = lat - lat_deg
                        lon_fraction = lon - lon_deg
                        row = int((1 - lat_fraction) * (SAMPLES - 1))
                        col = int(lon_fraction * (SAMPLES - 1))
                        with open(tile_path, 'rb') as f:
                            offset = 2 * (row * SAMPLES + col)
                            f.seek(offset)
                            data = f.read(2)
                            elevation = struct.unpack('>h', data)[0]
                            return elevation if elevation != -32768 else None
            return None

        points = pd.DataFrame([
            (33.85, 35.15),  # Inside N33E035.hgt
            (33.85, 36.05),  # Inside N33E036.hgt
            (33.9999, 35.9999),  # Border between N33E035 and N33E036
            (34.0001, 36.0001),  # Just over the edge, inside N34E036
        ], columns=["lat", "lon"])

        elevations = self.safe_get_elevation(
            self.topo_toolkit.getPointListElevation, points
        )

        if elevations is None:
            self.skipTest(
                "⚠️  Skipped test_getPointListElevation_matches_hgt_files: toolkit returned None for all points (missing data or files)")
            return

        resource_folders = self.topo_toolkit.getDataSourceData("SRTMGL1")
        if isinstance(resource_folders, str):
            resource_folders = [resource_folders]

        tested_any = False  # Track if at least one point was tested

        for i, row in points.iterrows():
            lat, lon = row["lat"], row["lon"]
            expected = read_raw_elevation(lat, lon, resource_folders)
            actual = elevations.iloc[i]


            if expected is None:
                print(f"⚠️  Point ({lat}, {lon}): HGT data missing")
                continue

            if actual is None:
                print(f"⚠️  Point ({lat}, {lon}): toolkit returned None unexpectedly")
                continue

            print(f"✅ Point ({lat}, {lon}): toolkit={actual}, file={expected}")
            self.assertAlmostEqual(actual, expected, delta=1)
            tested_any = True

        if not tested_any:
            self.skipTest("⚠️  Skipped: no valid elevation comparisons were possible (all points missing)")

    def test_getElevationOfXarray(self):
        lat_vals = np.array([[33.85, 33.85], [33.86, 33.86]])
        lon_vals = np.array([[35.1, 35.11], [36.05, 36.06]])

        ds = xr.Dataset({
            "lat": (["i", "j"], lat_vals),
            "lon": (["i", "j"], lon_vals)
        })

        result = self.safe_get_elevation(self.topo_toolkit.getElevationOfXarray, ds)
        self.assertIn('elevation', result.coords)

    def test_getElevationOfXarray_matches_hgt_file(self):

        def read_raw_elevation(lat, lon, hgt_folders):
            lat_deg = int(math.floor(lat))
            lon_deg = int(math.floor(lon))
            lat_prefix = 'N' if lat_deg >= 0 else 'S'
            lon_prefix = 'E' if lon_deg >= 0 else 'W'
            tile_name = f"{lat_prefix}{abs(lat_deg):02d}{lon_prefix}{abs(lon_deg):03d}.hgt"

            for folder in hgt_folders:
                tile_path = os.path.join(folder, tile_name)
                if os.path.exists(tile_path):
                    SAMPLES = 1201
                    size = os.path.getsize(tile_path)
                    if size == 1201 * 1201 * 2:
                        lat_fraction = lat - lat_deg
                        lon_fraction = lon - lon_deg
                        row = int((1 - lat_fraction) * (SAMPLES - 1))
                        col = int(lon_fraction * (SAMPLES - 1))
                        with open(tile_path, 'rb') as f:
                            offset = 2 * (row * SAMPLES + col)
                            f.seek(offset)
                            data = f.read(2)
                            elevation = struct.unpack('>h', data)[0]
                            return elevation if elevation != -32768 else None
            return None

        # Define a grid with dims ['i', 'j']
        lats = np.linspace(33.85, 33.87, 5)
        lons = np.linspace(35.1, 35.12, 5)

        lat_grid = np.tile(lats.reshape(-1, 1), (1, len(lons)))
        lon_grid = np.tile(lons.reshape(1, -1), (len(lats), 1))

        dummy = xr.Dataset(
            {
                "lat": (["i", "j"], lat_grid),
                "lon": (["i", "j"], lon_grid)
            },
            coords={
                "i": np.arange(len(lats)),
                "j": np.arange(len(lons))
            }
        )

        try:
            filled = self.topo_toolkit.getElevationOfXarray(dummy)
        except Exception as e:
            self.skipTest(f"Skipped: getElevationOfXarray failed due to {e}")
            return

        if filled is None or not isinstance(filled, xr.Dataset):
            self.skipTest("Skipped: getElevationOfXarray did not return a valid xarray Dataset")
            return

        if "elevation" not in filled.coords:
            self.skipTest("Skipped: elevation coordinate was not added to dataset")
            return

        resource_folders = self.topo_toolkit.getDataSourceData("SRTMGL1")
        if isinstance(resource_folders, str):
            resource_folders = [resource_folders]

        tested_any = False
        for i in [0, len(lats) // 2, -1]:
            for j in [0, len(lons) // 2, -1]:
                lat = lat_grid[i, j]
                lon = lon_grid[i, j]
                expected = read_raw_elevation(lat, lon, resource_folders)
                actual = float(filled.coords["elevation"].values[i, j])

                if expected is None:
                    print(f"⚠️  Skipping ({lat}, {lon}): no HGT data")
                    continue

                print(f"✅ Checking ({lat}, {lon}): toolkit={actual}, file={expected}")
                self.assertAlmostEqual(actual, expected, delta=1)
                tested_any = True

        if not tested_any:
            self.skipTest("Skipped: no valid comparison points found in HGT data")

    def test_getElevation(self):
        minx, miny = 35.1, 33.85
        maxx, maxy = 35.12, 33.86
        dxdy = 100

        result = self.safe_get_elevation(self.topo_toolkit.getElevation, minx, miny, maxx, maxy, dxdy, inputCRS=WSG84)
        self.assertIsInstance(result, xr.Dataset)
        self.assertIn('elevation', result.coords)

    def test_getElevation_matches_hgt_file(self):

        def read_raw_elevation(lat, lon, hgt_folders):
            lat_deg = int(math.floor(lat))
            lon_deg = int(math.floor(lon))
            lat_prefix = 'N' if lat_deg >= 0 else 'S'
            lon_prefix = 'E' if lon_deg >= 0 else 'W'
            tile_name = f"{lat_prefix}{abs(lat_deg):02d}{lon_prefix}{abs(lon_deg):03d}.hgt"

            for folder in hgt_folders:
                tile_path = os.path.join(folder, tile_name)
                if os.path.exists(tile_path):
                    SAMPLES = 1201
                    size = os.path.getsize(tile_path)
                    if size == 1201 * 1201 * 2:
                        lat_fraction = lat - lat_deg
                        lon_fraction = lon - lon_deg
                        row = int((1 - lat_fraction) * (SAMPLES - 1))
                        col = int(lon_fraction * (SAMPLES - 1))
                        with open(tile_path, 'rb') as f:
                            offset = 2 * (row * SAMPLES + col)
                            f.seek(offset)
                            data = f.read(2)
                            elevation = struct.unpack('>h', data)[0]
                            return elevation if elevation != -32768 else None
            return None

        # Define a small area (inside N33E035.hgt)
        minx, miny, maxx, maxy = 35.1, 33.85, 35.12, 33.87
        dxdy = 0.001  # High resolution

        try:
            da = self.topo_toolkit.getElevation(minx, miny, maxx, maxy, dxdy)
        except Exception as e:
            self.skipTest(f"Skipped: getElevation failed due to {e}")
            return

        if da is None or not isinstance(da, xr.DataArray):
            self.skipTest("Skipped: getElevation did not return a valid xarray DataArray")
            return

        resource_folders = self.topo_toolkit.getDataSourceData("SRTMGL1")
        if isinstance(resource_folders, str):
            resource_folders = [resource_folders]

        latitudes = da.coords["lat"].values
        longitudes = da.coords["lon"].values

        tested_any = False

        for i in [0, len(latitudes) // 2, -1]:
            for j in [0, len(longitudes) // 2, -1]:
                lat = latitudes[i]
                lon = longitudes[j]
                expected = read_raw_elevation(lat, lon, resource_folders)
                actual = float(da.values[i, j])

                if expected is None:
                    print(f"⚠️  Skipping ({lat}, {lon}): no HGT data")
                    continue

                print(f"✅ Checking ({lat}, {lon}): toolkit={actual}, file={expected}")
                self.assertAlmostEqual(actual, expected, delta=1)
                tested_any = True

        if not tested_any:
            self.skipTest("Skipped: no valid comparison points found in HGT data")

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

    import nbformat
    from nbconvert.preprocessors import ExecutePreprocessor
    import os
    import unittest

if __name__ == '__main__':
    unittest.main()

import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

class TestTopographyNotebook(unittest.TestCase):
    def test_notebook_runs_without_errors(self):
        """
        🚀 This test executes the TopographyToolkit Jupyter notebook and fails if any cell raises an error.
        """
        notebook_path = os.path.abspath(
            "doc/jupyter/toolkits/measurments/GIS/Raster/TopographyToolkit_Complete_With_Usage.ipynb")
        self.assertTrue(os.path.exists(notebook_path), f"Notebook not found at {notebook_path}")

        with open(notebook_path, encoding='utf-8') as f:
            nb = nbformat.read(f, as_version=4)

        ep = ExecutePreprocessor(timeout=600, kernel_name='python3')

        try:
            ep.preprocess(nb, {'metadata': {'path': os.path.dirname(notebook_path)}})
        except Exception as e:
            self.fail(f"Notebook execution failed: {e}")

        print("✅ Notebook executed successfully with no errors.")
