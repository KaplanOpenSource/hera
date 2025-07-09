import unittest
import os
import warnings
from shapely.geometry import Polygon
from hera.measurements.GIS.vector.demography import DemographyToolkit
import geopandas as gpd
from hera.toolkit import TOOLKIT_SAVEMODE_NOSAVE
from hera.datalayer.document.metadataDocument import nonDBMetadataFrame
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

# ✅ Suppress deprecation warnings from Shapely/NumPy etc.
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", message="An exception was ignored while fetching the attribute .*__array_interface__.*")

class TestDemographyToolkit(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        base_path = os.environ.get("HERA_DATA_PATH", "/home/ilay/hera/hera/tests")
        cls.population_path = os.path.join(base_path, "measurements", "GIS", "vector", "population_lamas.shp")
        cls.population_gdf = gpd.read_file(cls.population_path)

        cls.toolkit = DemographyToolkit(projectName="UNIT_TEST_DEMOGRAPHY")
        cls.toolkit.loadData("lamas_population", cls.population_path, overwrite=True)

    def test_calculatePopulationInPolygon_basic(self):
        base_polygon = self.population_gdf.iloc[0].geometry
        minx, miny, maxx, maxy = base_polygon.bounds
        buffer = 10
        test_polygon = Polygon([
            (minx - buffer, miny - buffer),
            (maxx + buffer, miny - buffer),
            (maxx + buffer, maxy + buffer),
            (minx - buffer, maxy + buffer),
            (minx - buffer, miny - buffer)
        ])

        result = self.toolkit.analysis.calculatePopulationInPolygon(
            shapelyPolygon=test_polygon,
            dataSourceOrData="lamas_population"
        )

        print(f"🧺as Number of intersecting features: {len(result)}")
        if result.empty:
            print("❌ No population data found for the test polygon.")
        else:
            print(result[["geometry", "areaFraction", "total_pop"]].head())

        self.assertFalse(result.empty, "Expected non-empty result from polygon intersection")
        self.assertIn("geometry", result.columns)
        self.assertIn("areaFraction", result.columns)

    def test_calculatePopulationInPolygon_partial_intersection(self):
        polygon1 = self.population_gdf.iloc[0].geometry
        polygon2 = self.population_gdf.iloc[1].geometry

        intersection_area = polygon1.union(polygon2).centroid.buffer(20)

        result = self.toolkit.analysis.calculatePopulationInPolygon(
            shapelyPolygon=intersection_area,
            dataSourceOrData="lamas_population"
        )

        print(f"🧺as Partial intersection result rows: {len(result)}")
        self.assertFalse(result.empty, "Expected result for partially intersecting polygon")
        self.assertGreaterEqual(len(result), 1, "Expected at least one intersecting region")
        self.assertIn("areaFraction", result.columns)

    def test_calculatePopulationInPolygon_outside_bounds(self):
        minx, miny, maxx, maxy = self.population_gdf.total_bounds
        far_polygon = Polygon([
            (maxx + 10000, maxy + 10000),
            (maxx + 10100, maxy + 10000),
            (maxx + 10100, maxy + 10100),
            (maxx + 10000, maxy + 10100),
            (maxx + 10000, maxy + 10000)
        ])

        result = self.toolkit.analysis.calculatePopulationInPolygon(
            shapelyPolygon=far_polygon,
            dataSourceOrData="lamas_population"
        )

        print(f"🧺as Outside bounds result rows: {len(result)}")
        self.assertTrue(result.empty, "Expected empty result for polygon outside data bounds")

    def test_calculatePopulationInPolygon_invalid_datasource(self):
        polygon = self.population_gdf.iloc[0].geometry

        with self.assertRaises(ValueError) as cm:
            self.toolkit.analysis.calculatePopulationInPolygon(
                shapelyPolygon=polygon,
                dataSourceOrData="non_existing_data_source"
            )

        print(f"🛑 Caught expected ValueError: {cm.exception}")

    def test_createNewArea_simple(self):
        print("\n🚀 Running test_createNewArea_simple (Testing createNewArea function)")

        bounds = self.population_gdf.total_bounds
        buffer = 10
        polygon = Polygon([
            (bounds[0] - buffer, bounds[1] - buffer),
            (bounds[2] + buffer, bounds[1] - buffer),
            (bounds[2] + buffer, bounds[3] + buffer),
            (bounds[0] - buffer, bounds[3] + buffer),
            (bounds[0] - buffer, bounds[1] - buffer)
        ])

        if not polygon.is_valid:
            polygon = polygon.buffer(0)

        result = self.toolkit.analysis.createNewArea(
            shapeNameOrData=polygon,
            dataSourceOrData="lamas_population",
            saveMode=TOOLKIT_SAVEMODE_NOSAVE
        )

        self.assertIsInstance(result, nonDBMetadataFrame)

        gdf = result.getData()
        self.assertIsInstance(gdf, gpd.GeoDataFrame)
        self.assertEqual(len(gdf), 1)
        self.assertIn("geometry", gdf.columns)
        self.assertIn("total_pop", gdf.columns)

        expected = self.population_gdf["total_pop"].sum()
        actual = gdf.iloc[0]["total_pop"]
        print(f"✅ Expected: {expected}, Got: {actual}")
        self.assertAlmostEqual(actual, expected, places=2)

    def test_setDefaultDirectory_creates_and_sets_path(self):
        print("🚀 Running test_setDefaultDirectory_creates_and_sets_path")
        import tempfile

        with tempfile.TemporaryDirectory() as tmpdirname:
            test_dir = os.path.join(tmpdirname, "demography_test_folder")
            self.toolkit.setDefaultDirectory(test_dir, create=True)

            saved_path = self.toolkit.filesDirectory
            self.assertTrue(os.path.exists(test_dir), "📁 Directory was not created")
            self.assertEqual(os.path.abspath(test_dir), saved_path, "🗭 Saved path in filesDirectory doesn't match")

    def test_calculatePopulationInPolygon_with_known_values(self):
        print("🚀 Running test_calculatePopulationInPolygon_with_known_values")

        geometry = [
            Polygon([(0, 0), (2, 0), (2, 2), (0, 2)]),
            Polygon([(1, 1), (3, 1), (3, 3), (1, 3)]),
            Polygon([(5, 5), (6, 5), (6, 6), (5, 6)])
        ]
        total_pop = [1000, 500, 200]
        gdf = gpd.GeoDataFrame({'total_pop': total_pop, 'geometry': geometry}, crs="EPSG:4326")

        test_poly = Polygon([(1, 1), (2.5, 1), (2.5, 2.5), (1, 2.5)])

        result = self.toolkit.analysis.calculatePopulationInPolygon(
            shapelyPolygon=test_poly,
            dataSourceOrData=gdf,
            populationTypes="total_pop"
        )

        self.assertFalse(result.empty)
        total_estimated = result["total_pop"].sum()
        print(f"✅ Total estimated population in polygon: {total_estimated}")
        self.assertGreater(total_estimated, 0)
        self.assertLess(total_estimated, 1500)


class TestDemographyNotebook(unittest.TestCase):
    def test_notebook_runs_without_errors(self):
        """
        🚀 This test executes the DemographyToolkit Jupyter notebook and fails if any cell raises an error.
        """
        notebook_path = os.path.abspath(
            "doc/jupyter/toolkits/measurments/GIS/Demography/Demography_Toolkit_And_Unit_Test.ipynb"
        )
        self.assertTrue(os.path.exists(notebook_path), f"Notebook not found at {notebook_path}")

        with open(notebook_path, encoding='utf-8') as f:
            nb = nbformat.read(f, as_version=4)

        ep = ExecutePreprocessor(timeout=600, kernel_name='python3')

        try:
            ep.preprocess(nb, {'metadata': {'path': os.path.dirname(notebook_path)}})
        except Exception as e:
            self.fail(f"Notebook execution failed: {e}")

        print("✅ Notebook executed successfully with no errors.")

if __name__ == '__main__':
    unittest.main()
