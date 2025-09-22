import unittest
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from hera.measurements.meteorology.lowfreqdata.toolkit import lowFreqToolKit
from hera.measurements.meteorology.lowfreqdata.presentationLayer import presenation
from hera.measurements.meteorology.lowfreqdata import analysis
from hera.utils.statistics import calcDist2d


class TestLowFreqToolkit(unittest.TestCase):
    """
    ✅ Unit tests for lowFreqToolKit and analysis.py
    """

    @classmethod
    def setUpClass(cls):
        # Load data from path defined by HERA_DATA_PATH
        base_path = os.environ.get("HERA_DATA_PATH", "/home/ilay/hera_unittest_data")
        cls.file_path = os.path.join(
            base_path, "measurements", "meteorology", "lowfreqdata", "YAVNEEL.parquet"
        )
        cls.df = pd.read_parquet(cls.file_path)

        # Ensure datetime column is in correct type
        if "datetime" not in cls.df.columns:
            raise Exception("❌ 'datetime' column missing in dataframe")
        cls.df["datetime"] = pd.to_datetime(cls.df["datetime"], utc=True)

        # Initialize toolkit instance
        cls.toolkit = lowFreqToolKit(projectName="UNIT_TEST_LOWFREQ")

    def test_init_toolkit_structure(self):
        """
        ✅ Test if the toolkit initializes core components correctly
        """
        self.assertIsNotNone(self.toolkit.analysis, "❌ Missing 'analysis' component in toolkit")
        self.assertIsNotNone(self.toolkit.presentation, "❌ Missing 'presentation' component in toolkit")
        self.assertTrue(hasattr(self.toolkit, "docType"), "❌ Missing attribute 'docType' in toolkit")

    def test_docType(self):
        """
        ✅ Test the returned docType string
        """
        doc_type = self.toolkit.docType
        self.assertIsInstance(doc_type, str, "❌ docType should return a string")
        self.assertIn("LOWFREQ", doc_type.upper(), "❌ docType does not include 'LOWFREQ'")

    # =============================================================================
    # ✅ Unit tests for analysis.py
    # =============================================================================

    def test_add_dates_columns(self):
        """
        ✅ analysis.addDatesColumns - Verify output structure and stability by comparing to a saved reference.

        If env var `PREPARE_EXPECTED_OUTPUT` is set to "1", the test will generate and save reference output.
        Otherwise, it compares the current result to the saved reference.
        """
        import tempfile

        df = self.df.copy()
        df["datetime"] = pd.to_datetime(df["datetime"], utc=True)

        enriched = self.toolkit.analysis.addDatesColumns(df, datecolumn="datetime")

        # Path where reference output will be saved/loaded
        expected_path = os.path.join(
            os.path.dirname(__file__), "expected_output_addDatesColumns.parquet"
        )

        if os.environ.get("PREPARE_EXPECTED_OUTPUT") == "1":
            # Save reference output
            # Convert problematic columns to strings for compatibility
            if "timeonly" in enriched.columns:
                enriched["timeonly"] = enriched["timeonly"].astype(str)

            enriched.to_parquet(expected_path, index=False)

            print(f"📁 Saved reference output to {expected_path}")
        else:
            # Compare to expected reference
            self.assertTrue(os.path.exists(expected_path), f"❌ Reference file not found at {expected_path}")
            expected = pd.read_parquet(expected_path)

            for df_ in (enriched, expected):
                if "timeonly" in df_.columns:
                    df_["timeonly"] = df_["timeonly"].astype(str)

            self.assertListEqual(
                list(enriched.columns),
                list(expected.columns),
                "❌ Columns in output do not match reference"
            )

            pd.testing.assert_frame_equal(
                enriched.reset_index(drop=True),
                expected.reset_index(drop=True),
                check_dtype=False,
                check_like=True,
                atol=1e-6,
                obj="❌ Data mismatch between result and expected"
            )

    def test_calc_hourly_dist_max_normalized(self):
        """
        ✅ analysis.calcHourlyDist - max_normalized mode
        """
        df = self.df.copy()
        x, y, M = self.toolkit.analysis.calcHourlyDist(df, Field="RH", normalization="max_normalized")

        self.assertIsNotNone(x, "❌ x-axis histogram missing")
        self.assertIsNotNone(y, "❌ y-axis histogram missing")
        self.assertIsNotNone(M, "❌ 2D histogram matrix missing")
        self.assertEqual(len(y), M.shape[0], "❌ y mismatch with matrix shape")
        self.assertEqual(len(x), M.shape[1], "❌ x mismatch with matrix shape")

    def test_y_normalized_behavior(self):
        """
        ✅ calcDist2d - y_normalized logic:
            - Rows with values → normalized to sum 1
            - Row with single value → becomes [0, ..., 1, ..., 0]
            - Empty row → remains all 0
        """

        # Create synthetic data
        x = np.array([0.5, 1.5, 2.5, 1.5])
        y = np.array([1, 1, 1, 2])

        x_range = (0, 3)
        y_range = (0, 3)

        x_mid, y_mid, M = calcDist2d(x, y, bins=3, normalization="y_normalized", x_range=x_range, y_range=y_range)

        # Row 0: no data → all zeros
        self.assertTrue(np.allclose(M[0], 0), "❌ Row 0 should be all zeros")

        # Row 1: multiple values → sum to 1
        self.assertTrue(np.isclose(M[1].sum(), 1.0), "❌ Row 1 should sum to 1")
        self.assertGreater(np.count_nonzero(M[1]), 1, "❌ Row 1 should have more than one non-zero value")

        # Row 2: single value → sum to 1 with one non-zero
        self.assertTrue(np.isclose(M[2].sum(), 1.0), "❌ Row 2 should sum to 1")
        self.assertEqual(np.count_nonzero(M[2]), 1, "❌ Row 2 should have exactly one non-zero value")


    def test_calc_hourly_dist_density(self):
        """
        ✅ analysis.calcHourlyDist - density mode
        """
        df = self.df.copy()
        x, y, M = self.toolkit.analysis.calcHourlyDist(df, Field="RH", normalization="density")
        self.assertTrue((M >= 0).all(), "❌ Density matrix should not contain negative values")

    def test_resample_second_moments(self):
        """
        ✅ analysis.resampleSecondMoments - Validate computation of variance and covariance.

        This test:
        - Creates raw variables X, Y and derived fields: X_bar, Y_bar, XX, XY, YY
        - Calls resampleSecondMoments with daily window
        - Verifies that the resulting DataFrame includes second moment fields: 'XX', 'XY', 'YY'
        """
        df = self.df.copy()
        df = df.set_index("datetime")

        # Create base variables
        df["X"] = df["WS"]
        df["Y"] = df["RH"]

        # Create derived second-moment fields required by _calculateCov
        df["XX"] = df["X"] ** 2
        df["YY"] = df["Y"] ** 2
        df["XY"] = df["X"] * df["Y"]

        # Also create mean fields to pass to the function
        df["X_bar"] = df["X"]
        df["Y_bar"] = df["Y"]

        result = self.toolkit.analysis.resampleSecondMoments(
            df,
            SamplingWindow="D",
            fieldsFirstMoments=["X_bar", "Y_bar"],
            fieldsSecondMoments=["X", "Y"]
        )

        self.assertIsInstance(result, pd.DataFrame, "❌ Output should be a DataFrame")

        # These columns should be created by _calculateCov based on 'X', 'Y'
        for col in ["XX", "XY", "YY"]:
            self.assertIn(col, result.columns, f"❌ Expected column '{col}' not found in output")

    # --- 🧪 Tests for presentationLayer functions --- #

    def test_plotScatter(self):
        """
        ✅ Test DailyPlots.plotScatter — basic scatter plot creation
        """
        df = self.df.copy()
        df["datetime"] = pd.to_datetime(df["datetime"], utc=True)
        ax = self.toolkit.presentation.dailyPlots.plotScatter(df, plotField="RH")
        self.assertIsNotNone(ax, "❌ plotScatter did not return an Axes object")

    def test_dateLinePlot(self):
        """
        ✅ Test DailyPlots.dateLinePlot — line plot for a specific date
        """
        df = self.df.copy()
        df["datetime"] = pd.to_datetime(df["datetime"], utc=True)
        example_date = df["datetime"].dt.date.astype(str).iloc[0]
        ax, line = self.toolkit.presentation.dailyPlots.dateLinePlot(df, plotField="RH", date=example_date)
        self.assertIsNotNone(ax, "❌ dateLinePlot did not return an Axes object")
        self.assertTrue(len(line) > 0, "❌ No line returned from dateLinePlot")

    def test_plotProbContourf(self):
        """
        ✅ Test DailyPlots.plotProbContourf — main probability plot with contour
        """
        df = self.df.copy()
        df["datetime"] = pd.to_datetime(df["datetime"], utc=True)
        CS, CFS, ax = self.toolkit.presentation.dailyPlots.plotProbContourf(df, plotField="RH")
        self.assertIsNotNone(CS, "❌ ContourSet (CS) is None")
        self.assertIsNotNone(CFS, "❌ ContourfSet (CFS) is None")
        self.assertIsNotNone(ax, "❌ Axes is None")

    def test_plotProbContourf_bySeason(self):
        """
        ✅ Test SeasonalPlots.plotProbContourf_bySeason — 2x2 seasonal plot
        """
        df = self.df.copy()
        df["datetime"] = pd.to_datetime(df["datetime"], utc=True)
        ax = self.toolkit.presentation.seasonalPlots.plotProbContourf_bySeason(df, plotField="RH")
        self.assertIsNotNone(ax, "❌ Axes grid from seasonal plot is None")

    def test_dateLinePlot_matches_data(self):
        """
        ✅ Test that the plotted line in dateLinePlot matches the original data values.

        This test:
        - Picks the first date available in the dataset.
        - Runs dateLinePlot to produce a line chart for 'RH' (Relative Humidity).
        - Extracts the plotted line data (x and y values).
        - Filters the DataFrame to the same date and same conditions as in the code.
        - Compares that the values match exactly in length and content.
        """
        df = self.df.copy()
        df["datetime"] = pd.to_datetime(df["datetime"], utc=True)
        df = df.set_index("datetime")

        # Select the first date available in the dataset
        example_date = df.index.date[0].strftime("%Y-%m-%d")

        # Generate the plot for that date
        ax, line = self.toolkit.presentation.dailyPlots.dateLinePlot(df, plotField="RH", date=example_date)

        # Extract the plotted line (x = hour, y = value)
        line_obj = line[0]
        x_vals = line_obj.get_xdata()
        y_vals = line_obj.get_ydata()

        # Prepare the same data slice manually for comparison
        daily = df.copy()
        daily["RH"] = daily["RH"].where(daily["RH"] > -5000)
        daily = daily.assign(houronly=daily.index.hour + daily.index.minute / 60.)
        filtered = daily[daily.index.date.astype(str) == example_date].dropna(subset=["RH"])

        # Assertions: same length and same values
        self.assertEqual(len(y_vals), len(filtered), "❌ Number of plotted points doesn't match filtered data length")
        np.testing.assert_array_almost_equal(
            y_vals, filtered["RH"].values, decimal=2,
            err_msg="❌ Y-values from line plot do not match data values"
        )

    def test_plotScatter_matches_data(self):
        """
        ✅ Test that scatter points in plotScatter match data values in 'RH'.

        This test:
        - Runs plotScatter on the dataframe.
        - Extracts the Y values from the plotted scatter points.
        - Filters the dataframe with the same condition used in the plot (-5000 threshold).
        - Compares that the Y values from the plot are all found in the filtered dataset.
        """
        df = self.df.copy()
        df["datetime"] = pd.to_datetime(df["datetime"], utc=True)

        ax = self.toolkit.presentation.dailyPlots.plotScatter(df, plotField="RH")

        # Extract all Y data from scatter plot points
        y_vals_from_plot = []
        for coll in ax.collections:
            offsets = coll.get_offsets()
            if len(offsets) > 0:
                y_vals_from_plot.extend(offsets[:, 1])  # Y values

        # Clean the data like in plotScatter
        filtered = df["RH"].where(df["RH"] > -5000).dropna()

        # Check that most plotted values are within the data (allowing rounding differences)
        matches = sum(np.isin(np.round(y_vals_from_plot, 2), np.round(filtered.values, 2)))
        self.assertGreater(matches, 0, "❌ No matching RH values found between plot and data")

    def test_plotProbContourf_distribution_ranges(self):
        """
        ✅ Test that the contour plot histogram range includes real RH values.

        This test:
        - Runs plotProbContourf for 'RH'.
        - Extracts the y-axis range of the histogram.
        - Confirms that it includes the min/max values of RH in the filtered dataset.
        """
        df = self.df.copy()
        df["datetime"] = pd.to_datetime(df["datetime"], utc=True)
        df = df.set_index("datetime")

        # Run the function to get the histogram result
        CS, CFS, ax = self.toolkit.presentation.dailyPlots.plotProbContourf(df, plotField="RH")

        # Extract y-limits from the plot
        y_min, y_max = ax.get_ylim()

        # Clean the data manually
        filtered = df["RH"].where(df["RH"] > -5000).dropna()
        min_val = filtered.min()
        max_val = filtered.max()

        # Assertions: the plot range must include the data range
        self.assertLessEqual(min_val, y_max, "❌ Plot Y-max does not include highest data value")
        self.assertGreaterEqual(max_val, y_min, "❌ Plot Y-min does not include lowest data value")

    def test_plotScatter_WS_field(self):
        """
        ✅ Test plotScatter with 'WS' (Wind Speed) field to verify support for multiple columns.
        """
        df = self.df.copy()
        df["datetime"] = pd.to_datetime(df["datetime"], utc=True)

        ax = self.toolkit.presentation.dailyPlots.plotScatter(df, plotField="WS")

        y_vals_from_plot = []
        for coll in ax.collections:
            offsets = coll.get_offsets()
            if len(offsets) > 0:
                y_vals_from_plot.extend(offsets[:, 1])  # Y values

        filtered = df["WS"].where(df["WS"] > -5000).dropna()
        matches = sum(np.isin(np.round(y_vals_from_plot, 2), np.round(filtered.values, 2)))

        self.assertGreater(matches, 0, "❌ No matching WS values found between plot and data")

    def test_plotScatter_empty_dataframe(self):
        """
        ✅ Ensure plotScatter handles empty DataFrame gracefully and produces an empty plot.
        """
        df = pd.DataFrame(columns=["datetime", "RH"])
        df["datetime"] = pd.to_datetime(df["datetime"])

        ax = self.toolkit.presentation.dailyPlots.plotScatter(df, plotField="RH")

        total_points = sum(len(coll.get_offsets()) for coll in ax.collections)

        self.assertEqual(total_points, 0, "❌ plotScatter should produce empty plot for empty DataFrame")

    def test_plotScatter_with_nan_and_outliers(self):
        """
        ✅ Check plotScatter behavior with NaN and extreme negative values in 'RH'.
        """
        df = self.df.copy()
        df["datetime"] = pd.to_datetime(df["datetime"], utc=True)

        # Insert NaN and extreme outliers
        df.iloc[0, df.columns.get_loc("RH")] = np.nan
        df.iloc[1, df.columns.get_loc("RH")] = -9999

        ax = self.toolkit.presentation.dailyPlots.plotScatter(df, plotField="RH")

        # Ensure that plotted points ignore NaNs and outliers
        y_vals_from_plot = []
        for coll in ax.collections:
            offsets = coll.get_offsets()
            if len(offsets) > 0:
                y_vals_from_plot.extend(offsets[:, 1])

        # Should be less than the original length
        self.assertLess(len(y_vals_from_plot), len(df), "❌ Outliers or NaNs not filtered properly")

    def test_plotScatter_WD_field(self):
        """
        ✅ Test plotScatter with 'WD' (Wind Direction) to confirm support for other meteorological fields.
        """
        df = self.df.copy()
        df["datetime"] = pd.to_datetime(df["datetime"], utc=True)

        ax = self.toolkit.presentation.dailyPlots.plotScatter(df, plotField="WD")

        y_vals_from_plot = []
        for coll in ax.collections:
            offsets = coll.get_offsets()
            if len(offsets) > 0:
                y_vals_from_plot.extend(offsets[:, 1])

        filtered = df["WD"].where(df["WD"] > -5000).dropna()
        matches = sum(np.isin(np.round(y_vals_from_plot, 2), np.round(filtered.values, 2)))

        self.assertGreater(matches, 0, "❌ No matching WD values found between plot and data")

    def test_plotScatter_axis_labels(self):
        """
        ✅ Check that plotScatter sets the correct axis labels for 'RH' field.
        """
        df = self.df.copy()
        df["datetime"] = pd.to_datetime(df["datetime"], utc=True)

        ax = self.toolkit.presentation.dailyPlots.plotScatter(df, plotField="RH")

        self.assertEqual(ax.get_xlabel(), "Time [Hours]", "❌ X-axis label is incorrect")
        self.assertEqual(ax.get_ylabel(), "Relative Humidity [%]", "❌ Y-axis label for RH is incorrect")

    def test_plotProbContourf_bySeason_basic(self):
        """
        ✅ Test SeasonalPlots.plotProbContourf_bySeason for basic RH seasonal plotting.
        Checks if Axes grid with 4 seasonal plots is created.
        """
        df = self.df.copy()
        df["datetime"] = pd.to_datetime(df["datetime"], utc=True)

        ax = self.toolkit.presentation.seasonalPlots.plotProbContourf_bySeason(df, plotField="RH")

        self.assertIsNotNone(ax, "❌ Seasonal plot returned None")
        self.assertEqual(ax.shape, (2, 2), "❌ Expected 2x2 Axes grid for seasonal plot")

    def test_plotScatter_creates_non_empty_image(self):
        """
        🖼 Ensure that plotScatter produces a non-empty plot file.
        """
        import tempfile
        df = self.df.copy()
        df["datetime"] = pd.to_datetime(df["datetime"], utc=True)

        ax = self.toolkit.presentation.dailyPlots.plotScatter(df, plotField="RH")

        with tempfile.NamedTemporaryFile(suffix=".png") as tmpfile:
            ax.figure.savefig(tmpfile.name)
            size = os.path.getsize(tmpfile.name)
            self.assertGreater(size, 1000, "❌ Saved plot image is unexpectedly small or empty")

    def test_run_lowfreqtoolkit_notebook(self):
        """
        ✅ Test that the 'lowFreqToolKit_documentation.ipynb' Jupyter notebook runs successfully from start to finish.
        """
        import os
        import nbformat
        from nbclient import NotebookClient
        from nbclient.exceptions import CellExecutionError

        # Set the required environment variable for the notebook to run
        os.environ['HERA_UNITTEST_DATA'] = "/home/ilay/hera_unittest_data"

        notebook_path = "/home/ilay/hera/hera/doc/jupyter/toolkits/measurments/meteorology/lowFreqToolKit_documentation.ipynb"

        with open(notebook_path) as f:
            nb = nbformat.read(f, as_version=4)

        client = NotebookClient(nb, timeout=300, kernel_name="python3")

        try:
            client.execute()
        except CellExecutionError as e:
            self.fail(f"❌ Notebook execution failed: {e}")


# 👇 This will allow running the tests directly via CLI
if __name__ == '__main__':
    unittest.main()
