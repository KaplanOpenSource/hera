import unittest
import os
import dask.dataframe as dd
import pandas as pd
import measurements.meteorology.highfreqdata.analysis.analysislayer as analysis



# ייבוא מחלקות ופונקציות לבדיקה מתוך ה- analysis
from .analysis.abstractcalculator import AbstractCalculator
from .analysis.meandatacalculator import MeanDataCalculator, AveragingCalculator
from .analysis.analysislayer import RawdataAnalysis
from .analysis.turbulencestatistics import singlePointTurbulenceStatistics



# שמור על ייבוא יחסי כמו בשאר הפרויקט
from .toolkit import HighFreqToolKit

def find_file_recursive(base_path, filename):
    """
    Recursively search for a file with a given name under base_path.
    Returns the full path if found, otherwise raises FileNotFoundError.
    """
    for root, dirs, files in os.walk(base_path):
        if filename in files:
            return os.path.join(root, filename)
    raise FileNotFoundError(f"File {filename} not found under {base_path}")

class TestHighFreqToolKit(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print("\n[SETUP] Setting up HighFreqToolKit for tests")

        # Require HERA_DATA_PATH to be set
        cls.HERA_DATA_PATH = os.getenv("HERA_DATA_PATH")
        if not cls.HERA_DATA_PATH:
            raise EnvironmentError(
                "Environment variable 'HERA_DATA_PATH' must be set. "
                "Please export it before running tests."
            )

        # 🧠 חיפוש חכם של קבצים בתוך HERA_DATA_PATH
        cls.sonic_file = find_file_recursive(cls.HERA_DATA_PATH, "slicedYamim_sonic.parquet")
        cls.trh_file = find_file_recursive(cls.HERA_DATA_PATH, "slicedYamim_TRH.parquet")

        # Debug print
        print(f"🔵 Found Sonic file at: {cls.sonic_file}")
        print(f"🔵 Found TRH file at: {cls.trh_file}")

        # Initialize the toolkit
        cls.toolkit = HighFreqToolKit(projectName="unittest_project")

    def test_docType_property(self):
        """✅ Test that the docType property returns the correct string"""
        print("\n[TEST] Running test_docType_property")
        expected = "highFreqMeteorology_HighFreqData"
        actual = self.toolkit.docType
        self.assertEqual(expected, actual)

    def test_read_sonic_data(self):
        """✅ Test reading the sonic parquet file into a Dask DataFrame"""
        print("\n[TEST] Running test_read_sonic_data")
        df = dd.read_parquet(self.sonic_file)
        self.assertIsInstance(df, dd.DataFrame)
        self.assertFalse(df.columns.empty)
        print(f"Loaded sonic data columns: {list(df.columns)}")

    def test_read_trh_data(self):
        """✅ Test reading the TRH parquet file into a Dask DataFrame"""
        print("\n[TEST] Running test_read_trh_data")
        df = dd.read_parquet(self.trh_file)
        self.assertIsInstance(df, dd.DataFrame)
        self.assertFalse(df.columns.empty)
        print(f"Loaded TRH data columns: {list(df.columns)}")

    def test_read_nonexistent_file(self):
        """✅ Test that trying to read a nonexistent file raises an appropriate error"""
        print("\n[TEST] Running test_read_nonexistent_file")
        fake_path = os.path.join(self.HERA_DATA_PATH, "nonexistent.parquet")
        with self.assertRaises((FileNotFoundError, AttributeError, OSError)):
            dd.read_parquet(fake_path)

    def test_sonic_time_range(self):
        """✅ Test that sonic data has a valid time range with timestamp or Time column"""
        print("\n[TEST] Running test_sonic_time_range")
        df = dd.read_parquet(self.sonic_file)
        df_pd = df.compute()
        time_column = None
        if "timestamp" in df_pd.columns:
            time_column = "timestamp"
        elif "Time" in df_pd.columns:
            time_column = "Time"
        self.assertIsNotNone(time_column, "No time column ('timestamp' or 'Time') found in sonic data.")
        start_time = df_pd[time_column].min()
        end_time = df_pd[time_column].max()
        self.assertLess(start_time, end_time)
        print(f"Sonic data time range: {start_time} -> {end_time}")

    def test_trh_time_range(self):
        """✅ Test that TRH data has a valid time range with timestamp or Time column"""
        print("\n[TEST] Running test_trh_time_range")
        df = dd.read_parquet(self.trh_file)
        df_pd = df.compute()
        time_column = None
        if "timestamp" in df_pd.columns:
            time_column = "timestamp"
        elif "Time" in df_pd.columns:
            time_column = "Time"
        self.assertIsNotNone(time_column, "No time column ('timestamp' or 'Time') found in TRH data.")
        start_time = df_pd[time_column].min()
        end_time = df_pd[time_column].max()
        self.assertLess(start_time, end_time)
        print(f"TRH data time range: {start_time} -> {end_time}")

    def test_campbelToParquet_with_nonexistent_file(self):
        """✅ Test that campbelToParquet handles nonexistent binary files gracefully"""
        print("\n[TEST] Running test_campbelToParquet_with_nonexistent_file")
        with self.assertRaises(ValueError):
            self.toolkit.campbelToParquet(binaryFile="/path/to/nonexistent_file.dat")

    def test_asciiToParquet_with_nonexistent_file(self):
        """✅ Test that asciiToParquet handles nonexistent ascii files gracefully"""
        print("\n[TEST] Running test_asciiToParquet_with_nonexistent_file")
        with self.assertRaises(FileNotFoundError):
            self.toolkit.asciiToParquet(path="/path/to/nonexistent_file.txt")

    def test_sonic_specific_point(self):
        """✅ Test specific known values from the sonic dataset."""
        print("\n[TEST] Running test_sonic_specific_point")
        df = dd.read_parquet(self.sonic_file).compute()

        first_row = df.iloc[0]

        self.assertAlmostEqual(first_row["u"], -2.65, places=2)
        self.assertAlmostEqual(first_row["v"], 3.05, places=2)
        self.assertAlmostEqual(first_row["w"], 1.96, places=2)  # 🔵 עדכון כאן
        self.assertAlmostEqual(first_row["T"], 24.93, places=2)

    def test_trh_specific_point(self):
        """✅ Test specific known values from the TRH dataset."""
        print("\n[TEST] Running test_trh_specific_point")
        df = dd.read_parquet(self.trh_file).compute()

        first_row = df.iloc[0]

        self.assertAlmostEqual(first_row["TC_T"], 24.71, places=2)
        self.assertAlmostEqual(first_row["RH"], 70.2, places=2)

        # --- טסטים חדשים ---

    def test_AbstractCalculator_init(self):
        """✅ Test initialization of AbstractCalculator"""
        print("\n[TEST] Running test_AbstractCalculator_init")
        import pandas as pd
        from measurements.meteorology.highfreqdata.analysis.abstractcalculator import AbstractCalculator

        df = pd.DataFrame({"u": [1, 2], "v": [3, 4], "w": [5, 6], "T": [7, 8]})
        ac = AbstractCalculator(rawData=df, metadata={"samplingWindow": 10})
        self.assertIsNotNone(ac)
        self.assertIsInstance(ac.RawData, pd.DataFrame)

    def test_MeanDataCalculator_calculate_mean(self):
        """✅ Test basic mean calculation with MeanDataCalculator"""
        print("\n[TEST] Running test_MeanDataCalculator_calculate_mean")
        df = dd.read_parquet(self.sonic_file).compute()

        # תיקון אינדקס
        if not isinstance(df.index, pd.DatetimeIndex):
            if "timestamp" in df.columns:
                df["timestamp"] = pd.to_datetime(df["timestamp"])
                df = df.set_index("timestamp")
            elif "Time" in df.columns:
                df["Time"] = pd.to_datetime(df["Time"])
                df = df.set_index("Time")
            else:
                raise ValueError("No timestamp or Time column found for setting datetime index.")

        metadata = {
            "samplingWindow": "10S",
            "isMissingData": False,
            "start": df.index.min(),
            "end": df.index.max(),
        }

        turb_stats = singlePointTurbulenceStatistics(df, metadata=metadata)
        mc = MeanDataCalculator(rawData=df, metadata=metadata, TurbCalcOrData=turb_stats)

        self.assertIsNotNone(mc.MeanData)
        self.assertFalse(mc.MeanData.empty)
        print("✅ MeanData calculated successfully.")

    def test_RawdataAnalysis_process(self):
        """✅ Test basic functionality of RawdataAnalysis"""
        print("\n[TEST] Running test_RawdataAnalysis_process")
        df = dd.read_parquet(self.sonic_file).compute()

        metadata = {
            "samplingWindow": "10S",
            "isMissingData": False,
            "start": df.index.min(),
            "end": df.index.max(),
        }

        # נבנה נכון את analysis
        analysis = RawdataAnalysis(datalayer=self.toolkit)

        # נשתמש בפונקציה מתוך analysis, לא ישירות ב-init
        turb_calc = analysis.singlePointTurbulenceStatistics(
            sonicData=df,
            samplingWindow=metadata["samplingWindow"],
            start=metadata["start"],
            end=metadata["end"],
            height=10,
            buildingHeight=0,
            averagedHeight=10,
            inmemory=True,
            isMissingData=metadata["isMissingData"]
        )

        self.assertIsNotNone(turb_calc)
        print("✅ RawdataAnalysis singlePointTurbulenceStatistics worked successfully.")

    def test_singlePointTurbulenceStatistics(self):
        """✅ Test basic initialization of singlePointTurbulenceStatistics"""
        print("\n[TEST] Running test_singlePointTurbulenceStatistics")
        import pandas as pd
        from measurements.meteorology.highfreqdata.analysis.turbulencestatistics import singlePointTurbulenceStatistics

        df = pd.DataFrame({"u": [0.5, 1.0], "v": [0.1, 0.2], "w": [0.3, 0.4], "T": [10, 11]})
        stats = singlePointTurbulenceStatistics(df, metadata={"samplingWindow": 10})
        self.assertIsNotNone(stats)

    # =====================================
    # 🔵 AbstractCalculator unit tests
    # =====================================

    def test_AbstractCalculator_init_basic(self):
        """✅ Test that AbstractCalculator initializes correctly with rawData and metadata."""
        print("\n[TEST] Running test_AbstractCalculator_init_basic")
        df = dd.read_parquet(self.sonic_file).compute()
        ac = AbstractCalculator(rawData=df, metadata={"samplingWindow": "10S", "isMissingData": False})
        self.assertIsNotNone(ac.RawData)
        self.assertEqual(ac.metaData["samplingWindow"], "10S")

    def test_AbstractCalculator_sampling_window(self):
        """✅ Test that SamplingWindow property returns the correct sampling window."""
        print("\n[TEST] Running test_AbstractCalculator_sampling_window")
        df = dd.read_parquet(self.sonic_file).compute()
        ac = AbstractCalculator(rawData=df, metadata={"samplingWindow": "10S", "isMissingData": False})
        self.assertEqual(ac.SamplingWindow, "10S")

    def test_AbstractCalculator_compute_methods_exist(self):
        """✅ Test that compute and _compute methods exist in AbstractCalculator."""
        print("\n[TEST] Running test_AbstractCalculator_compute_methods_exist")
        df = dd.read_parquet(self.sonic_file).compute()
        ac = AbstractCalculator(rawData=df, metadata={"samplingWindow": "10S", "isMissingData": False})
        self.assertTrue(hasattr(ac, "compute"))
        self.assertTrue(hasattr(ac, "_compute"))

    def test_AbstractCalculator_set_save_properties(self):
        """✅ Test that set_saveProperties either returns a dict or None (if data missing)."""
        print("\n[TEST] Running test_AbstractCalculator_set_save_properties")
        df = dd.read_parquet(self.sonic_file).compute()
        ac = AbstractCalculator(rawData=df, metadata={"samplingWindow": "10S"})

        props = ac.set_saveProperties(dataFormat="testFormat")

        # Allow both dict or None
        self.assertTrue(props is None or isinstance(props, dict), "set_saveProperties should return dict or None")

    # 📂 Starting Tests for MeanDataCalculator (Real Functions, Fixed)

    import pandas as pd
    import dask.dataframe as dd
    from measurements.meteorology.highfreqdata.analysis.meandatacalculator import MeanDataCalculator
    from measurements.meteorology.highfreqdata.analysis.turbulencestatistics import singlePointTurbulenceStatistics

    # 🧪 Test 1: hour and timeWithinDay
    def test_MeanDataCalculator_hour_and_timeWithinDay(self):
        """✅ Test that 'hour' and 'timeWithinDay' columns are created."""
        print("\n[TEST] Running test_MeanDataCalculator_hour_and_timeWithinDay")
        df = dd.read_parquet(self.sonic_file).compute()
        df = df.set_index("Time")
        df["TC_T_bar"] = df["T"]  # ✅ Adding required column temporarily
        metadata = {"samplingWindow": "10S", "isMissingData": False, "start": df.index.min(), "end": df.index.max()}
        turb_stats = singlePointTurbulenceStatistics(df, metadata=metadata)
        calc = MeanDataCalculator(TurbCalcOrData=turb_stats)
        calc.hour()
        calc.timeWithinDay()
        self.assertIn("hour", calc.MeanData.columns)
        self.assertIn("timeWithinDay", calc.MeanData.columns)

    # 🧪 Test 2: horizontalSpeed
    def test_MeanDataCalculator_horizontalSpeed(self):
        """✅ Test creation of 'horizontal_speed_bar'."""
        print("\n[TEST] Running test_MeanDataCalculator_horizontalSpeed")
        df = dd.read_parquet(self.sonic_file).compute()
        df = df.set_index("Time")
        df["TC_T_bar"] = df["T"]
        metadata = {"samplingWindow": "10S", "isMissingData": False, "start": df.index.min(), "end": df.index.max()}
        turb_stats = singlePointTurbulenceStatistics(df, metadata=metadata)
        calc = MeanDataCalculator(TurbCalcOrData=turb_stats)
        calc.horizontalSpeed()
        self.assertIn("horizontal_speed_bar", calc.MeanData.columns)

    # 🧪 Test 3: sigma and sigmaH
    def test_MeanDataCalculator_sigma_sigmaH(self):
        """✅ Test creation of sigma columns."""
        print("\n[TEST] Running test_MeanDataCalculator_sigma_sigmaH")
        df = dd.read_parquet(self.sonic_file).compute()
        df = df.set_index("Time")
        df["TC_T_bar"] = df["T"]
        metadata = {"samplingWindow": "10S", "isMissingData": False, "start": df.index.min(), "end": df.index.max()}
        turb_stats = singlePointTurbulenceStatistics(df, metadata=metadata)
        calc = MeanDataCalculator(TurbCalcOrData=turb_stats)
        calc.sigma()
        calc.sigmaH()
        for col in ["sigmaU", "sigmaV", "sigmaW", "sigmaH"]:
            self.assertIn(col, calc.MeanData.columns)

    # 🧪 Test 4: Ustar and uStarOverWindSpeed
    def test_MeanDataCalculator_Ustar_and_uStarOverWindSpeed(self):
        """✅ Test creation of Ustar and uStarOverWindSpeed."""
        print("\n[TEST] Running test_MeanDataCalculator_Ustar_and_uStarOverWindSpeed")
        df = dd.read_parquet(self.sonic_file).compute()
        df = df.set_index("Time")
        df["TC_T_bar"] = df["T"]
        metadata = {"samplingWindow": "10S", "isMissingData": False, "start": df.index.min(), "end": df.index.max()}
        turb_stats = singlePointTurbulenceStatistics(df, metadata=metadata)
        calc = MeanDataCalculator(TurbCalcOrData=turb_stats)
        calc.Ustar()
        calc.uStarOverWindSpeed()
        self.assertIn("Ustar", calc.MeanData.columns)
        self.assertIn("uStarOverWindSpeed", calc.MeanData.columns)

    def test_TKE_adds_column(self):
        print("\n[TEST] Running test_TKE_adds_column")
        from .analysis.meandatacalculator import MeanDataCalculator
        from .analysis.turbulencestatistics import singlePointTurbulenceStatistics

        df = dd.read_parquet(self.sonic_file).compute().head(100)
        df["Time"] = pd.to_datetime(df["Time"], errors="coerce")
        df = df.dropna(subset=["Time"]).set_index("Time").sort_index()

        metadata = {
            "samplingWindow": "10S",
            "isMissingData": False,
            "start": df.index.min(),
            "end": df.index.max()

        }

        turb = singlePointTurbulenceStatistics(df, metadata)
        turb.secondMoments()

        calc = MeanDataCalculator(TurbCalcOrData=turb, metadata=metadata)
        calc.TKE()

        self.assertIn("TKE", calc.MeanData.columns)

    def test_MOLength_adds_column(self):
        print("\n[TEST] Running test_MOLength_adds_column")
        from .analysis.meandatacalculator import MeanDataCalculator
        from .analysis.turbulencestatistics import singlePointTurbulenceStatistics

        df = dd.read_parquet(self.sonic_file).compute().head(100)
        df["Time"] = pd.to_datetime(df["Time"], errors="coerce")
        df = df.dropna(subset=["Time"]).set_index("Time").sort_index()

        metadata = {
            "samplingWindow": "10S",
            "isMissingData": False,
            "start": df.index.min(),
            "end": df.index.max()
        }

        turb = singlePointTurbulenceStatistics(df, metadata)
        turb.secondMoments()

        calc = MeanDataCalculator(TurbCalcOrData=turb, metadata=metadata)

        calc.MeanData["TC_T_bar"] = 20
        calc.MeanData["Ustar"] = 0.3

        calc.MOLength()
        self.assertIn("L", calc.MeanData.columns)

    # 🧪 Test 7: compute method returns DataFrame
    def test_MeanDataCalculator_compute_returns_dataframe(self):
        """✅ Test compute method returns a pandas DataFrame."""
        print("\n[TEST] Running test_MeanDataCalculator_compute_returns_dataframe")
        df = dd.read_parquet(self.sonic_file).compute()
        df = df.set_index("Time")
        df["TC_T_bar"] = df["T"]
        metadata = {"samplingWindow": "10S", "isMissingData": False, "start": df.index.min(), "end": df.index.max()}
        turb_stats = singlePointTurbulenceStatistics(df, metadata=metadata)
        calc = MeanDataCalculator(TurbCalcOrData=turb_stats)
        result = calc.compute()
        self.assertIsInstance(result, pd.DataFrame)

    # ✅ End of real tests for MeanDataCalculator

    # ============================================
    # 🔎 Unit Tests for RawdataAnalysis (analysislayer.py)
    # ============================================

    def test_singlePointTurbulenceStatistics_returns_instance(self):
        """
        ✅ Test that singlePointTurbulenceStatistics returns a valid result
        when provided with correct sonic data and metadata.
        """
        print("\n[TEST] Running test_singlePointTurbulenceStatistics_returns_instance")
        df = dd.read_parquet(self.sonic_file).compute()
        df["Time"] = pd.to_datetime(df["Time"])
        df = df.set_index("Time")

        metadata = {
            "samplingWindow": "10S",
            "start": df.index.min(),
            "end": df.index.max(),
            "height": 10,
            "buildingHeight": 5,
            "averagedHeight": 7,
            "isMissingData": False,
        }

        analysis = RawdataAnalysis(datalayer=self.toolkit)
        result = analysis.singlePointTurbulenceStatistics(
            sonicData=df,
            samplingWindow=metadata["samplingWindow"],
            start=metadata["start"],
            end=metadata["end"],
            height=metadata["height"],
            buildingHeight=metadata["buildingHeight"],
            averagedHeight=metadata["averagedHeight"],
            inmemory=True,
            isMissingData=metadata["isMissingData"]
        )
        self.assertIsNotNone(result)

    def test_singlePointTurbulenceStatistics_raises_on_invalid_input(self):
        """
        ❌ Test that singlePointTurbulenceStatistics raises an error on invalid input.
        """
        print("\n[TEST] Running test_singlePointTurbulenceStatistics_raises_on_invalid_input")
        analysis = RawdataAnalysis(datalayer=self.toolkit)
        with self.assertRaises(TypeError):
            analysis.singlePointTurbulenceStatistics(sonicData="not_a_dataframe")

    # === ANALYSISLAYER: Tests for AveragingCalculator ===

    def test_AveragingCalculator_returns_instance(self):
        """✅ Test that AveragingCalculator returns a valid DataFrame when passed correct parameters."""
        print("\n[TEST] Running test_AveragingCalculator_returns_instance")

        df = pd.read_parquet(self.sonic_file)

        # Try to parse Time column
        if "Time" in df.columns:
            try:
                df["Time"] = pd.to_datetime(df["Time"])
            except Exception:
                df.drop(columns=["Time"], inplace=True)

        if "Time" not in df.columns:
            df["Time"] = pd.date_range("2020-01-01", periods=len(df), freq="1S")

        df = df.set_index("Time").sort_index()

        # יצירת מופע של RawdataAnalysis עם datalayer מה-toolkit
        analysis_layer = RawdataAnalysis(datalayer=self.toolkit)

        # קריאה לפי סדר בלבד (positional), כי זו ההגדרה של הפונקציה
        calc = analysis_layer.AveragingCalculator(
            df,
            "5min",
            df.index[0],
            df.index[-1],
            2,
            10,
            8
        )

        result = calc.compute()  # מאחר וזה מחזיר Calculator – נדרשת קריאה ל־compute

        self.assertIsInstance(result, pd.DataFrame)
        self.assertFalse(result.empty)

    def test_AveragingCalculator_raises_on_invalid_input(self):
        """
        ❌ Test that AveragingCalculator raises an error on invalid input.
        """
        print("\n[TEST] Running test_AveragingCalculator_raises_on_invalid_input")
        analysis = RawdataAnalysis(datalayer=self.toolkit)
        with self.assertRaises(TypeError):
            analysis.AveragingCalculator(sonicData=None)


    # =============================
    # 🔬 Unit tests for turbulencestatistics.py
    # =============================

import unittest
import pandas as pd
from measurements.meteorology.highfreqdata.analysis.turbulencestatistics import singlePointTurbulenceStatistics

class TestSinglePointTurbulenceStatisticsMethods(unittest.TestCase):
    def setUp(self):
        # Dummy DataFrame with values
        self.df = pd.DataFrame({
            "u": [1.0, 2.0, 3.0],
            "v": [0.5, 1.5, 2.5],
            "w": [0.1, 0.2, 0.3],
            "T": [20.0, 21.0, 22.0]
        })

        # Metadata for testing
        self.metadata = {
            "samplingWindow": "10S",
            "start": pd.Timestamp("2020-01-01 00:00:00"),
            "end": pd.Timestamp("2020-01-01 00:00:02"),
            "height": 10,
            "buildingHeight": 2,
            "averagedHeight": 6,
            "isMissingData": False
        }

        # Add a datetime index to support .resample()
        self.df["Time"] = pd.date_range(start=self.metadata["start"], periods=len(self.df), freq="1S")
        self.df.set_index("Time", inplace=True)

    def test_instantiation(self):
        calc = singlePointTurbulenceStatistics(self.df, self.metadata)
        self.assertIsNotNone(calc)
        self.assertTrue(hasattr(calc, "RawData"))
        self.assertTrue(hasattr(calc, "metaData"))

    def test_invalid_input_type(self):
        """❌ Expect ValueError when input is not a DataFrame"""
        with self.assertRaises(ValueError):
            singlePointTurbulenceStatistics("not_a_dataframe", self.metadata)

    def test_fluctuations_outputs(self):
        calc = singlePointTurbulenceStatistics(self.df.copy(), self.metadata)
        calc.fluctuations()
        for col in ["u_bar", "v_bar", "w_bar", "T_bar", "up", "vp", "wp", "Tp"]:
            self.assertIn(col, calc.RawData.columns)

    def test_secondMoments_outputs_exist(self):
        calc = singlePointTurbulenceStatistics(self.df.copy(), self.metadata)
        calc.secondMoments()
        for col in ["uu", "vv", "ww", "uw", "uv", "uT", "vT", "wT", "TT", "vw"]:
            self.assertIn(col, calc.TemporaryData.columns)

    def test_sigma_computation(self):
        calc = singlePointTurbulenceStatistics(self.df.copy(), self.metadata)
        calc.sigma()
        for col in ["sigmaU", "sigmaV", "sigmaW"]:
            self.assertIn(col, calc.TemporaryData.columns)

    def test_horizontal_speed(self):
        calc = singlePointTurbulenceStatistics(self.df.copy(), self.metadata)
        calc.horizontalSpeed()
        self.assertIn("horizontal_speed_bar", calc.TemporaryData.columns)

    def test_Ustar_computation(self):
        calc = singlePointTurbulenceStatistics(self.df.copy(), self.metadata)
        calc.Ustar()
        self.assertIn("Ustar", calc.TemporaryData.columns)

    def test_TKE_computation(self):
        calc = singlePointTurbulenceStatistics(self.df.copy(), self.metadata)
        calc.uu().vv().ww().TKE()
        self.assertIn("TKE", calc.TemporaryData.columns)

    def test_MOLength_Sonic_computation(self):
        calc = singlePointTurbulenceStatistics(self.df.copy(), self.metadata)
        calc.wT().Ustar().MOLength_Sonic()
        self.assertIn("L_Sonic", calc.TemporaryData.columns)


if __name__ == "__main__":
    unittest.main()
