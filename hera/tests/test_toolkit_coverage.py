"""
Baseline I/O regression tests for primary Hera toolkits.

Covers: TopographyToolkit, LSMTemplate, OFToolkit, BriggsRural (Gaussian),
        LowFreq meteorology analysis — all pure-logic paths, no MongoDB.

All classes are @pytest.mark.unit: no external dependencies, no network,
no MongoDB connections.
"""

from __future__ import annotations

import glob
import os
import sys
import tempfile
import types
from unittest.mock import MagicMock

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Mock unavailable optional packages before any hera.simulations.openFoam
# import so the module-level imports in toolkit.py don't crash.
# ---------------------------------------------------------------------------
def _mock_pkg(name: str) -> types.ModuleType:
    mod = types.ModuleType(name)
    mod.__path__ = []
    mod.__package__ = name
    sys.modules[name] = mod
    return mod


for _pkg in ["PyFoam", "PyFoam.RunDictionary", "PyFoam.Basics", "paraview", "paraview.simple"]:
    if _pkg not in sys.modules:
        _mock_pkg(_pkg)
for _sub in [
    "PyFoam.RunDictionary.ParsedParameterFile",
    "PyFoam.RunDictionary.BoundaryDict",
    "PyFoam.Basics.DataStructures",
]:
    if _sub not in sys.modules:
        sys.modules[_sub] = MagicMock()

for _ext in ["evtk", "evtk.hl"]:
    if _ext not in sys.modules:
        sys.modules[_ext] = MagicMock()

if "dask.distributed" not in sys.modules:
    sys.modules["dask.distributed"] = MagicMock()

# hermes: add repo path then try to import; fall back to mock
_HERMES_PATH = os.path.join(os.path.dirname(__file__), "..", "..", "Hermes")
if os.path.isdir(_HERMES_PATH) and _HERMES_PATH not in sys.path:
    sys.path.insert(0, os.path.abspath(_HERMES_PATH))

# dask.dataframe must be imported before analysis.addDatesColumns is called
import dask.dataframe  # noqa: E402 — order matters

# ---------------------------------------------------------------------------
# Imports from hera (safe after mocking)
# ---------------------------------------------------------------------------
from hera.measurements.GIS.raster.topography import TopographyToolkit, topographyAnalysis  # noqa: E402
from hera.measurements.meteorology.lowfreqdata.analysis import analysis as MeteoAnalysis  # noqa: E402
from hera.simulations.LSM.template import LSMTemplate  # noqa: E402
from hera.simulations.gaussian.Sigma import BriggsRural  # noqa: E402
from hera.simulations.openFoam.toolkit import OFToolkit, Presentation  # noqa: E402

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
WSG84 = 4326
ITM = 2039


# ===========================================================================
# 1. TopographyToolkit — convertPointsCRS
# ===========================================================================

@pytest.mark.unit
class TestTopographyCRS:
    """convertPointsCRS: pure geopandas CRS transformation, no MongoDB."""

    def test_returns_geodataframe(self):
        import geopandas as gpd
        result = TopographyToolkit.convertPointsCRS(None, [(34.78, 32.08)], WSG84, ITM)
        assert isinstance(result, gpd.GeoDataFrame)
        assert "geometry" in result.columns

    def test_wgs84_to_itm_single_point(self):
        """Tel Aviv (lon=34.78, lat=32.08) should land in the ITM domain."""
        result = TopographyToolkit.convertPointsCRS(None, [(34.78, 32.08)], WSG84, ITM)
        pt = result.geometry.iloc[0]
        # ITM easting for Israel is roughly 100,000 – 300,000 m
        assert 100_000 < pt.x < 300_000, f"ITM easting out of range: {pt.x}"
        # ITM northing for Israel is roughly 400,000 – 800,000 m
        assert 400_000 < pt.y < 800_000, f"ITM northing out of range: {pt.y}"

    def test_multiple_points_length_preserved(self):
        """All input points should appear in the output GeoDataFrame."""
        points = [(34.0, 31.0), (35.0, 32.0), (36.0, 33.0)]
        result = TopographyToolkit.convertPointsCRS(None, points, WSG84, ITM)
        assert len(result) == len(points)

    def test_roundtrip_itm_to_wgs84(self):
        """Converting ITM → WGS84 then WGS84 → ITM must recover original coords."""
        original = [(200_000.0, 650_000.0)]
        in_wgs = TopographyToolkit.convertPointsCRS(None, original, ITM, WSG84)
        back = TopographyToolkit.convertPointsCRS(
            None,
            [(in_wgs.geometry.iloc[0].x, in_wgs.geometry.iloc[0].y)],
            WSG84,
            ITM,
        )
        recovered = back.geometry.iloc[0]
        assert abs(recovered.x - original[0][0]) < 1.0
        assert abs(recovered.y - original[0][1]) < 1.0

    def test_numpy_array_input(self):
        """2-D numpy array of shape (N, 2) must be accepted."""
        pts = np.array([[34.78, 32.08], [35.0, 32.5]])
        result = TopographyToolkit.convertPointsCRS(None, pts, WSG84, ITM)
        assert len(result) == 2


# ===========================================================================
# 2. TopographyToolkit — analysis.calculateStastics
# ===========================================================================

@pytest.mark.unit
class TestTopographyStatistics:
    """calculateStastics: pure numpy statistics over a tabular elevation dataset."""

    @pytest.fixture
    def elevation_df(self):
        return pd.DataFrame(
            {
                "X": [1.0, 2.0, 3.0, 1.0, 2.0, 3.0],
                "Y": [10.0, 10.0, 10.0, 20.0, 20.0, 20.0],
                "Elevation": [100.0, 150.0, 200.0, 110.0, 160.0, 210.0],
            }
        )

    def test_all_keys_present(self, elevation_df):
        expected_keys = {
            "xmin", "xmax", "ymin", "ymax", "size",
            "mean", "std", "domainmax", "domainmaxlocation",
            "domainmin", "domainminlocation",
        }
        topo_analysis = topographyAnalysis(None)
        stats = topo_analysis.calculateStastics(elevation_df)
        assert set(stats.keys()) == expected_keys

    def test_mean_is_correct(self, elevation_df):
        topo_analysis = topographyAnalysis(None)
        stats = topo_analysis.calculateStastics(elevation_df)
        expected_mean = np.mean([100, 150, 200, 110, 160, 210])
        assert abs(stats["mean"] - expected_mean) < 1e-9

    def test_domain_size(self, elevation_df):
        topo_analysis = topographyAnalysis(None)
        stats = topo_analysis.calculateStastics(elevation_df)
        expected = (3.0 - 1.0) * (20.0 - 10.0)
        assert abs(stats["size"] - expected) < 1e-9

    def test_max_min_values(self, elevation_df):
        topo_analysis = topographyAnalysis(None)
        stats = topo_analysis.calculateStastics(elevation_df)
        assert stats["domainmax"] == 210.0
        assert stats["domainmin"] == 100.0

    def test_max_location_correct(self, elevation_df):
        """The max elevation (210) is at X=3, Y=20."""
        topo_analysis = topographyAnalysis(None)
        stats = topo_analysis.calculateStastics(elevation_df)
        xloc, yloc = stats["domainmaxlocation"]
        assert xloc == 3.0
        assert yloc == 20.0


# ===========================================================================
# 3. LSMTemplate — prepareParams
# ===========================================================================

@pytest.mark.unit
class TestLSMPrepareParams:
    """prepareParams: static method, pure unit-handling with pint."""

    def test_passthrough_when_desc_is_none(self):
        """No units desc → params returned unchanged (plain scalars)."""
        params = {"windspeed": 5.0, "height": 10.0, "name": "test"}
        result = LSMTemplate.prepareParams(None, dict(params))
        assert result["windspeed"] == 5.0
        assert result["height"] == 10.0
        assert result["name"] == "test"

    def test_duration_stays_numeric(self):
        """Duration key is wrapped as minutes then stripped back to magnitude."""
        params = {"duration": 30}
        result = LSMTemplate.prepareParams({"units": {"duration": "minutes"}}, dict(params))
        # stripConfigurationUnits returns the magnitude, duration=30 (minutes) is not
        # re-standardised because "duration" is in ignoreStandardization
        assert isinstance(result["duration"], (int, float))
        assert result["duration"] == 30

    def test_topoxn_cast_to_int(self):
        """TopoXn/TopoYn float values are silently cast to int."""
        params = {"TopoXn": 50.9, "TopoYn": 30.1}
        result = LSMTemplate.prepareParams(None, dict(params))
        assert isinstance(result["TopoXn"], int)
        assert result["TopoXn"] == 50
        assert isinstance(result["TopoYn"], int)
        assert result["TopoYn"] == 30

    def test_unrecognised_keys_preserved(self):
        """Keys not in desc['units'] pass through unmodified."""
        params = {"duration": 10, "extra_field": "hello"}
        result = LSMTemplate.prepareParams({"units": {"duration": "minutes"}}, dict(params))
        assert result["extra_field"] == "hello"


# ===========================================================================
# 4. OFToolkit — processorList
# ===========================================================================

@pytest.mark.unit
class TestOFProcessorList:
    """processorList: filesystem glob — processor* directories only."""

    def test_detects_processor_dirs(self, tmp_path):
        (tmp_path / "processor0").mkdir()
        (tmp_path / "processor1").mkdir()
        (tmp_path / "processor2").mkdir()
        (tmp_path / "constant").mkdir()      # must NOT appear
        (tmp_path / "0").mkdir()             # must NOT appear
        result = sorted(OFToolkit.processorList(None, str(tmp_path)))
        assert result == ["processor0", "processor1", "processor2"]

    def test_empty_case_dir(self, tmp_path):
        result = OFToolkit.processorList(None, str(tmp_path))
        assert result == []

    def test_returns_basenames_only(self, tmp_path):
        (tmp_path / "processor0").mkdir()
        result = OFToolkit.processorList(None, str(tmp_path))
        for name in result:
            assert os.sep not in name, f"Expected basename but got path: {name}"


# ===========================================================================
# 5. OFToolkit — Presentation.to_paraview_CSV
# ===========================================================================

@pytest.mark.unit
class TestOFParaviewCSV:
    """to_paraview_CSV: writes per-timestep CSV files from a particles DataFrame."""

    @pytest.fixture
    def presentation(self):
        p = Presentation.__new__(Presentation)
        p.datalayer = MagicMock()
        p.analysis = MagicMock()
        return p

    @pytest.fixture
    def particles_df(self):
        return pd.DataFrame(
            {
                "time": [1.0, 1.0, 2.0, 2.0],
                "globalX": [1.0, 2.0, 3.0, 4.0],
                "globalY": [5.0, 6.0, 7.0, 8.0],
                "globalZ": [9.0, 10.0, 11.0, 12.0],
            }
        )

    def test_creates_one_file_per_timestep(self, tmp_path, presentation, particles_df):
        presentation.to_paraview_CSV(particles_df, str(tmp_path), "particles")
        files = sorted(os.listdir(str(tmp_path)))
        assert len(files) == 2  # time 1.0 and time 2.0

    def test_filenames_contain_base(self, tmp_path, presentation, particles_df):
        presentation.to_paraview_CSV(particles_df, str(tmp_path), "mycloud")
        files = os.listdir(str(tmp_path))
        assert all("mycloud" in f for f in files)

    def test_csv_contains_correct_columns(self, tmp_path, presentation, particles_df):
        presentation.to_paraview_CSV(particles_df, str(tmp_path), "pts")
        for f in sorted(os.listdir(str(tmp_path))):
            content = (tmp_path / f).read_text()
            first_line = content.splitlines()[0]
            assert "globalX" in first_line
            assert "globalY" in first_line
            assert "globalZ" in first_line

    def test_each_csv_has_correct_row_count(self, tmp_path, presentation, particles_df):
        """Each timestep has 2 particles → each CSV should have 2 data rows + header."""
        presentation.to_paraview_CSV(particles_df, str(tmp_path), "pts")
        for f in sorted(os.listdir(str(tmp_path))):
            lines = (tmp_path / f).read_text().strip().splitlines()
            assert len(lines) == 3  # 1 header + 2 data rows


# ===========================================================================
# 6. BriggsRural — getSigma
# ===========================================================================

@pytest.mark.unit
class TestBriggsRuralSigma:
    """BriggsRural.getSigma: pure Briggs dispersion formula, no external deps."""

    @pytest.fixture
    def briggs(self):
        return BriggsRural()

    def test_returns_all_required_keys(self, briggs):
        result = briggs.getSigma(1000, "A", units=False)
        assert set(result.keys()) == {"sigmaX", "sigmaY", "sigmaZ", "distance"}

    def test_sigma_y_stability_a_at_1000m(self, briggs):
        """
        Stability A: A=0.22, B=1e-4, C=-0.5
        sigmaY = 0.22 * 1000 * (1 + 1e-4 * 1000)^(-0.5)
        """
        result = briggs.getSigma(1000, "A", units=False)
        expected = 0.22 * 1000 * (1 + 1e-4 * 1000) ** (-0.5)
        assert abs(result["sigmaY"][0] - expected) < 0.01

    def test_all_stability_classes_return_positive_sigmas(self, briggs):
        for stab in ["A", "B", "C", "D", "E", "F"]:
            result = briggs.getSigma(500, stab, units=False)
            assert result["sigmaY"][0] > 0, f"sigmaY <= 0 for stability {stab}"
            assert result["sigmaZ"][0] > 0, f"sigmaZ <= 0 for stability {stab}"

    def test_sigma_increases_with_distance(self, briggs):
        """Dispersion sigma must grow monotonically with downwind distance."""
        distances = [100, 500, 1000, 5000]
        result = briggs.getSigma(distances, "D", units=False)
        sigmas = list(result["sigmaY"])
        assert sigmas == sorted(sigmas), f"sigmaY not monotonic: {sigmas}"

    def test_units_flag_returns_pint_quantities(self, briggs):
        from pint import Quantity
        result = briggs.getSigma(1000, "B", units=True)
        assert isinstance(result["sigmaY"], Quantity)

    def test_distance_column_matches_input(self, briggs):
        x = [200.0, 400.0]
        result = briggs.getSigma(x, "C", units=False)
        np.testing.assert_array_almost_equal(result["distance"], x)


# ===========================================================================
# 7. LowFreq meteorology — addDatesColumns
# ===========================================================================

@pytest.mark.unit
class TestAddDatesColumns:
    """addDatesColumns: pure pandas date feature engineering, no MongoDB."""

    @pytest.fixture
    def meteo_analysis(self):
        return MeteoAnalysis(None)

    @pytest.fixture
    def hourly_df(self):
        return pd.DataFrame({"ts": pd.date_range("2023-07-15 06:00", periods=4, freq="h")})

    def test_output_contains_date_columns(self, meteo_analysis, hourly_df):
        result = meteo_analysis.addDatesColumns(hourly_df, datecolumn="ts")
        for col in ["yearonly", "monthonly", "dayonly", "season"]:
            assert col in result.columns, f"Missing column: {col}"

    def test_year_extracted_correctly(self, meteo_analysis, hourly_df):
        result = meteo_analysis.addDatesColumns(hourly_df, datecolumn="ts")
        assert (result["yearonly"] == 2023).all()

    def test_month_extracted_correctly(self, meteo_analysis, hourly_df):
        result = meteo_analysis.addDatesColumns(hourly_df, datecolumn="ts")
        assert (result["monthonly"] == 7).all()

    def test_july_maps_to_summer(self, meteo_analysis, hourly_df):
        result = meteo_analysis.addDatesColumns(hourly_df, datecolumn="ts")
        assert (result["season"] == "Summer").all(), f"Expected Summer, got {result['season'].unique()}"

    def test_january_maps_to_winter(self, meteo_analysis):
        df = pd.DataFrame({"ts": pd.date_range("2023-01-10", periods=3, freq="D")})
        result = meteo_analysis.addDatesColumns(df, datecolumn="ts")
        assert (result["season"] == "Winter").all()

    def test_time_column_is_hhmm(self, meteo_analysis):
        """Time column should be hour*100 + minute."""
        df = pd.DataFrame({"ts": [pd.Timestamp("2023-06-01 14:30")]})
        result = meteo_analysis.addDatesColumns(df, datecolumn="ts")
        assert result["Time"].iloc[0] == 1430

    def test_index_used_when_no_datecolumn(self, meteo_analysis):
        """When datecolumn=None, the datetime index is used."""
        df = pd.DataFrame(
            {"value": [1, 2]},
            index=pd.DatetimeIndex(["2023-03-01", "2023-09-01"]),
        )
        result = meteo_analysis.addDatesColumns(df, datecolumn=None)
        assert "season" in result.columns
        seasons = result["season"].tolist()
        assert "Spring" in seasons
        assert "Autumn" in seasons
