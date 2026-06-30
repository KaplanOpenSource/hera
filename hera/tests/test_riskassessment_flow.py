"""
Tests for the hera risk assessment mock flow (no live database required).

Uses a synthetic Gaussian concentration field to verify that the core risk
assessment pipeline — Agent instantiation, toxic load accumulation, and injury
contour extraction — remains functional as dependencies evolve.

Issue: https://github.com/KaplanOpenSource/hera/issues/882
"""

import numpy as np
import pandas as pd
import pytest

from hera.utils import ureg
from hera.riskassessment.agents.Agents import Agent


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Two-severity Lognormal10 agent with TenBerge (n=1) calculator.
# n=1 is equivalent to Haber's law, keeping the math straightforward.
_AGENT_DESCRIPTOR = {
    "name": "mock_agent",
    "effectParameters": {"tenbergeCoefficient": 1},
    "effects": {
        "inhalation": {
            "type": "Lognormal10",
            "calculator": {"TenBerge": {}},
            "units": "mg/m**3",
            "parameters": {
                "type": "Lognormal10DoseResponse",
                "levels": ["Severe", "Moderate"],
                "parameters": {
                    "Severe":   {"TL_50": 400, "sigma": 0.5},
                    "Moderate": {"TL_50": 100, "sigma": 0.5},
                },
            },
        }
    },
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_concentration_field(C_max=50.0, grid_size=41, nt=10):
    """Return a synthetic xr.Dataset that mimics an LSM dispersion output.

    Generates a time-invariant Gaussian plume centred at (1000, 1000) m on a
    2 km × 2 km grid.  With C_max=50 mg/m³ and nt=10 min, the peak cumulative
    toxic load is 500, comfortably above both TL_50 values (100 and 400).
    """
    import xarray as xr

    x = np.linspace(0, 2000, grid_size)
    y = np.linspace(0, 2000, grid_size)
    X, Y = np.meshgrid(x, y, indexing="ij")
    sigma = 400.0  # m — wide enough for clear polygons at 50 m grid spacing
    C_spatial = C_max * np.exp(-(((X - 1000.0) ** 2 + (Y - 1000.0) ** 2) / (2 * sigma ** 2)))

    times = pd.date_range("2024-01-01 00:00", periods=nt, freq="1min")
    C_data = np.broadcast_to(C_spatial[np.newaxis], (nt, grid_size, grid_size)).copy()

    return xr.Dataset(
        {
            "C": xr.DataArray(
                C_data,
                dims=["datetime", "x", "y"],
                coords={"datetime": times, "x": x, "y": y},
            )
        },
        attrs={"C": 1 * ureg.mg / ureg.m**3, "dt": 1 * ureg.min},
    )


# ---------------------------------------------------------------------------
# TestAgentInstantiation
# ---------------------------------------------------------------------------

class TestAgentInstantiation:
    """Verify that Agent objects can be built from plain JSON descriptors."""

    def test_agent_created_from_descriptor(self):
        agent = Agent(_AGENT_DESCRIPTOR)
        assert agent.name == "mock_agent"

    def test_agent_exposes_inhalation_effect(self):
        agent = Agent(_AGENT_DESCRIPTOR)
        assert "inhalation" in agent.effectNames

    def test_injury_has_both_severity_levels(self):
        agent = Agent(_AGENT_DESCRIPTOR)
        effect = agent["inhalation"]
        assert "Severe" in effect.levelNames
        assert "Moderate" in effect.levelNames

    def test_tenberge_coefficient_default(self):
        agent = Agent(_AGENT_DESCRIPTOR)
        assert agent.tenbergeCoefficient == 1

    def test_tenberge_coefficient_setter_rebuilds_effects(self):
        agent = Agent(_AGENT_DESCRIPTOR)
        agent.tenbergeCoefficient = 2
        assert agent.tenbergeCoefficient == 2
        # Effect must still be accessible after rebuild
        assert "inhalation" in agent.effectNames


# ---------------------------------------------------------------------------
# TestInjuryLevelMath
# ---------------------------------------------------------------------------

class TestInjuryLevelMath:
    """Unit-test the dose-response math for individual InjuryLevel objects."""

    @pytest.fixture(scope="class")
    def moderate_level(self):
        return Agent(_AGENT_DESCRIPTOR)["inhalation"]["Moderate"]

    def test_percent_at_tl50_is_half(self, moderate_level):
        tl50 = moderate_level.TL_50.magnitude
        pct = moderate_level.getPercent(tl50)
        assert abs(pct - 0.5) < 0.01, f"Expected ~0.5, got {pct}"

    def test_percent_monotonically_increasing(self, moderate_level):
        loads = [10, 50, 100, 200, 500]
        percents = [moderate_level.getPercent(tl) for tl in loads]
        assert percents == sorted(percents)

    def test_get_toxic_load_at_50pct_equals_tl50(self, moderate_level):
        tl50_ref = moderate_level.TL_50.magnitude
        tl50_calc = moderate_level.getToxicLoad(0.50)
        assert abs(tl50_calc - tl50_ref) / tl50_ref < 0.01

    def test_severe_has_higher_tl50_than_moderate(self):
        agent = Agent(_AGENT_DESCRIPTOR)
        tl50_severe = agent["inhalation"]["Severe"].TL_50.magnitude
        tl50_moderate = agent["inhalation"]["Moderate"].TL_50.magnitude
        assert tl50_severe > tl50_moderate


# ---------------------------------------------------------------------------
# TestToxicLoadCalculation
# ---------------------------------------------------------------------------

class TestToxicLoadCalculation:
    """Verify the cumulative toxic load produced by the TenBerge calculator."""

    @pytest.fixture(scope="class")
    def concentration_ds(self):
        return _make_concentration_field(C_max=50.0, grid_size=41, nt=10)

    def test_peak_toxic_load_exceeds_both_tl50s(self, concentration_ds):
        effect = Agent(_AGENT_DESCRIPTOR)["inhalation"]
        tl = effect.calculateToxicLoads(
            concentration_ds, time="datetime", field="C"
        ).compute()
        peak = float(tl.isel(datetime=-1).max())
        assert peak > 400, f"Peak TL ({peak:.1f}) must exceed Severe TL_50 (400)"

    def test_toxic_load_non_decreasing_at_peak(self, concentration_ds):
        effect = Agent(_AGENT_DESCRIPTOR)["inhalation"]
        tl = effect.calculateToxicLoads(
            concentration_ds, time="datetime", field="C"
        ).compute()
        x_mid = float(concentration_ds.x[len(concentration_ds.x) // 2])
        y_mid = float(concentration_ds.y[len(concentration_ds.y) // 2])
        tl_centre = tl.sel(x=x_mid, y=y_mid, method="nearest").values
        assert np.all(np.diff(tl_centre) >= 0)

    def test_toxic_load_near_zero_at_grid_corners(self, concentration_ds):
        effect = Agent(_AGENT_DESCRIPTOR)["inhalation"]
        tl = effect.calculateToxicLoads(
            concentration_ds, time="datetime", field="C"
        ).compute()
        corner = float(tl.sel(x=0.0, y=0.0, method="nearest").isel(datetime=-1))
        # Corner is 1414 m from centre; with sigma=400 m, concentration is negligible
        assert corner < 5.0, f"Corner TL ({corner:.3f}) should be near zero"

    def test_higher_n_raises_peak_for_high_concentration(self):
        """For C > 1, TenBerge with n=2 must give a higher toxic load than n=1."""
        ds = _make_concentration_field(C_max=10.0, grid_size=21, nt=5)

        desc_n1 = dict(_AGENT_DESCRIPTOR)
        desc_n1 = {**_AGENT_DESCRIPTOR, "effectParameters": {"tenbergeCoefficient": 1}}
        desc_n2 = {**_AGENT_DESCRIPTOR, "effectParameters": {"tenbergeCoefficient": 2}}

        tl1 = Agent(desc_n1)["inhalation"].calculateToxicLoads(ds, time="datetime", field="C").compute()
        tl2 = Agent(desc_n2)["inhalation"].calculateToxicLoads(ds, time="datetime", field="C").compute()

        peak1 = float(tl1.isel(datetime=-1).max())
        peak2 = float(tl2.isel(datetime=-1).max())
        assert peak2 > peak1, f"n=2 peak ({peak2:.1f}) must exceed n=1 peak ({peak1:.1f})"


# ---------------------------------------------------------------------------
# TestRiskAssessmentFlow
# ---------------------------------------------------------------------------

class TestRiskAssessmentFlow:
    """Full end-to-end: concentration field → thresholdGeoDataFrame."""

    @pytest.fixture(scope="class")
    def injury_contours(self):
        pytest.importorskip("geopandas")
        pytest.importorskip("shapely")
        ds = _make_concentration_field(C_max=50.0, grid_size=41, nt=10)
        agent = Agent(_AGENT_DESCRIPTOR)
        return agent["inhalation"].calculateRegionOfInjured(
            ds, field="C", time="datetime"
        )

    def test_result_is_geodataframe(self, injury_contours):
        import geopandas
        assert isinstance(injury_contours, geopandas.GeoDataFrame)

    def test_result_not_empty(self, injury_contours):
        assert len(injury_contours) > 0

    def test_expected_columns_present(self, injury_contours):
        for col in ("severity", "ToxicLoad", "TotalPolygon", "ThresholdPolygon", "percentEffected"):
            assert col in injury_contours.columns, f"Missing column: {col}"

    def test_both_severity_levels_present(self, injury_contours):
        severities = set(injury_contours["severity"].unique())
        assert "Severe" in severities
        assert "Moderate" in severities

    def test_all_total_polygons_have_positive_area(self, injury_contours):
        areas = injury_contours.set_geometry("TotalPolygon")["TotalPolygon"].area
        assert (areas > 0).all()

    def test_moderate_max_area_larger_than_severe(self, injury_contours):
        """Moderate (lower TL_50=100) must reach a larger footprint than Severe (TL_50=400)."""
        by_severity = injury_contours.set_geometry("TotalPolygon").groupby("severity")
        severe_max = by_severity.get_group("Severe")["TotalPolygon"].area.max()
        moderate_max = by_severity.get_group("Moderate")["TotalPolygon"].area.max()
        assert moderate_max >= severe_max, (
            f"Moderate max area ({moderate_max:.0f} m²) should be ≥ Severe ({severe_max:.0f} m²)"
        )

    def test_percent_effected_in_range(self, injury_contours):
        pct = injury_contours["percentEffected"]
        assert (pct >= 0).all() and (pct <= 1).all()
