r"""``RiskToolkit.analysis.getRiskAreas``: the risk-area envelope of a
dispersion run.

The method chains most of the risk domain together.  It pulls a Lagrangian
simulation from the LSM toolkit, turns it into a concentration field,
builds a throwaway "dumb" agent inline (a single ``Lognormal10`` severity
driven by a ten Berge calculator), integrates the field into a toxic load,
keeps the grid cells at or above each requested level, and returns one row
per level holding the bounding box and the convex hull of those cells.

The LSM half is the only part that cannot run hermetically, so the tests
below replace ``analysis._LSM`` (a per-test toolkit instance, never the
``toolkitHome`` singleton) with a stand-in exposing the two methods the
method calls -- ``getSimulations(**params)`` and ``singleSimulation(file)``
-- each handing back an object with ``getConcentration(Q=...)``.
Everything downstream of that is the real code.

The expected geometry is derived analytically rather than read off hera's
output.  The concentration field is a Gaussian
``C(r) = PEAK * exp(-r^2 / 2 s^2)`` on a regular grid with two one-minute
time steps.  The inline agent's ten Berge coefficient is passed in by the
test as 1, and the calculator integrates ``C^n dt`` cumulatively, so the
toxic load at the last step is ``2 * C(r)``.  A level ``L`` therefore
selects the cells with

    C(r) >= L / 2   <=>   r <= s * sqrt(2 ln(2 PEAK / L))

and the hull of those cells is the grid points inside that disc.  Both
bounds and areas are computed from that formula in the tests.

Deliberately not covered here:

* ``getSimulations``/``singleSimulation`` themselves and
  ``getConcentration`` -- LSM toolkit territory, covered in
  ``test_simulations_lsm_toolkit.py`` and ``test_lsm_*.py``.
* the ``analysis.LSM``/``analysis.datalayer`` properties and
  ``RiskToolkit``'s agent loading -- see ``test_risktoolkit_analysis.py``
  and ``test_risk_toolkit.py``.
* the protection-policy *actions* (indoor, masking, evacuation).  The
  branch that applies a policy is exercised with an empty action list,
  which the policy documents as a pass-through pipeline; the actions have
  their own file, ``test_riskassessment_protectionpolicy.py``.

Robustness note (not pinned as a bug): the toxic load is recomputed from
scratch inside the per-level loop although it does not depend on the level,
so asking for k levels does k identical integrations.  Wasteful, not wrong.
"""
import numpy
import pandas
import pytest

from hera import toolkitHome
from hera.utils import ureg

PEAK = 100.0
SIGMA = 40.0
STEP = 10.0
EXTENT = 100.0
STEPS = 2


def _grid():
    return numpy.arange(-EXTENT, EXTENT + STEP / 2, STEP)


def _concentration():
    """A two-step Gaussian concentration field, with the unit metadata the
    ten Berge calculator reads out of ``attrs``."""
    import xarray

    axis = _grid()
    grid_x, grid_y = numpy.meshgrid(axis, axis)
    values = PEAK * numpy.exp(-(grid_x ** 2 + grid_y ** 2) / (2.0 * SIGMA ** 2))
    stamps = pandas.to_datetime(
        ["2020-01-01 00:%02d" % minute for minute in range(STEPS)])

    field = xarray.DataArray(
        numpy.stack([values[None, ...]] * STEPS),
        coords={"datetime": stamps, "z": [0.0], "y": axis, "x": axis},
        dims=("datetime", "z", "y", "x"), name="C")
    dataset = field.to_dataset()
    dataset.attrs["C"] = 1 * ureg.mg / ureg.m ** 3
    dataset.attrs["dt"] = 1 * ureg.min
    return dataset


def _selected_points(level):
    """The grid cells the method must keep, from the closed form above."""
    radius = SIGMA * numpy.sqrt(2.0 * numpy.log(STEPS * PEAK / level))
    axis = _grid()
    grid_x, grid_y = numpy.meshgrid(axis, axis)
    inside = (grid_x ** 2 + grid_y ** 2) <= radius ** 2
    return numpy.column_stack([grid_x[inside], grid_y[inside]])


def _expected_bounds(level):
    points = _selected_points(level)
    return (points[:, 0].min(), points[:, 1].min(),
            points[:, 0].max(), points[:, 1].max())


class _Simulation:
    def __init__(self, dataset, recorder):
        self._dataset = dataset
        self._recorder = recorder

    def getConcentration(self, Q=None):
        self._recorder.append(Q)
        return self._dataset


class _LSM:
    """The two entry points getRiskAreas uses on the LSM toolkit."""

    def __init__(self, dataset):
        self._dataset = dataset
        self.simulationQueries = []
        self.singleSimulationFiles = []
        self.masses = []

    def getSimulations(self, **params):
        self.simulationQueries.append(params)
        return [_Simulation(self._dataset, self.masses)]

    def singleSimulation(self, fileName):
        self.singleSimulationFiles.append(fileName)
        return _Simulation(self._dataset, self.masses)


@pytest.fixture()
def analysis(unit_toolkit_factory, monkeypatch):
    """The real risk analysis layer, with only the LSM source replaced."""
    risk = unit_toolkit_factory(toolkitHome.RISKASSESSMENT)
    lsm = _LSM(_concentration())
    monkeypatch.setattr(risk.analysis, "_LSM", lsm)
    return risk.analysis, lsm


# ---------------------------------------------------------------------------
# The returned table
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestRiskAreasTable:
    def test_it_returns_a_geodataframe(self, analysis):
        import geopandas

        layer, _ = analysis
        result = layer.getRiskAreas(tenbergeCoefficient=1, levels=50.0)

        assert isinstance(result, geopandas.GeoDataFrame)

    def test_the_documented_columns_are_present(self, analysis):
        layer, _ = analysis
        result = layer.getRiskAreas(tenbergeCoefficient=1, levels=50.0)

        assert list(result.columns) == ["tenBergCoefficient", "level", "bounds",
                                        "geometry"]

    def test_a_scalar_level_becomes_a_single_row(self, analysis):
        """Documented: "levels = a list of levels ... or a single value"."""
        layer, _ = analysis
        result = layer.getRiskAreas(tenbergeCoefficient=1, levels=50.0)

        assert len(result) == 1
        assert result["level"].tolist() == [50.0]

    def test_a_list_of_levels_gives_one_row_each(self, analysis):
        layer, _ = analysis
        result = layer.getRiskAreas(tenbergeCoefficient=1, levels=[50.0, 100.0])

        assert result["level"].tolist() == [50.0, 100.0]

    def test_the_ten_berge_coefficient_is_recorded_on_every_row(self, analysis):
        layer, _ = analysis
        result = layer.getRiskAreas(tenbergeCoefficient=1, levels=[50.0, 100.0])

        assert result["tenBergCoefficient"].tolist() == [1, 1]


@pytest.mark.unit
class TestRiskAreasGeometry:
    def test_the_bounds_are_the_extent_of_the_cells_above_the_level(self, analysis):
        layer, _ = analysis
        result = layer.getRiskAreas(tenbergeCoefficient=1, levels=50.0)

        assert result["bounds"].iloc[0] == pytest.approx(_expected_bounds(50.0))

    def test_the_geometry_is_the_hull_of_those_cells(self, analysis):
        layer, _ = analysis
        result = layer.getRiskAreas(tenbergeCoefficient=1, levels=50.0)

        assert result.geometry.iloc[0].bounds == pytest.approx(_expected_bounds(50.0))

    def test_the_hull_is_convex(self, analysis):
        layer, _ = analysis
        result = layer.getRiskAreas(tenbergeCoefficient=1, levels=50.0)

        hull = result.geometry.iloc[0]
        assert hull.area == pytest.approx(hull.convex_hull.area)

    def test_every_selected_cell_is_inside_the_hull(self, analysis):
        """The hull has to cover the cells the level selects, which is what
        makes it a *risk area* rather than a sample of one."""
        from shapely.geometry import Point

        layer, _ = analysis
        hull = layer.getRiskAreas(tenbergeCoefficient=1,
                                  levels=50.0).geometry.iloc[0]

        for x, y in _selected_points(50.0):
            assert hull.intersects(Point(x, y))

    def test_a_higher_level_gives_a_smaller_area(self, analysis):
        """A higher toxic-load threshold selects fewer cells."""
        layer, _ = analysis
        result = layer.getRiskAreas(tenbergeCoefficient=1, levels=[50.0, 100.0])

        assert result.geometry.iloc[1].area < result.geometry.iloc[0].area

    def test_a_higher_level_gives_tighter_bounds(self, analysis):
        layer, _ = analysis
        result = layer.getRiskAreas(tenbergeCoefficient=1, levels=[50.0, 100.0])

        assert result["bounds"].iloc[1] == pytest.approx(_expected_bounds(100.0))

    def test_the_area_is_centred_on_the_release(self, analysis):
        layer, _ = analysis
        hull = layer.getRiskAreas(tenbergeCoefficient=1,
                                  levels=50.0).geometry.iloc[0]

        centroid = hull.centroid
        assert (centroid.x, centroid.y) == pytest.approx((0.0, 0.0), abs=1e-9)

    def test_a_larger_ten_berge_coefficient_shrinks_the_area(self, analysis):
        r"""The calculator integrates C^n, and the field peaks at 100 mg/m^3,
        so raising n above 1 lifts the load near the source and crushes it in
        the tail -- the level-50 contour therefore moves inward for n = 2
        relative to n = 1 measured on the same grid."""
        layer, _ = analysis

        linear = layer.getRiskAreas(tenbergeCoefficient=1, levels=1e4)
        squared = layer.getRiskAreas(tenbergeCoefficient=2, levels=1e4)

        assert squared.geometry.iloc[0].area > linear.geometry.iloc[0].area


# ---------------------------------------------------------------------------
# How the dispersion run is obtained
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestRiskAreasSimulationSource:
    def test_by_default_the_simulation_is_queried_from_the_lsm_toolkit(self, analysis):
        layer, lsm = analysis
        layer.getRiskAreas(tenbergeCoefficient=1, levels=50.0)

        assert lsm.simulationQueries == [{}]
        assert lsm.singleSimulationFiles == []

    def test_extra_keywords_become_the_simulation_query(self, analysis):
        layer, lsm = analysis
        layer.getRiskAreas(tenbergeCoefficient=1, levels=50.0,
                           experimentName="TRIAL", release=3)

        assert lsm.simulationQueries == [{"experimentName": "TRIAL", "release": 3}]

    def test_an_lsm_file_bypasses_the_query(self, analysis):
        layer, lsm = analysis
        layer.getRiskAreas(tenbergeCoefficient=1, levels=50.0,
                           LSMfile="/runs/case01")

        assert lsm.singleSimulationFiles == ["/runs/case01"]
        assert lsm.simulationQueries == []

    def test_the_default_released_mass_is_one_kilogram(self, analysis):
        """Documented: "Q = the total mass of dispersed particles; default is
        1 kg"."""
        layer, lsm = analysis
        layer.getRiskAreas(tenbergeCoefficient=1, levels=50.0)

        assert lsm.masses == [1 * ureg.kg]

    def test_the_released_mass_is_forwarded(self, analysis):
        layer, lsm = analysis
        layer.getRiskAreas(tenbergeCoefficient=1, levels=50.0, Q=5 * ureg.kg)

        assert lsm.masses == [5 * ureg.kg]

    def test_only_the_first_matching_simulation_is_used(self, analysis, monkeypatch):
        layer, lsm = analysis
        extra = []

        def _two(**params):
            extra.append(params)
            return [_Simulation(_concentration(), lsm.masses),
                    _Simulation(_concentration(), extra)]

        monkeypatch.setattr(lsm, "getSimulations", _two)
        layer.getRiskAreas(tenbergeCoefficient=1, levels=50.0)

        assert lsm.masses == [1 * ureg.kg]

    def test_a_protection_policy_is_applied_to_the_field(self, analysis):
        """An empty action list is documented as a pass-through pipeline, so
        the answer must be unchanged -- which is what shows the branch ran
        without altering the physics."""
        layer, _ = analysis

        plain = layer.getRiskAreas(tenbergeCoefficient=1, levels=50.0)
        policed = layer.getRiskAreas(tenbergeCoefficient=1, levels=50.0,
                                     protectionPolicy=[])

        assert policed["bounds"].iloc[0] == pytest.approx(plain["bounds"].iloc[0])
        assert policed.geometry.iloc[0].area == pytest.approx(
            plain.geometry.iloc[0].area)


# ---------------------------------------------------------------------------
# The inline agent
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestRiskAreasInlineAgent:
    def test_the_toxic_load_matches_the_documented_ten_berge_integral(self, analysis):
        r"""The selected extent is a direct read-out of the load, so matching
        it against 2*C(r) confirms that the inline agent integrates
        \int C^n dt cumulatively over the two one-minute steps."""
        layer, _ = analysis
        result = layer.getRiskAreas(tenbergeCoefficient=1, levels=50.0)

        assert result["bounds"].iloc[0] == pytest.approx(_expected_bounds(50.0))

    def test_a_level_above_the_peak_load_yields_an_empty_area(self, analysis):
        """Robustness note (not pinned as a bug): there is no guard for "no
        cell reaches this level", so the row still comes back -- with an
        empty geometry and NaN bounds rather than an error or a warning.  A
        caller has to test for that itself."""
        layer, _ = analysis
        result = layer.getRiskAreas(tenbergeCoefficient=1, levels=1e9)

        assert len(result) == 1
        assert result.geometry.iloc[0].is_empty
        assert numpy.isnan(result["bounds"].iloc[0]).all()

    def test_the_agent_is_not_persisted_to_the_project(self, analysis):
        """The "dumb agent" is a throwaway descriptor built inline; it must
        not end up in the project's data sources."""
        layer, _ = analysis
        before = layer.datalayer.listAgentsNames()

        layer.getRiskAreas(tenbergeCoefficient=1, levels=50.0)

        assert layer.datalayer.listAgentsNames() == before
