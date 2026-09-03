"""``riskAreaAlgorithm_Sweep``'s actual sweep: the bounding-box mesh, the
per-release-point calculation, and the ``calculate`` driver.

The algorithm answers "if the release happened *here*, how many casualties
would there be?" for every point on a regular mesh.  ``_findBoundingBox``
builds that mesh, ``_doCalculation`` evaluates one point, and ``calculate``
walks the mesh (serially, or through ``multiprocessing`` when
``parallel=True``).

``_findBoundingBox``'s geometry is derivable in closed form, and all the
expected numbers below are derived that way rather than read off hera's
output.  For a demography whose convex hull spans ``[0, 300]^2`` and a
mathematical angle of 0 (no rotation), with the largest-area severity's
total polygon spanning ``x in [0, 200]``, ``y in [-100, 100]``:

    length = max(x) = 200            width = max(y) - min(y) = 200
    minX = 0   - (length + outlayers*dxdy)
    maxX = 300 + outlayers*dxdy
    minY = 0   - (width  + outlayers*dxdy)
    maxY = 300 + width + outlayers*dxdy

and the mesh is ``numpy.arange(minX, maxX, dxdy) x
numpy.arange(minY, maxY, dxdy)``.  With ``dxdy=150`` and ``outlayers=1``
that is 6 x 7 = 42 points.  The asymmetry (``maxY`` gets an extra
``width``, ``maxX`` does not) is the algorithm's own, and is asserted as
such.

Two bugs are pinned below.  B263: ``_doCalculation`` calls a method that
does not exist, so no sweep can complete.  B266: ``_findBoundingBox``
cannot even size its mesh when the final time step has a single isopleth
row, which is the normal shape for a single-severity injury.
``calculate``'s parallel branch is
therefore only reachable up to the same failure; it is not exercised here
because ``multiprocessing.Pool`` inside a hermetic unit test would fork the
mongomock-patched interpreter, and the serial branch already covers the
same driver logic.

Deliberately not covered here: ``getRiskAreaAlgorithm`` dispatch, the
constructor defaults and the ``dxdy``/``parallel``/``workerCount``
property setters -- see ``test_risk_riskareas.py``.
"""
import numpy
import pandas
import pytest
from shapely.geometry import Polygon

from hera.riskassessment.agents.effects.thresholdGeoDataFrame import thresholdGeoDataFrame
from hera.riskassessment.analysis.riskAreas import riskAreaAlgorithm_Sweep

DXDY = 150.0
OUTLAYERS = 1

# The demography hull and the largest total polygon, as described above.
DEMOG_SPAN = 300.0
ISOPLETH_LENGTH = 200.0
ISOPLETH_WIDTH = 200.0

EXPECTED_MIN_X = -(ISOPLETH_LENGTH + OUTLAYERS * DXDY)
EXPECTED_MAX_X = DEMOG_SPAN + OUTLAYERS * DXDY
EXPECTED_MIN_Y = -(ISOPLETH_WIDTH + OUTLAYERS * DXDY)
EXPECTED_MAX_Y = DEMOG_SPAN + ISOPLETH_WIDTH + OUTLAYERS * DXDY


def _rect(x0, y0, width, height):
    return Polygon([(x0, y0), (x0 + width, y0), (x0 + width, y0 + height),
                    (x0, y0 + height)])


def _isopleths(timestamps=("2020-01-01 00:00",)):
    """Two severities per time step; 'Light' always has the larger area."""
    rows = []
    for stamp in timestamps:
        rows.append({"severity": "Severe", "datetime": stamp,
                     "TotalPolygon": _rect(0, -50, 100, 100),
                     "ThresholdPolygon": _rect(0, -50, 100, 100),
                     "percentEffected": 1.0, "ToxicLoad": 10.0})
        rows.append({"severity": "Light", "datetime": stamp,
                     "TotalPolygon": _rect(0, -ISOPLETH_WIDTH / 2,
                                           ISOPLETH_LENGTH, ISOPLETH_WIDTH),
                     "ThresholdPolygon": _rect(0, -ISOPLETH_WIDTH / 2,
                                               ISOPLETH_LENGTH, ISOPLETH_WIDTH),
                     "percentEffected": 0.5, "ToxicLoad": 5.0})
    frame = pandas.DataFrame(rows)
    frame["datetime"] = pandas.to_datetime(frame["datetime"])
    return thresholdGeoDataFrame(frame, geometry="ThresholdPolygon")


def _demography():
    import geopandas

    return geopandas.GeoDataFrame(
        {"total_pop": [100.0], "geometry": [_rect(0, 0, DEMOG_SPAN, DEMOG_SPAN)]},
        geometry="geometry", crs=2039)


@pytest.fixture()
def sweep():
    return riskAreaAlgorithm_Sweep(dxdy=DXDY, outlayers=OUTLAYERS, parallel=False)


# ---------------------------------------------------------------------------
# _findBoundingBox
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestFindBoundingBox:
    @pytest.fixture()
    def mesh(self, sweep):
        return sweep._findBoundingBox(_isopleths(), _demography(),
                                      mathematical_angle=0.0)

    def test_the_mesh_spans_the_derived_bounding_box(self, mesh):
        minX, minY, maxX, maxY = mesh.total_bounds
        expectedXs = numpy.arange(EXPECTED_MIN_X, EXPECTED_MAX_X, DXDY)
        expectedYs = numpy.arange(EXPECTED_MIN_Y, EXPECTED_MAX_Y, DXDY)

        assert (minX, maxX) == (expectedXs[0], expectedXs[-1])
        assert (minY, maxY) == (expectedYs[0], expectedYs[-1])

    def test_the_mesh_has_one_point_per_grid_cell(self, mesh):
        expected = (len(numpy.arange(EXPECTED_MIN_X, EXPECTED_MAX_X, DXDY))
                    * len(numpy.arange(EXPECTED_MIN_Y, EXPECTED_MAX_Y, DXDY)))

        assert len(mesh) == expected

    def test_the_spacing_between_neighbours_is_dxdy(self, mesh):
        ys = sorted({point.y for point in mesh})

        assert numpy.diff(ys) == pytest.approx(DXDY)

    def test_the_upwind_margin_covers_the_isopleth_length(self, mesh):
        """A release has to be allowed far enough upwind that its plume can
        still reach the demography, hence the extra `length` on minX."""
        assert mesh.total_bounds[0] <= -ISOPLETH_LENGTH

    def test_the_downwind_margin_is_only_the_outlayers(self, mesh):
        """Characterisation of the algorithm's asymmetry: maxX gets
        outlayers*dxdy, while maxY additionally gets the plume width."""
        maxX = mesh.total_bounds[2]
        maxY = mesh.total_bounds[3]

        assert maxX < DEMOG_SPAN + ISOPLETH_WIDTH
        assert maxY >= DEMOG_SPAN + ISOPLETH_WIDTH - DXDY

    def test_more_outlayers_grow_the_mesh(self, sweep):
        wider = riskAreaAlgorithm_Sweep(dxdy=DXDY, outlayers=5, parallel=False)

        narrow = sweep._findBoundingBox(_isopleths(), _demography(),
                                        mathematical_angle=0.0)
        broad = wider._findBoundingBox(_isopleths(), _demography(),
                                       mathematical_angle=0.0)

        assert len(broad) > len(narrow)
        assert broad.total_bounds[0] < narrow.total_bounds[0]

    def test_a_finer_spacing_gives_more_points(self):
        coarse = riskAreaAlgorithm_Sweep(dxdy=300, outlayers=0, parallel=False)
        fine = riskAreaAlgorithm_Sweep(dxdy=150, outlayers=0, parallel=False)

        assert len(fine._findBoundingBox(_isopleths(), _demography(),
                                         mathematical_angle=0.0)) > \
            len(coarse._findBoundingBox(_isopleths(), _demography(),
                                        mathematical_angle=0.0))

    def test_the_largest_severity_sets_the_mesh_extent(self, sweep):
        """Documented behaviour: the isopleth used for sizing is the severity
        with the largest dissolved area.  Shrinking the *smaller* severity
        must therefore leave the mesh alone."""
        shrunk = _isopleths()
        shrunk.loc[shrunk["severity"] == "Severe", "TotalPolygon"] = _rect(0, -5, 10, 10)

        assert sweep._findBoundingBox(shrunk, _demography(),
                                      mathematical_angle=0.0).total_bounds.tolist() == \
            sweep._findBoundingBox(_isopleths(), _demography(),
                                   mathematical_angle=0.0).total_bounds.tolist()

    def test_only_the_last_time_step_sizes_the_mesh(self, sweep):
        """The plume is largest at the final time step, so that is the one the
        bounding box is built from -- adding earlier steps must not change it.
        """
        oneStep = sweep._findBoundingBox(_isopleths(("2020-01-01 00:05",)),
                                         _demography(), mathematical_angle=0.0)
        twoSteps = sweep._findBoundingBox(
            _isopleths(("2020-01-01 00:00", "2020-01-01 00:05")),
            _demography(), mathematical_angle=0.0)

        assert oneStep.total_bounds.tolist() == twoSteps.total_bounds.tolist()

    def test_rotating_the_wind_rotates_the_mesh(self, sweep):
        """The mesh is built in the wind frame and rotated back, so a 90 deg
        wind gives a mesh of the same size in a different orientation."""
        straight = sweep._findBoundingBox(_isopleths(), _demography(),
                                          mathematical_angle=0.0)
        turned = sweep._findBoundingBox(_isopleths(), _demography(),
                                        mathematical_angle=90.0)

        assert len(turned) == len(straight)
        assert turned.total_bounds.tolist() != straight.total_bounds.tolist()

    def test_a_full_turn_returns_the_same_mesh(self, sweep):
        straight = sweep._findBoundingBox(_isopleths(), _demography(),
                                          mathematical_angle=0.0)
        turned = sweep._findBoundingBox(_isopleths(), _demography(),
                                        mathematical_angle=360.0)

        assert turned.total_bounds == pytest.approx(straight.total_bounds)

    def test_it_returns_a_geoseries_not_the_documented_geodataframe(self, mesh):
        """Robustness note (not pinned as a bug): the docstring promises a
        GeoDataFrame, but the final `.rotate(...)` returns a GeoSeries.  Every
        in-tree caller only iterates it and reads ``.x``/``.y``, which works
        on both, so this is a documentation defect rather than a functional
        one -- recorded here so a future change of return type is a decision.
        """
        import geopandas

        assert isinstance(mesh, geopandas.GeoSeries)
        assert not isinstance(mesh, geopandas.GeoDataFrame)

    def test_the_severity_column_name_is_configurable(self, sweep):
        renamed = _isopleths().rename(columns={"severity": "level"})
        renamed = thresholdGeoDataFrame(renamed, geometry="ThresholdPolygon")

        mesh = sweep._findBoundingBox(renamed, _demography(),
                                      mathematical_angle=0.0,
                                      severityColumn="level")

        assert len(mesh) == len(sweep._findBoundingBox(_isopleths(), _demography(),
                                                       mathematical_angle=0.0))


# ---------------------------------------------------------------------------
# _doCalculation / calculate
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestASingleIsoplethRowBreaksTheMesh:
    """B266: see the reasons below."""

    @staticmethod
    def _single():
        """One severity, one time step -- an injury with a single level."""
        frame = _isopleths()
        frame = frame[frame["severity"] == "Light"]
        return thresholdGeoDataFrame(frame, geometry="ThresholdPolygon")

    REASON = (
        "B266: _findBoundingBox narrows the isopleths to the last time step "
        "with effectIsopleths.set_index('datetime').loc[maxdatetime].  When "
        "only one row carries that timestamp, pandas' .loc on a scalar key "
        "returns a Series rather than a one-row frame, and the next line "
        "(bounds_effectIsopleths.set_geometry(geometryColumn)) raises "
        "AttributeError: 'Series' object has no attribute 'set_geometry'.  "
        "One row per final time step is the normal shape for an injury with a "
        "single severity level, so the sweep cannot size its mesh for such an "
        "agent.  Selecting with .loc[[maxdatetime]] (or a boolean mask) keeps "
        "the frame.  See the consolidated findings issue."
    )

    @pytest.mark.xfail(strict=True, reason=REASON)
    def test_a_single_severity_should_still_size_the_mesh(self, sweep):
        mesh = sweep._findBoundingBox(self._single(), _demography(),
                                      mathematical_angle=0.0)

        expected = (len(numpy.arange(EXPECTED_MIN_X, EXPECTED_MAX_X, DXDY))
                    * len(numpy.arange(EXPECTED_MIN_Y, EXPECTED_MAX_Y, DXDY)))
        assert len(mesh) == expected

    def test_a_single_severity_currently_raises(self, sweep):
        """Characterisation of B266."""
        with pytest.raises(AttributeError,
                           match="'Series' object has no attribute 'set_geometry'"):
            sweep._findBoundingBox(self._single(), _demography(),
                                   mathematical_angle=0.0)

    def test_adding_a_second_row_at_the_same_time_step_makes_it_work(self, sweep):
        """Characterisation of B266: isolates the row count as the trigger.
        The extra row is a duplicate of the first, so the mesh it produces is
        the one the single-row case should have produced."""
        doubled = self._single()
        doubled = thresholdGeoDataFrame(
            pandas.concat([doubled, doubled], ignore_index=True),
            geometry="ThresholdPolygon")

        mesh = sweep._findBoundingBox(doubled, _demography(), mathematical_angle=0.0)

        expected = (len(numpy.arange(EXPECTED_MIN_X, EXPECTED_MAX_X, DXDY))
                    * len(numpy.arange(EXPECTED_MIN_Y, EXPECTED_MAX_Y, DXDY)))
        assert len(mesh) == expected

    def test_the_scalar_loc_selection_is_what_degrades(self):
        """Characterisation of B266's mechanism, derived from pandas rather
        than from hera."""
        frame = self._single().set_index("datetime")
        stamp = frame.index.max()

        assert isinstance(frame.loc[stamp], pandas.Series)
        assert isinstance(frame.loc[[stamp]], pandas.DataFrame)


@pytest.mark.unit
class TestSweepCannotRun:
    """B263: see the reasons below."""

    @staticmethod
    def _params():
        return {"effectIsopleths": _isopleths(),
                "demog": _demography(),
                "rotate_angle": 0.0,
                "valueColumn": "effectedtotal_pop"}

    REASON = (
        "B263: riskAreaAlgorithm_Sweep._doCalculation projects the "
        "isopleths with effectIsopleths.datalayer(demog, releaseLoc, "
        "mathematical_angle=...), but thresholdGeoDataFrame has no "
        "'datalayer' attribute -- its projection method is 'project', with "
        "exactly that signature.  Every release point therefore raises "
        "AttributeError: 'thresholdGeoDataFrame' object has no attribute "
        "'datalayer', so riskAreaAlgorithm_Sweep.calculate can never "
        "complete and the only implemented risk-area algorithm is unusable. "
        "See the consolidated findings issue."
    )

    @pytest.mark.xfail(strict=True, reason=REASON)
    def test_one_release_point_should_yield_casualties_per_severity(self, sweep):
        result = sweep._doCalculation((0.0, 0.0), self._params())[0]

        assert sorted(result["severity"].unique()) == ["Light", "Severe"]

    def test_one_release_point_currently_raises(self, sweep):
        """Characterisation of B263."""
        with pytest.raises(AttributeError, match="no attribute 'datalayer'"):
            sweep._doCalculation((0.0, 0.0), self._params())

    def test_the_projection_method_it_should_be_calling_does_exist(self):
        """Characterisation of B263: names the intended call site, so the
        fix is unambiguous."""
        assert not hasattr(_isopleths(), "datalayer")
        assert callable(_isopleths().project)

    @pytest.mark.xfail(strict=True, reason=REASON)
    def test_a_serial_sweep_should_return_one_row_per_point_and_severity(self, sweep):
        result = sweep.calculate(_isopleths(), _demography(), mathematical_angle=0.0)

        assert set(result.columns) >= {"x", "y", "severity", "datetime"}

    def test_a_serial_sweep_currently_raises_on_the_first_point(self, sweep, capsys):
        """Characterisation of B263: the driver itself is fine -- it reaches
        the first mesh point and prints its progress before dying."""
        with pytest.raises(AttributeError, match="no attribute 'datalayer'"):
            sweep.calculate(_isopleths(), _demography(), mathematical_angle=0.0)

        assert "Processing point 0 out of 42" in capsys.readouterr().out

    def test_the_mesh_is_built_before_the_failure(self, sweep, monkeypatch):
        """Characterisation of B263: _findBoundingBox is reached with the
        converted angle, which is why it is asserted rather than xfailed."""
        seen = {}
        original = type(sweep)._findBoundingBox

        def _spy(self, effectIsopleths, demog, mathematical_angle, **kwargs):
            seen["angle"] = mathematical_angle
            seen["kwargs"] = kwargs
            return original(self, effectIsopleths, demog, mathematical_angle, **kwargs)

        monkeypatch.setattr(type(sweep), "_findBoundingBox", _spy)
        with pytest.raises(AttributeError):
            sweep.calculate(_isopleths(), _demography(), mathematical_angle=0.0)

        assert seen["angle"] == 0.0
        assert seen["kwargs"]["geometryColumn"] == "TotalPolygon"

    def test_a_meteorological_angle_is_converted_before_the_sweep(self, sweep,
                                                                  monkeypatch):
        """Characterisation of B263, and the one piece of `calculate`'s own
        logic that is observable: the angle conversion."""
        from hera.utils import toMathematicalAngle

        seen = {}
        monkeypatch.setattr(type(sweep), "_findBoundingBox",
                            lambda self, *a, **kw: seen.setdefault(
                                "angle", kw["mathematical_angle"]) and [])
        with pytest.raises(Exception):
            sweep.calculate(_isopleths(), _demography(), meteorological_angle=90.0)

        assert seen["angle"] == toMathematicalAngle(90.0)
