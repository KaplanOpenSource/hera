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

Two bugs were found here, both now fixed.  B263: ``_doCalculation`` called
a method that does not exist (``datalayer`` instead of ``project``), so no
sweep could complete.  B266: ``_findBoundingBox`` could not even size its
mesh when the final time step had a single isopleth row, which is the normal
shape for a single-severity injury.  A real (unmocked) projection is still
blocked one line past B263 by a geometry column reaching ``groupby().sum()``
-- a separate, later defect, documented with the two xfails below.
``calculate``'s parallel branch is not exercised here because
``multiprocessing.Pool`` inside a hermetic unit test would fork the
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
class TestASingleIsoplethRow:
    """B266, fixed: see the reasons below."""

    @staticmethod
    def _single():
        """One severity, one time step -- an injury with a single level."""
        frame = _isopleths()
        frame = frame[frame["severity"] == "Light"]
        return thresholdGeoDataFrame(frame, geometry="ThresholdPolygon")

    def test_a_single_severity_still_sizes_the_mesh(self, sweep):
        """B266, fixed: _findBoundingBox narrowed the isopleths to the last
        time step with .loc[maxdatetime].  When only one row carried that
        timestamp -- the normal shape for an injury with a single severity
        level -- pandas returned a Series rather than a one-row frame, and
        the next line raised AttributeError: 'Series' object has no
        attribute 'set_geometry'.  .loc[[maxdatetime]] keeps the frame."""
        mesh = sweep._findBoundingBox(self._single(), _demography(),
                                      mathematical_angle=0.0)

        expected = (len(numpy.arange(EXPECTED_MIN_X, EXPECTED_MAX_X, DXDY))
                    * len(numpy.arange(EXPECTED_MIN_Y, EXPECTED_MAX_Y, DXDY)))
        assert len(mesh) == expected

    def test_a_second_row_at_the_same_time_step_gives_the_same_mesh(self, sweep):
        """The row count used to be the trigger.  The extra row is a
        duplicate of the first, so it must produce the same mesh."""
        doubled = self._single()
        doubled = thresholdGeoDataFrame(
            pandas.concat([doubled, doubled], ignore_index=True),
            geometry="ThresholdPolygon")

        mesh = sweep._findBoundingBox(doubled, _demography(), mathematical_angle=0.0)

        expected = (len(numpy.arange(EXPECTED_MIN_X, EXPECTED_MAX_X, DXDY))
                    * len(numpy.arange(EXPECTED_MIN_Y, EXPECTED_MAX_Y, DXDY)))
        assert len(mesh) == expected

    def test_the_scalar_loc_selection_is_what_used_to_degrade(self):
        """B266's mechanism, derived from pandas rather than from hera."""
        frame = self._single().set_index("datetime")
        stamp = frame.index.max()

        assert isinstance(frame.loc[stamp], pandas.Series)
        assert isinstance(frame.loc[[stamp]], pandas.DataFrame)


@pytest.mark.unit
class TestSweepCannotRun:
    """B263, fixed; a second, later defect still blocks the sweep -- see
    the reasons below."""

    @staticmethod
    def _params():
        return {"effectIsopleths": _isopleths(),
                "demog": _demography(),
                "rotate_angle": 0.0,
                "valueColumn": "effectedtotal_pop"}

    @staticmethod
    def _noCasualties(monkeypatch, calls=None):
        """Replace the projection itself, so what is exercised is the sweep
        driver and its call site rather than the geopandas projection."""
        def _project(self, demographic, loc, **kwargs):
            if calls is not None:
                calls.append((demographic, loc, kwargs))
            return None

        monkeypatch.setattr(thresholdGeoDataFrame, "project", _project)

    REASON = (
        "Not B263 (that one -- _doCalculation calling the non-existent "
        "effectIsopleths.datalayer(...) instead of .project(...) -- is "
        "fixed, and the call site is asserted below).  What now fails is "
        "the next line: projected.groupby([...]).sum() sums a frame that "
        "still carries its geometry column, which pandas 2.2 / geopandas "
        "1.1 reject with \"TypeError: agg function failed [how->sum, "
        "dtype->geometry]\" -- older pandas dropped non-numeric columns "
        "silently.  Only the numeric columns are read out of that sum, so "
        "the intent is a numeric-only aggregation, but that is a separate "
        "defect from B263 and is not fixed here."
    )

    @pytest.mark.xfail(strict=True, reason=REASON)
    def test_one_release_point_should_yield_casualties_per_severity(self, sweep):
        result = sweep._doCalculation((0.0, 0.0), self._params())[0]

        assert sorted(result["severity"].unique()) == ["Light", "Severe"]

    def test_the_release_point_is_projected_onto_the_demography(self, sweep, monkeypatch):
        """B263, fixed: _doCalculation used to call the non-existent
        effectIsopleths.datalayer(...).  It now calls project(), with the
        demography, the release location and the rotated angle."""
        calls = []
        self._noCasualties(monkeypatch, calls)
        params = self._params()

        sweep._doCalculation((1.0, 2.0), params)

        demographic, loc, kwargs = calls[0]
        assert demographic is params["demog"]
        assert loc == (1.0, 2.0)
        assert kwargs["mathematical_angle"] == 0.0

    def test_a_release_point_that_affects_nobody_reports_zero_per_severity(self, sweep, monkeypatch):
        self._noCasualties(monkeypatch)

        data = sweep._doCalculation((1.0, 2.0), self._params())[0]

        assert sorted(data["severity"]) == ["Light", "Severe"]
        assert data["effectedPopulation"].tolist() == [0, 0]
        assert data["x"].tolist() == [1.0, 1.0]

    @pytest.mark.xfail(strict=True, reason=REASON)
    def test_a_serial_sweep_should_return_one_row_per_point_and_severity(self, sweep):
        result = sweep.calculate(_isopleths(), _demography(), mathematical_angle=0.0)

        assert set(result.columns) >= {"x", "y", "severity", "datetime"}

    def test_a_serial_sweep_walks_every_mesh_point(self, sweep, monkeypatch, capsys):
        """The driver itself is fine: with the projection replaced, it walks
        the whole 6 x 7 mesh and reports both severities at every point."""
        self._noCasualties(monkeypatch)

        result = sweep.calculate(_isopleths(), _demography(), mathematical_angle=0.0)

        assert len(result) == 42 * 2
        assert set(result.columns) >= {"x", "y", "severity", "datetime"}
        assert "Processing point 0 out of 42" in capsys.readouterr().out

    def test_the_mesh_is_built_with_the_converted_angle(self, sweep, monkeypatch):
        seen = {}
        original = type(sweep)._findBoundingBox

        def _spy(self, effectIsopleths, demog, mathematical_angle, **kwargs):
            seen["angle"] = mathematical_angle
            seen["kwargs"] = kwargs
            return original(self, effectIsopleths, demog, mathematical_angle, **kwargs)

        monkeypatch.setattr(type(sweep), "_findBoundingBox", _spy)
        self._noCasualties(monkeypatch)
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
