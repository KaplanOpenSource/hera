"""thresholdGeoDataFrame: shifting/rotating threshold polygons, and projecting
them onto a demographic layer.

B65, resolved independently: an earlier pass found ``project()`` opening
with ``isinstance(meteorological_angle, collections.Iterable)``, an
attribute removed from the ``collections`` module in Python 3.10 (moved to
``collections.abc.Iterable`` in 3.3). That was a real, accurate finding at
the time; commit db405cec ("fix: repair Python 3.10+ and typo crashes in
risk assessment modules") fixed it to ``collections.abc.Iterable``,
unrelated to this test-expansion effort. ``project()``'s single/list angle
dispatch is exercised below against a monkeypatched ``_project`` (the real
one needs a full geopandas demographic layer with population columns and
a live CRS conversion -- out of scope for a hermetic unit test).

B97: ``shiftLocationAndAngle`` silently loses the ``thresholdGeoDataFrame``
subclass. It calls ``self.copy()``, and under the pinned geopandas
(1.0.1), ``thresholdGeoDataFrame`` -- a ``geopandas.GeoDataFrame``
subclass with no ``_constructor`` override -- has its type erased by
``.copy()`` down to a plain ``GeoDataFrame``. Chaining another
``thresholdGeoDataFrame``-only method (``shiftLocationAndAngle`` again,
or ``project``) onto the result fails with ``AttributeError``.

B108: ``_calculatePopulationInPolygon`` always raises. It computes
``demography.loc[not demography["geometry"].intersection(poly).is_empty]``
-- ``not`` on a multi-element pandas Series raises
``ValueError: The truth value of a Series is ambiguous``, and pandas
raises that same error even for a length-1 Series (unlike a numpy
scalar). Every call fails, regardless of how many rows demography has;
the elementwise negation it clearly needs is ``~`` (or ``.apply(...)``),
not ``not``.
"""
import pandas
import pytest
from shapely.geometry import Point, Polygon

from hera.riskassessment.agents.effects.thresholdGeoDataFrame import thresholdGeoDataFrame


def _single_polygon_frame(poly=None):
    poly = poly or Polygon([(0, 0), (10, 0), (10, 5), (0, 5)])
    return thresholdGeoDataFrame(
        {"severity": ["Severe"], "ToxicLoad": [1.0], "ThresholdPolygon": [poly]},
        geometry="ThresholdPolygon",
    )


@pytest.mark.unit
class TestShiftLocationAndAngleLosesTheSubclass:
    """B97: see the module docstring."""

    @pytest.mark.xfail(
        strict=True,
        reason="B97: shiftLocationAndAngle calls self.copy(), which under "
               "the pinned geopandas (1.0.1) erases the thresholdGeoDataFrame "
               "subclass (no _constructor override) down to a plain "
               "GeoDataFrame. See the consolidated findings issue.",
    )
    def test_the_result_should_still_be_a_threshold_geo_data_frame(self):
        shifted = _single_polygon_frame().shiftLocationAndAngle(loc=(0, 0), meteorological_angle=0)
        assert isinstance(shifted, thresholdGeoDataFrame)

    def test_the_result_is_currently_a_plain_geodataframe(self):
        """Characterisation of B97."""
        import geopandas

        shifted = _single_polygon_frame().shiftLocationAndAngle(loc=(0, 0), meteorological_angle=0)
        assert type(shifted) is geopandas.GeoDataFrame

    def test_chaining_shift_again_on_the_result_currently_fails(self):
        """Characterisation of B97: losing the subclass breaks chaining."""
        shifted = _single_polygon_frame().shiftLocationAndAngle(loc=(0, 0), meteorological_angle=0)
        with pytest.raises(AttributeError, match="shiftLocationAndAngle"):
            shifted.shiftLocationAndAngle(loc=(1, 1), meteorological_angle=0)


@pytest.mark.unit
class TestShiftLocationAndAngle:
    def test_a_pure_translation_moves_every_vertex_by_the_offset(self):
        """meteorological_angle=270 is mathematical_angle=0 -- no rotation."""
        original = _single_polygon_frame()
        shifted = original.shiftLocationAndAngle(loc=(100, 200), meteorological_angle=270)
        assert shifted["ThresholdPolygon"].iloc[0].bounds == pytest.approx((100.0, 200.0, 110.0, 205.0))

    def test_either_angle_convention_is_accepted(self):
        original = _single_polygon_frame()
        by_met = original.shiftLocationAndAngle(loc=(0, 0), meteorological_angle=270)
        by_math = original.shiftLocationAndAngle(loc=(0, 0), mathematical_angle=0)
        assert by_met["ThresholdPolygon"].iloc[0].equals_exact(by_math["ThresholdPolygon"].iloc[0], tolerance=1e-9)

    def test_without_any_angle_it_raises(self):
        with pytest.raises(ValueError, match="met_angle or math_angle"):
            _single_polygon_frame().shiftLocationAndAngle(loc=(0, 0))

    def test_the_original_frame_is_left_unchanged(self):
        original = _single_polygon_frame()
        original_bounds = original["ThresholdPolygon"].iloc[0].bounds
        original.shiftLocationAndAngle(loc=(100, 200), meteorological_angle=90)
        assert original["ThresholdPolygon"].iloc[0].bounds == original_bounds


@pytest.mark.unit
class TestProjectAngleDispatch:
    """B65 is resolved (see the module docstring) -- these exercise the
    single-vs-list angle dispatch against a monkeypatched _project, since
    the real one needs a live geopandas demographic layer."""

    def test_a_single_meteorological_angle_calls_project_once(self, monkeypatch):
        calls = []

        def fake_project(self, **kwargs):
            calls.append(kwargs)
            return pandas.DataFrame({"casualties": [1]})

        monkeypatch.setattr(thresholdGeoDataFrame, "_project", fake_project)
        result = _single_polygon_frame().project(demographic="demog", loc=(0, 0), meteorological_angle=90)
        assert len(calls) == 1
        assert calls[0]["meteorological_angle"] == 90
        assert len(result) == 1

    def test_a_list_of_meteorological_angles_calls_project_once_per_angle(self, monkeypatch):
        calls = []

        def fake_project(self, **kwargs):
            calls.append(kwargs["meteorological_angle"])
            return pandas.DataFrame({"casualties": [1]})

        monkeypatch.setattr(thresholdGeoDataFrame, "_project", fake_project)
        result = _single_polygon_frame().project(demographic="demog", loc=(0, 0), meteorological_angle=[0, 90, 180])
        assert calls == [0, 90, 180]
        assert len(result) == 3

    def test_the_list_dispatch_records_both_angle_conventions_per_row(self, monkeypatch):
        monkeypatch.setattr(
            thresholdGeoDataFrame, "_project",
            lambda self, **kwargs: pandas.DataFrame({"casualties": [1]}),
        )
        result = _single_polygon_frame().project(demographic="demog", loc=(0, 0), meteorological_angle=[90])
        assert result["meteorological_angle"].iloc[0] == 90
        assert "mathematical_angle" in result.columns

    def test_a_list_of_mathematical_angles_calls_project_once_per_angle(self, monkeypatch):
        calls = []

        def fake_project(self, **kwargs):
            calls.append(kwargs["mathematical_angle"])
            return pandas.DataFrame({"casualties": [1]})

        monkeypatch.setattr(thresholdGeoDataFrame, "_project", fake_project)
        result = _single_polygon_frame().project(demographic="demog", loc=(0, 0), mathematical_angle=[0, 45])
        assert calls == [0, 45]
        assert len(result) == 2


@pytest.mark.unit
class TestCalculatePopulationInPolygonIsBroken:
    """B108: see the module docstring."""

    @staticmethod
    def _demography(n_rows=2):
        import geopandas

        return geopandas.GeoDataFrame({
            "geometry": [Point(i * 5, i * 5).buffer(1) for i in range(n_rows)],
            "total_pop": [100 * (i + 1) for i in range(n_rows)],
        })

    @pytest.mark.xfail(
        strict=True,
        reason="B108: `not <Series>` raises ValueError for any "
               "multi-element pandas Series -- the elementwise negation "
               "needed here is `~`, not `not`. See the consolidated "
               "findings issue.",
    )
    def test_it_should_return_the_intersecting_rows_with_population(self):
        poly = Point(0, 0).buffer(2)
        result = thresholdGeoDataFrame._calculatePopulationInPolygon(
            self._demography(), poly, ["total_pop"],
        )
        assert len(result) == 1

    def test_it_currently_raises_for_multiple_rows(self):
        """Characterisation of B108."""
        poly = Point(0, 0).buffer(2)
        with pytest.raises(ValueError, match="ambiguous"):
            thresholdGeoDataFrame._calculatePopulationInPolygon(
                self._demography(), poly, ["total_pop"],
            )

    def test_it_currently_raises_even_for_a_single_row(self):
        """Characterisation of B108: pandas' Series.__bool__ raises the
        same way regardless of length, unlike a numpy scalar."""
        poly = Point(0, 0).buffer(2)
        with pytest.raises(ValueError, match="ambiguous"):
            thresholdGeoDataFrame._calculatePopulationInPolygon(
                self._demography(n_rows=1), poly, ["total_pop"],
            )
