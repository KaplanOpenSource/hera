"""thresholdGeoDataFrame: shifting/rotating threshold polygons, and projecting
them onto a demographic layer.

``project`` cannot run at all on this Python:

* B65: it starts with ``isinstance(meteorological_angle,
  collections.Iterable)``, but ``collections.Iterable`` was removed in
  Python 3.10 (moved to ``collections.abc.Iterable`` in 3.3, alias dropped
  in 3.10). Every call raises ``AttributeError`` before touching any of its
  arguments -- including the demographic data, so a single-angle call
  fails exactly like a list-of-angles call.
"""
import pytest
from shapely.geometry import Polygon

from hera.riskassessment.agents.effects.thresholdGeoDataFrame import thresholdGeoDataFrame


def _single_polygon_frame(poly=None):
    poly = poly or Polygon([(0, 0), (10, 0), (10, 5), (0, 5)])
    return thresholdGeoDataFrame(
        {"severity": ["Severe"], "ToxicLoad": [1.0], "ThresholdPolygon": [poly]},
        geometry="ThresholdPolygon",
    )


@pytest.mark.unit
class TestShiftLocationAndAngle:
    def test_the_result_is_still_a_threshold_geo_data_frame(self):
        shifted = _single_polygon_frame().shiftLocationAndAngle(loc=(0, 0), meteorological_angle=0)
        assert isinstance(shifted, thresholdGeoDataFrame)

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
class TestProjectIsUnusableOnThisPython:
    @pytest.mark.xfail(
        strict=True,
        reason="B65: project() opens with `isinstance(meteorological_angle, "
               "collections.Iterable)`, an attribute removed from the "
               "collections module in Python 3.10 (it lives at "
               "collections.abc.Iterable). Every call raises AttributeError "
               "immediately, regardless of arguments. "
               "See the consolidated findings issue.",
    )
    def test_project_with_a_single_angle_should_return_casualty_estimates(self):
        _single_polygon_frame().project(demographic=None, loc=(0, 0), meteorological_angle=90)

    def test_project_currently_raises_for_any_input_at_all(self):
        """Characterisation of B65: even a deliberately-broken demographic
        argument (None) never gets used -- the crash happens first."""
        with pytest.raises(AttributeError, match="Iterable"):
            _single_polygon_frame().project(demographic=None, loc=(0, 0), meteorological_angle=90)

    def test_project_raises_the_same_way_with_a_list_of_angles(self):
        with pytest.raises(AttributeError, match="Iterable"):
            _single_polygon_frame().project(demographic=None, loc=(0, 0), meteorological_angle=[0, 90])
