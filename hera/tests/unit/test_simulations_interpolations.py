"""Inverse-distance spatial interpolation between measuring stations.

The weights are 1/(1 + r^2/a^2) horizontally and 1/(1 + C' dz^2/a^2)
vertically, normalised to sum to one.  That makes the interpolant a convex
combination of the station values, which gives three properties worth pinning:

* at a station's own coordinates the station's value comes back exactly
* the result is bracketed by the smallest and largest station value
* identical stations give that common value regardless of geometry
"""
import pytest

from hera.simulations.utils.interpolations import spatialInterpolate


@pytest.fixture()
def interpolator():
    return spatialInterpolate()


def station(lat, lon, elevation, value):
    """A station as interp expects it: [lat, lon, elevation, value]."""
    return [lat, lon, elevation, value]


@pytest.mark.unit
class TestExactAtStations:
    def test_a_point_on_a_station_returns_that_station_value(self, interpolator):
        """The early-exit branch: an exact coordinate match short-circuits."""
        stations = [station(0, 0, 0, 10.0), station(100, 100, 0, 20.0)]
        assert interpolator.interp([0, 0, 0], stations) == 10.0

    def test_it_matches_the_second_station_too(self, interpolator):
        stations = [station(0, 0, 0, 10.0), station(100, 100, 0, 20.0)]
        assert interpolator.interp([100, 100, 0], stations) == 20.0

    def test_the_match_requires_all_three_coordinates(self, interpolator):
        """Same lat/lon but a different elevation is not an exact match."""
        stations = [station(0, 0, 0, 10.0), station(100, 100, 0, 20.0)]
        assert interpolator.interp([0, 0, 50], stations) != 10.0


@pytest.mark.unit
class TestConvexCombination:
    def test_the_result_is_bracketed_by_the_station_values(self, interpolator):
        stations = [station(0, 0, 0, 10.0), station(100, 100, 0, 20.0)]
        value = interpolator.interp([50, 50, 0], stations)
        assert 10.0 <= value <= 20.0

    def test_identical_stations_give_their_common_value(self, interpolator):
        """Normalised weights sum to one, so the answer cannot drift."""
        stations = [
            station(0, 0, 0, 7.0),
            station(100, 0, 0, 7.0),
            station(0, 100, 0, 7.0),
        ]
        assert interpolator.interp([33, 44, 0], stations) == pytest.approx(7.0)

    def test_the_midpoint_of_two_equal_stations_is_their_mean(self, interpolator):
        """Symmetric geometry means equal weights."""
        stations = [station(-100, 0, 0, 10.0), station(100, 0, 0, 20.0)]
        assert interpolator.interp([0, 0, 0], stations) == pytest.approx(15.0)

    def test_a_single_station_dominates_everywhere(self, interpolator):
        stations = [station(0, 0, 0, 42.0)]
        assert interpolator.interp([500, 500, 0], stations) == pytest.approx(42.0)


@pytest.mark.unit
class TestDistanceWeighting:
    def test_the_nearer_station_pulls_harder(self, interpolator):
        stations = [station(0, 0, 0, 10.0), station(1000, 0, 0, 20.0)]
        near = interpolator.interp([10, 0, 0], stations)
        assert near < 15.0, "the far station should not dominate a near point"

    def test_moving_towards_a_station_moves_the_value_towards_it(self, interpolator):
        stations = [station(0, 0, 0, 10.0), station(1000, 0, 0, 20.0)]
        values = [
            interpolator.interp([offset, 0, 0], stations)
            for offset in (100, 300, 500, 700, 900)
        ]
        assert values == sorted(values), "the interpolant is not monotonic along the line"

    def test_a_larger_cell_widens_the_influence(self, interpolator):
        """a^2 = (dx^2 + dy^2)/4, so a bigger cell flattens the weighting."""
        stations = [station(0, 0, 0, 10.0), station(1000, 0, 0, 20.0)]
        tight = interpolator.interp([100, 0, 0], stations, dx=20, dy=20)
        broad = interpolator.interp([100, 0, 0], stations, dx=2000, dy=2000)
        assert broad > tight


@pytest.mark.unit
class TestVerticalWeighting:
    def test_an_elevation_difference_reduces_a_station_weight(self, interpolator):
        """Two stations at the same place, one far above: the near one wins."""
        stations = [station(0, 0, 0, 10.0), station(0, 0, 1000, 20.0)]
        assert interpolator.interp([0, 0, 1, ], stations) < 15.0

    def test_topography_switches_to_the_near_ground_coefficient(self, interpolator):
        """With topography given, C' follows a tanh ramp from D to C."""
        stations = [station(0, 0, 0, 10.0), station(0, 0, 500, 20.0)]
        without = interpolator.interp([0, 0, 100], stations)
        with_ground = interpolator.interp([0, 0, 100], stations, topography=0)
        assert without != with_ground

    def test_the_result_stays_bracketed_with_topography(self, interpolator):
        stations = [station(0, 0, 0, 10.0), station(0, 0, 500, 20.0)]
        value = interpolator.interp([0, 0, 100], stations, topography=0)
        assert 10.0 <= value <= 20.0


@pytest.mark.unit
class TestVectorValues:
    def test_a_vector_station_returns_a_vector(self, interpolator):
        stations = [
            station(0, 0, 0, [1.0, 2.0]),
            station(100, 0, 0, [3.0, 4.0]),
        ]
        result = interpolator.interp([50, 0, 0], stations)
        assert isinstance(result, list)
        assert len(result) == 2

    def test_each_component_is_interpolated_independently(self, interpolator):
        stations = [
            station(-100, 0, 0, [10.0, 100.0]),
            station(100, 0, 0, [20.0, 200.0]),
        ]
        first, second = interpolator.interp([0, 0, 0], stations)
        assert first == pytest.approx(15.0)
        assert second == pytest.approx(150.0)

    def test_a_vector_station_hit_exactly_returns_its_own_vector(self, interpolator):
        stations = [
            station(0, 0, 0, [1.0, 2.0]),
            station(100, 0, 0, [3.0, 4.0]),
        ]
        assert interpolator.interp([0, 0, 0], stations) == [1.0, 2.0]

    def test_components_stay_bracketed(self, interpolator):
        stations = [
            station(0, 0, 0, [1.0, 50.0]),
            station(100, 100, 0, [3.0, 60.0]),
        ]
        first, second = interpolator.interp([50, 50, 0], stations)
        assert 1.0 <= first <= 3.0
        assert 50.0 <= second <= 60.0
