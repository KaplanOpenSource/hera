"""simulations/utils/interpolations.py: spatialInterpolate, a
distance-weighted station interpolator, and two bugs in it.

B93 (fixed): ``windprofile`` was defined inside the ``spatialInterpolate``
class body but its signature had no ``self``, so calling it on an instance
bound the instance itself to ``z``. It is a ``staticmethod`` now, which is
what the signature always meant, and both calling forms work.

B94, retracted: an earlier pass found ``interpPandas`` writing each
interpolated value with ``points["interpulation"][i] = self.interp(...)``
-- a chained assignment -- and, verified against a locally-drifted pandas
(3.0.2) where copy-on-write is mandatory, that assignment silently never
updated the DataFrame. Under the pandas actually pinned in
requirements.txt (2.2.3), copy-on-write is opt-in, not yet the default:
the chained assignment still works, with only a FutureWarning/
SettingWithCopyWarning. ``interpPandas`` does fill in the "interpulation"
column correctly today. This bug would resurface if/when pandas is
upgraded past 2.x with copy-on-write made the default -- worth knowing,
not worth pinning as a currently-live defect.

B95 (fixed): ``interpArray`` called ``self.interpPandas(points=...,
stations=..., topography=newPoints["topography"], ...)``, but
``interpPandas`` has no ``topography`` parameter -- only a ``columnNames``
dict that already maps ``"topography"`` -- so every call raised
``TypeError``.
"""
import numpy
import pandas
import pytest

from hera.simulations.utils.interpolations import spatialInterpolate


@pytest.fixture()
def si():
    return spatialInterpolate()


@pytest.mark.unit
class TestInterp:
    def test_an_exact_station_match_returns_its_value_directly(self, si):
        stations = [[0, 0, 0, 10.0], [10, 0, 0, 20.0]]
        assert si.interp([0, 0, 0], stations) == 10.0

    def test_a_midpoint_is_weighted_between_two_equidistant_stations(self, si):
        stations = [[0, 0, 0, 10.0], [10, 0, 0, 20.0]]
        result = si.interp([5, 0, 0], stations)
        assert result == pytest.approx(15.0)

    def test_a_point_closer_to_one_station_leans_toward_its_value(self, si):
        stations = [[0, 0, 0, 10.0], [10, 0, 0, 20.0]]
        near_first = si.interp([2, 0, 0], stations)
        near_second = si.interp([8, 0, 0], stations)
        assert near_first < near_second

    def test_vector_valued_stations_interpolate_elementwise(self, si):
        stations = [[0, 0, 0, [10.0, 100.0]], [10, 0, 0, [20.0, 200.0]]]
        result = si.interp([5, 0, 0], stations)
        assert result == pytest.approx([15.0, 150.0])

    def test_a_topography_offset_still_produces_a_finite_result(self, si):
        stations = [[0, 0, 0, 10.0], [10, 0, 0, 20.0]]
        result = si.interp([5, 0, 5], stations, topography=2)
        assert numpy.isfinite(result)


@pytest.mark.unit
class TestCheckInterpulation:
    def test_it_returns_one_percentage_difference_list_per_station(self, si):
        stations = [[0, 0, 0, [10.0]], [10, 0, 0, [20.0]], [0, 10, 0, [15.0]]]
        differences = si.checkInterpulation(stations)
        assert len(differences) == len(stations)
        assert all(len(d) == 1 for d in differences)

    def test_a_station_matching_its_neighbors_perfectly_has_zero_difference(self, si):
        stations = [[0, 0, 0, [10.0]], [10, 0, 0, [10.0]], [5, 0, 0, [10.0]]]
        differences = si.checkInterpulation(stations)
        assert differences[2][0] == pytest.approx(0.0, abs=1e-6)


@pytest.mark.unit
class TestInterpPandas:
    """B94 is retracted (see the module docstring) -- under the pinned
    pandas (2.2.3), the chained assignment inside interpPandas still
    works."""

    def test_it_fills_in_the_interpolated_value(self, si):
        stations = [[0, 0, 0, 10.0], [10, 0, 0, 20.0]]
        points = pandas.DataFrame({"x": [5], "y": [0], "z": [0]})
        result = si.interpPandas(points, stations)
        assert result["interpulation"].iloc[0] == pytest.approx(15.0)

@pytest.mark.unit
class TestInterpArray:
    """B95: interpArray used to pass topography=... to interpPandas, which
    has no such parameter -- the topography column reaches interp through
    the columnNames mapping instead."""

    def test_it_returns_an_interpolated_column(self, si):
        stations = [[0, 0, 0, 10.0], [10, 0, 0, 20.0]]
        points = pandas.DataFrame({"x": [5], "y": [0], "z": [0], "topography": [0]})
        result = si.interpArray(points, stations)
        assert result["interpulation"].iloc[0] == pytest.approx(15.0)


@pytest.mark.unit
class TestWindprofile:
    def test_calling_it_on_an_instance_computes_a_velocity(self, si):
        assert si.windprofile(24) == pytest.approx(3.0, abs=1e-6)

    def test_it_also_works_when_called_on_the_class_directly(self):
        """It is a staticmethod, so both calling forms reach the same code."""
        result = spatialInterpolate.windprofile(24)
        assert result == pytest.approx(3.0, abs=1e-6)
