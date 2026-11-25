"""InjuryLevel.getContoursFromLevels: matplotlib contour extraction from a
2D toxic-load field."""
import numpy
import pytest
import xarray

from hera.riskassessment.agents.effects.InjuryLevel import InjuryLevel


@pytest.fixture()
def il():
    return InjuryLevel.__new__(InjuryLevel)


def _field(dims):
    x = numpy.linspace(0, 10, 5)
    y = numpy.linspace(0, 10, 5)
    data = numpy.random.RandomState(0).rand(5, 5)
    return xarray.DataArray(data, dims=dims, coords={"x": x, "y": y})


@pytest.mark.unit
class TestGetContoursFromLevels:
    def test_it_returns_a_contour_set_for_yx_dims(self, il):
        result = il.getContoursFromLevels(_field(["y", "x"]), "x", "y", 0.5)
        import matplotlib.contour

        assert isinstance(result, matplotlib.contour.QuadContourSet)

    def test_it_transposes_xy_dims_to_match(self, il):
        import matplotlib.contour

        result = il.getContoursFromLevels(_field(["x", "y"]), "x", "y", 0.5)
        assert isinstance(result, matplotlib.contour.QuadContourSet)

    def test_a_third_dimension_that_does_not_squeeze_away_raises(self, il):
        field = xarray.DataArray(
            numpy.zeros((5, 5, 2)), dims=["y", "x", "z"],
            coords={"x": range(5), "y": range(5), "z": range(2)},
        )
        with pytest.raises(Exception, match="Cannot extract contours"):
            il.getContoursFromLevels(field, "x", "y", 0.5)
