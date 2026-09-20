"""gaussian/Sigma.py: ``AbstractSigma.getVirtualDistance`` as a base-class
algorithm in its own right.

test_gaussian_sigma.py covers it only through ``BriggsRural``, where the
Briggs coefficient table and the virtual-source solution are entangled.
The method itself is generic: given an initial cloud size it inverts the
subclass's own sigma curves by root finding, once per axis.  The tests
below pin that contract against a deliberately trivial sigma model whose
inverse is known in closed form, so the base class is checked without
Briggs in the way, and confirm that the unimplemented base ``getSigma``
propagates out of the solver rather than being swallowed.

Not covered here: ``BriggsRural.getSigma`` / ``__call__`` / the
coefficient tables and the Briggs virtual distances -- all in
test_gaussian_sigma.py.  No new defects surfaced.
"""
import numpy
import pytest

from hera.simulations.gaussian.Sigma import AbstractSigma, BriggsRural
from hera.utils.unitHandler import ureg


class LinearSigma(AbstractSigma):
    """sigma_i = slope_i * x, so the virtual distance is sigma0_i / slope_i."""

    SLOPES = {"sigmaX": 0.1, "sigmaY": 0.2, "sigmaZ": 0.5}

    def __init__(self):
        self.calls = []

    def getSigma(self, x, stability, sigma0=None, units=True):
        self.calls.append((x, stability, sigma0, units))
        x = numpy.atleast_1d(numpy.asarray(x, dtype=float))
        result = {name: slope * x for name, slope in self.SLOPES.items()}
        result["distance"] = x
        return result


@pytest.mark.unit
class TestVirtualDistanceOnTheAbstractBase:
    def test_the_unimplemented_base_sigma_propagates_out_of_the_solver(self):
        """The root finder must not hide a missing implementation."""
        with pytest.raises(NotImplementedError, match="not implemented for abstractSigma"):
            AbstractSigma().getVirtualDistance((10 * ureg.m,) * 3, "D")

    def test_it_fails_for_every_stability_class(self):
        for stability in ("A", "B", "C", "D", "E", "F"):
            with pytest.raises(NotImplementedError):
                AbstractSigma().getVirtualDistance((1 * ureg.m,) * 3, stability)


@pytest.mark.unit
class TestVirtualDistanceInvertsTheSubclassCurves:
    @pytest.fixture()
    def model(self):
        return LinearSigma()

    def test_each_axis_is_solved_from_its_own_curve(self, model):
        """sigma0 = (1, 2, 5) m against slopes (0.1, 0.2, 0.5) all invert to 10 m."""
        distances = model.getVirtualDistance(
            (1 * ureg.m, 2 * ureg.m, 5 * ureg.m), "D"
        )
        for distance in distances:
            assert distance.m_as(ureg.m) == pytest.approx(10.0, rel=1e-6)

    def test_the_solution_reproduces_the_initial_size(self, model):
        """The defining property, checked against the closed-form inverse."""
        sigma0 = (3 * ureg.m, 4 * ureg.m, 5 * ureg.m)
        alongwind, crosswind, vertical = model.getVirtualDistance(sigma0, "D")
        assert alongwind.m_as(ureg.m) == pytest.approx(3.0 / 0.1, rel=1e-6)
        assert crosswind.m_as(ureg.m) == pytest.approx(4.0 / 0.2, rel=1e-6)
        assert vertical.m_as(ureg.m) == pytest.approx(5.0 / 0.5, rel=1e-6)

    def test_it_returns_three_lengths_in_metres(self, model):
        distances = model.getVirtualDistance((1 * ureg.m,) * 3, "D")
        assert len(distances) == 3
        for distance in distances:
            assert distance.units == ureg.m

    def test_the_solution_is_linear_in_the_initial_size(self, model):
        single = model.getVirtualDistance((1 * ureg.m,) * 3, "D")[0].m_as(ureg.m)
        tenfold = model.getVirtualDistance((10 * ureg.m,) * 3, "D")[0].m_as(ureg.m)
        assert tenfold == pytest.approx(10 * single, rel=1e-6)

    def test_a_zero_initial_size_gives_a_point_source(self, model):
        distances = model.getVirtualDistance((0 * ureg.m,) * 3, "D")
        for distance in distances:
            assert distance.m_as(ureg.m) == pytest.approx(0.0, abs=1e-9)

    def test_the_initial_size_may_be_given_in_any_length_unit(self, model):
        """Each component goes through tonumber(..., m), so units convert."""
        mixed = model.getVirtualDistance(
            (0.01 * ureg.km, 1000 * ureg.cm, 10 * ureg.m), "D"
        )
        plain = model.getVirtualDistance((10 * ureg.m,) * 3, "D")
        for fromMixed, fromPlain in zip(mixed, plain):
            assert fromMixed.m_as(ureg.m) == pytest.approx(
                fromPlain.m_as(ureg.m), rel=1e-6
            )

    def test_the_stability_class_is_handed_to_the_subclass_unchanged(self, model):
        model.getVirtualDistance((1 * ureg.m,) * 3, "E")
        assert {call[1] for call in model.calls} == {"E"}

    def test_the_subclass_is_asked_for_a_point_source_curve(self, model):
        """sigma0 must not be fed back in, or the inversion would recurse."""
        model.getVirtualDistance((1 * ureg.m,) * 3, "D")
        assert all(call[2] is None for call in model.calls)

    def test_the_subclass_is_asked_for_unitless_values(self, model):
        model.getVirtualDistance((1 * ureg.m,) * 3, "D")
        assert all(call[3] is False for call in model.calls)

    def test_all_three_axes_are_solved(self, model):
        """One newton run per axis, so the subclass is called for each."""
        model.getVirtualDistance((1 * ureg.m, 2 * ureg.m, 3 * ureg.m), "D")
        assert len(model.calls) >= 3


@pytest.mark.unit
class TestVirtualDistanceAgreesWithTheBriggsSubclass:
    def test_the_briggs_solution_reproduces_its_own_initial_size(self):
        """Cross-check: the same base-class algorithm, a real sigma model."""
        sigma = BriggsRural()
        alongwind, crosswind, vertical = sigma.getVirtualDistance(
            (7 * ureg.m, 7 * ureg.m, 3 * ureg.m), "D"
        )
        atVirtualSource = sigma.getSigma(
            x=[
                alongwind.m_as(ureg.m),
                crosswind.m_as(ureg.m),
                vertical.m_as(ureg.m),
            ],
            stability="D",
            units=False,
        )
        assert atVirtualSource["sigmaX"][0] == pytest.approx(7.0, rel=1e-6)
        assert atVirtualSource["sigmaY"][1] == pytest.approx(7.0, rel=1e-6)
        assert atVirtualSource["sigmaZ"][2] == pytest.approx(3.0, rel=1e-6)
