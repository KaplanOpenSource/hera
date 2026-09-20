r"""Briggs dispersion coefficients, checked against the published formulas.

Briggs (1973) open-country ("rural") sigmas, as they appear in the standard
references:

    stability   sigma_y                       sigma_z
    A           0.22 x (1 + 1e-4 x)^-1/2      0.20 x
    B           0.16 x (1 + 1e-4 x)^-1/2      0.12 x
    C           0.11 x (1 + 1e-4 x)^-1/2      0.08 x (1 + 2e-4 x)^-1/2
    D           0.08 x (1 + 1e-4 x)^-1/2      0.06 x (1 + 1.5e-3 x)^-1/2
    E           0.06 x (1 + 1e-4 x)^-1/2      0.03 x (1 + 3e-4 x)^-1
    F           0.04 x (1 + 1e-4 x)^-1/2      0.016 x (1 + 3e-4 x)^-1

The tests below evaluate those expressions independently and compare.  That is
the point: an assertion copied from the implementation's own output would pass
even if the coefficient table were wrong.
"""
import numpy as np
import pytest

from hera.simulations.gaussian.Sigma import AbstractSigma, BriggsRural, briggsRural
from hera.utils.unitHandler import ureg

STABILITIES = ["A", "B", "C", "D", "E", "F"]
DISTANCES = [10.0, 100.0, 500.0, 1000.0, 5000.0, 10000.0]


def published_sigma_y(x, stability):
    coefficient = {"A": 0.22, "B": 0.16, "C": 0.11, "D": 0.08, "E": 0.06, "F": 0.04}
    return coefficient[stability] * x * (1 + 1e-4 * x) ** -0.5


def published_sigma_z(x, stability):
    if stability == "A":
        return 0.20 * x
    if stability == "B":
        return 0.12 * x
    if stability == "C":
        return 0.08 * x * (1 + 2e-4 * x) ** -0.5
    if stability == "D":
        return 0.06 * x * (1 + 1.5e-3 * x) ** -0.5
    if stability == "E":
        return 0.03 * x * (1 + 3e-4 * x) ** -1
    return 0.016 * x * (1 + 3e-4 * x) ** -1


@pytest.fixture()
def sigma():
    return BriggsRural()


@pytest.mark.unit
class TestAgainstPublishedFormulas:
    @pytest.mark.parametrize("stability", STABILITIES)
    @pytest.mark.parametrize("distance", DISTANCES)
    def test_crosswind_sigma(self, sigma, stability, distance):
        computed = sigma.getSigma(distance, stability, units=False)["sigmaY"][0]
        assert computed == pytest.approx(published_sigma_y(distance, stability), rel=1e-9)

    @pytest.mark.parametrize("stability", STABILITIES)
    @pytest.mark.parametrize("distance", DISTANCES)
    def test_vertical_sigma(self, sigma, stability, distance):
        computed = sigma.getSigma(distance, stability, units=False)["sigmaZ"][0]
        assert computed == pytest.approx(published_sigma_z(distance, stability), rel=1e-9)

    def test_the_returned_distance_is_the_input(self, sigma):
        assert sigma.getSigma(1234.0, "D", units=False)["distance"][0] == pytest.approx(
            1234.0
        )


@pytest.mark.unit
class TestCoefficientTables:
    """The tables themselves, so a typo is caught at its source."""

    def test_crosswind_a_coefficients(self, sigma):
        assert sigma._coeffX["A"].tolist() == [0.22, 0.16, 0.11, 0.08, 0.06, 0.04]

    def test_crosswind_b_and_c_are_uniform(self, sigma):
        assert sigma._coeffX["B"].tolist() == [1e-4] * 6
        assert sigma._coeffX["C"].tolist() == [-0.5] * 6

    def test_vertical_a_coefficients(self, sigma):
        assert sigma._coeffZ["A"].tolist() == [0.2, 0.12, 0.08, 0.06, 0.03, 0.016]

    def test_vertical_b_coefficients(self, sigma):
        assert sigma._coeffZ["B"].tolist() == [0, 0, 2e-4, 1.5e-3, 3e-4, 3e-4]

    def test_vertical_c_coefficients(self, sigma):
        assert sigma._coeffZ["C"].tolist() == [1, 1, -0.5, -0.5, -1, -1]

    def test_tables_are_indexed_by_stability_class(self, sigma):
        assert sigma._coeffX.index.tolist() == STABILITIES
        assert sigma._coeffZ.index.tolist() == STABILITIES

    def test_unstable_classes_have_no_vertical_correction(self, sigma):
        """A and B are pure linear growth: sigma_z = a x, no (1+bx) factor."""
        for stability in ("A", "B"):
            assert sigma._coeffZ.loc[stability]["B"] == 0


@pytest.mark.unit
class TestPhysicalProperties:
    """Properties any dispersion model must have, whatever the coefficients."""

    @pytest.mark.parametrize("stability", STABILITIES)
    def test_a_point_source_has_no_initial_spread(self, sigma, stability):
        result = sigma.getSigma(0.0, stability, units=False)
        assert result["sigmaY"][0] == pytest.approx(0.0)
        assert result["sigmaZ"][0] == pytest.approx(0.0)

    @pytest.mark.parametrize("stability", STABILITIES)
    def test_spread_grows_with_distance(self, sigma, stability):
        values = [
            sigma.getSigma(x, stability, units=False)["sigmaY"][0] for x in DISTANCES
        ]
        assert values == sorted(values), "sigma_y is not monotonic in x"

    @pytest.mark.parametrize("stability", STABILITIES)
    def test_vertical_spread_grows_with_distance(self, sigma, stability):
        values = [
            sigma.getSigma(x, stability, units=False)["sigmaZ"][0] for x in DISTANCES
        ]
        assert values == sorted(values), "sigma_z is not monotonic in x"

    @pytest.mark.parametrize("distance", [100.0, 1000.0, 10000.0])
    def test_more_unstable_air_disperses_more(self, sigma, distance):
        """A is the most unstable class and F the most stable, so sigma must
        decrease monotonically from A to F at any fixed distance."""
        values = [
            sigma.getSigma(distance, stability, units=False)["sigmaY"][0]
            for stability in STABILITIES
        ]
        assert values == sorted(values, reverse=True)

    @pytest.mark.parametrize("distance", [100.0, 1000.0, 10000.0])
    def test_vertical_spread_also_orders_by_stability(self, sigma, distance):
        values = [
            sigma.getSigma(distance, stability, units=False)["sigmaZ"][0]
            for stability in STABILITIES
        ]
        assert values == sorted(values, reverse=True)

    @pytest.mark.parametrize("stability", STABILITIES)
    def test_all_sigmas_are_positive_downwind(self, sigma, stability):
        for x in DISTANCES:
            result = sigma.getSigma(x, stability, units=False)
            assert result["sigmaY"][0] > 0
            assert result["sigmaZ"][0] > 0

    def test_alongwind_equals_crosswind(self, sigma):
        """Both are computed from the crosswind table, so sigma_x == sigma_y.

        Pinned because it is a modelling choice, not an accident: _coeffX holds
        the published sigma_y coefficients and is used for both.
        """
        result = sigma.getSigma(1000.0, "D", units=False)
        assert result["sigmaX"][0] == pytest.approx(result["sigmaY"][0])


@pytest.mark.unit
class TestVirtualSource:
    """sigma0 is modelled by shifting the origin upwind to a virtual source."""

    def test_the_virtual_distance_reproduces_the_initial_size(self, sigma):
        """The defining property: at x=0 the spread must equal sigma0 exactly."""
        initial = (10 * ureg.m,) * 3
        result = sigma.getSigma(0.0, "D", sigma0=initial, units=False)
        assert result["sigmaY"][0] == pytest.approx(10.0, rel=1e-6)

    @pytest.mark.parametrize("size", [1.0, 10.0, 50.0])
    def test_it_holds_for_several_initial_sizes(self, sigma, size):
        initial = (size * ureg.m,) * 3
        result = sigma.getSigma(0.0, "D", sigma0=initial, units=False)
        assert result["sigmaY"][0] == pytest.approx(size, rel=1e-6)
        assert result["sigmaZ"][0] == pytest.approx(size, rel=1e-6)

    def test_an_initial_size_increases_the_spread_downwind(self, sigma):
        initial = (10 * ureg.m,) * 3
        with_cloud = sigma.getSigma(100.0, "D", sigma0=initial, units=False)["sigmaY"][0]
        point = sigma.getSigma(100.0, "D", units=False)["sigmaY"][0]
        assert with_cloud > point

    def test_the_virtual_distance_is_positive(self, sigma):
        distances = sigma.getVirtualDistance((10 * ureg.m,) * 3, "D")
        for value in distances:
            assert value.to(ureg.m).magnitude > 0

    def test_a_bigger_cloud_sits_further_upwind(self, sigma):
        small = sigma.getVirtualDistance((5 * ureg.m,) * 3, "D")[0]
        large = sigma.getVirtualDistance((50 * ureg.m,) * 3, "D")[0]
        assert large.to(ureg.m).magnitude > small.to(ureg.m).magnitude

    def test_alongwind_and_crosswind_virtual_distances_coincide(self, sigma):
        """Both are solved from the crosswind curve, so they must be equal."""
        alongwind, crosswind, _ = sigma.getVirtualDistance((10 * ureg.m,) * 3, "D")
        assert alongwind.to(ureg.m).magnitude == pytest.approx(
            crosswind.to(ureg.m).magnitude
        )


@pytest.mark.unit
class TestInputHandling:
    def test_a_scalar_becomes_a_one_element_result(self, sigma):
        assert len(sigma.getSigma(100.0, "D", units=False)["sigmaY"]) == 1

    def test_a_list_is_evaluated_elementwise(self, sigma):
        result = sigma.getSigma([100.0, 200.0, 300.0], "D", units=False)
        assert len(result["sigmaY"]) == 3
        for index, x in enumerate([100.0, 200.0, 300.0]):
            assert result["sigmaY"][index] == pytest.approx(published_sigma_y(x, "D"))

    def test_a_numpy_array_is_accepted(self, sigma):
        result = sigma.getSigma(np.array([100.0, 200.0]), "D", units=False)
        assert len(result["sigmaY"]) == 2

    def test_a_length_in_kilometres_is_converted(self, sigma):
        """CLAUDE.md requires pint for dimensional values; 0.1 km is 100 m."""
        from_km = sigma.getSigma(0.1 * ureg.km, "D", units=False)["sigmaY"][0]
        from_m = sigma.getSigma(100.0, "D", units=False)["sigmaY"][0]
        assert from_km == pytest.approx(from_m)

    def test_units_true_returns_quantities_in_metres(self, sigma):
        result = sigma.getSigma(100.0, "D")
        assert result["sigmaY"].units == ureg.m
        assert result["sigmaZ"].units == ureg.m
        assert result["distance"].units == ureg.m

    def test_units_true_and_false_agree_on_the_magnitude(self, sigma):
        with_units = sigma.getSigma(100.0, "D")["sigmaY"].magnitude[0]
        without = sigma.getSigma(100.0, "D", units=False)["sigmaY"][0]
        assert with_units == pytest.approx(without)

    def test_calling_the_instance_is_the_same_as_getSigma(self, sigma):
        assert sigma(100.0, "D")["sigmaY"].magnitude[0] == pytest.approx(
            sigma.getSigma(100.0, "D")["sigmaY"].magnitude[0]
        )

    def test_the_module_singleton_is_usable(self):
        assert isinstance(briggsRural, BriggsRural)
        assert briggsRural.getSigma(100.0, "D", units=False)["sigmaY"][
            0
        ] == pytest.approx(published_sigma_y(100.0, "D"))

    def test_every_expected_key_is_present(self, sigma):
        assert set(sigma.getSigma(100.0, "D", units=False)) == {
            "sigmaX",
            "sigmaY",
            "sigmaZ",
            "distance",
        }


@pytest.mark.unit
class TestErrorHandling:
    def test_the_abstract_base_refuses_to_compute(self):
        with pytest.raises(NotImplementedError, match="not implemented for abstractSigma"):
            AbstractSigma().getSigma(100.0, "D")

    @pytest.mark.parametrize("stability", ["G", "Z", "d", "", "AB"])
    def test_an_unknown_stability_class_raises(self, sigma, stability):
        """Documented as A-F in capitals; anything else must not be guessed at."""
        with pytest.raises(KeyError):
            sigma.getSigma(100.0, stability, units=False)

    @pytest.mark.xfail(
        strict=True,
        reason="B24: getSigma does not validate the distance. A negative x yields "
               "a negative standard deviation (x=-100, D -> sigma_y = -8.04 m), and "
               "a large negative x yields NaN with a numpy RuntimeWarning, because "
               "(1 + b x) goes negative under a fractional power. "
               "See the consolidated findings issue.",
    )
    @pytest.mark.parametrize("distance", [-1.0, -100.0, -20000.0])
    def test_an_upwind_distance_is_rejected(self, sigma, distance):
        """A standard deviation cannot be negative, and NaN is worse than an error.

        Either the distance should be rejected or the result should stay a
        positive real number; silently returning -8.04 m or NaN lets a bad
        input propagate into a concentration field.
        """
        with pytest.raises((ValueError, KeyError)):
            sigma.getSigma(distance, "D", units=False)


@pytest.mark.unit
class TestAbstractSigmaGetSigmaIsUnimplemented:
    def test_calling_it_on_the_base_class_raises_notimplementederror(self):
        with pytest.raises(NotImplementedError, match="not implemented for abstractSigma"):
            AbstractSigma().getSigma(100.0, "A")
