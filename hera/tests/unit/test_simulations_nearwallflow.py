"""Near-wall flow: rough-channel relations from Schlichting, chapters 16-17.

The quantities here are dimensionless groups, so the assertions are about
dimensional consistency and known limits as much as about numbers:

* a Reynolds number is dimensionless, by definition
* the technical roughness is 3.5 Ra (Schlichting p.534)
* C+ is 5 in the hydraulically smooth limit and 8 in the fully rough one
"""
import numpy as np
import pytest

from hera.simulations.hydrodynamics.nearWallFlow import (
    channelFlow,
    couetteFlow,
    functionG,
    nearWallFlow,
)
from hera.utils.unitHandler import ureg

AIR_VISCOSITY = 1.5e-5 * ureg.m**2 / ureg.s


@pytest.fixture()
def channel():
    """A 2 m channel with 0.1 mm roughness in air."""
    return channelFlow(
        Ra=1e-4 * ureg.m, nu=AIR_VISCOSITY, channelHeight=2 * ureg.m
    )


@pytest.mark.unit
class TestRoughness:
    def test_technical_roughness_is_three_and_a_half_times_ra(self, channel):
        """Schlichting p.534: k = 3.5 Ra."""
        assert channel.technical_roughness.to(ureg.m).magnitude == pytest.approx(3.5e-4)

    def test_it_keeps_length_units(self, channel):
        assert channel.technical_roughness.check("[length]")

    def test_the_dimensionless_roughness_is_a_bare_number(self, channel):
        """k+ = u* k / nu carries no units by construction."""
        value = channel.technical_roughness_plus(0.3 * ureg.m / ureg.s)
        assert isinstance(value, float)
        assert value == pytest.approx(0.3 * 3.5e-4 / 1.5e-5)

    def test_it_grows_with_friction_velocity(self, channel):
        slow = channel.technical_roughness_plus(0.1 * ureg.m / ureg.s)
        fast = channel.technical_roughness_plus(1.0 * ureg.m / ureg.s)
        assert fast > slow

    def test_a_rougher_wall_gives_a_larger_dimensionless_roughness(self):
        smooth = channelFlow(Ra=1e-6 * ureg.m, nu=AIR_VISCOSITY, channelHeight=2 * ureg.m)
        rough = channelFlow(Ra=1e-3 * ureg.m, nu=AIR_VISCOSITY, channelHeight=2 * ureg.m)
        ustar = 0.3 * ureg.m / ureg.s
        assert rough.technical_roughness_plus(ustar) > smooth.technical_roughness_plus(
            ustar
        )


@pytest.mark.unit
class TestRoughnessCorrection:
    """C+ (eqn 17.40) has two known plateaus and a log ramp between them."""

    def test_the_hydraulically_smooth_limit_is_five(self):
        smooth = channelFlow(
            Ra=1e-9 * ureg.m, nu=AIR_VISCOSITY, channelHeight=2 * ureg.m
        )
        assert smooth.Cplus(0.3 * ureg.m / ureg.s) == 5

    def test_the_fully_rough_limit_is_eight(self):
        rough = channelFlow(Ra=1.0 * ureg.m, nu=AIR_VISCOSITY, channelHeight=2 * ureg.m)
        assert rough.Cplus(0.3 * ureg.m / ureg.s) == 8

    def test_the_transition_follows_the_documented_log_law(self, channel):
        """Between k+ = 5 and 70: C+ = 8 - log(3.4 + k+)/kappa."""
        ustar = 0.3 * ureg.m / ureg.s
        kplus = channel.technical_roughness_plus(ustar)
        assert 5 < kplus < 70, "fixture no longer sits in the transition band"

        expected = 8 - (1 / 0.41) * np.log(3.4 + kplus)
        assert channel.Cplus(ustar) == pytest.approx(expected)

    def test_the_transition_is_bounded_by_the_two_plateaus(self, channel):
        value = channel.Cplus(0.3 * ureg.m / ureg.s)
        assert min(5, 8) - 5 <= value <= 8


@pytest.mark.unit
class TestFrictionReynoldsNumber:
    def test_it_is_dimensionless(self, channel):
        """Re_tau uses the channel height as a Quantity, so units cancel."""
        assert isinstance(channel.Re_tau(0.3 * ureg.m / ureg.s), float)

    def test_it_matches_its_definition(self, channel):
        """Re_tau = u* (H/2) / nu, with H the full channel height."""
        expected = 0.3 * (2.0 / 2) / 1.5e-5
        assert channel.Re_tau(0.3 * ureg.m / ureg.s) == pytest.approx(expected)

    def test_it_scales_linearly_with_friction_velocity(self, channel):
        one = channel.Re_tau(0.2 * ureg.m / ureg.s)
        two = channel.Re_tau(0.4 * ureg.m / ureg.s)
        assert two / one == pytest.approx(2.0)


@pytest.mark.unit
class TestBulkReynoldsNumber:
    @pytest.mark.xfail(
        strict=True,
        reason="B34: ReynoldsUm divides by a viscosity that carries units while "
               "multiplying by hydraulicHeight, which does not -- the height "
               "property returns _channelHeight.m_as(m), a bare float. The result "
               "comes back as 80000/meter instead of a dimensionless Reynolds "
               "number. Re_tau in the same class avoids it by using "
               "_channelHeight directly. See the consolidated findings issue.",
    )
    def test_a_reynolds_number_is_dimensionless(self, channel):
        value = channel.ReynoldsUm(0.3 * ureg.m / ureg.s, 5 * ureg.m / ureg.s)
        assert not hasattr(value, "units") or value.dimensionless

    @pytest.mark.xfail(
        strict=True,
        reason="B35: ReynoldsUm is documented as 'the reynolds based on the Um and "
               "hydraulic diameter' but its body never reads Um -- it computes "
               "ustar * hydraulicHeight / nu. Changing Um leaves the answer "
               "identical. See the consolidated findings issue.",
    )
    def test_the_bulk_velocity_changes_the_answer(self, channel):
        ustar = 0.3 * ureg.m / ureg.s
        slow = channel.ReynoldsUm(ustar, 1 * ureg.m / ureg.s)
        fast = channel.ReynoldsUm(ustar, 100 * ureg.m / ureg.s)
        assert slow != fast

    def test_the_height_property_strips_its_units(self, channel):
        """Characterisation of B34's mechanism, so the cause is unambiguous."""
        assert not hasattr(channel.height, "units"), (
            "height carries units after all; B34's diagnosis needs revisiting"
        )
        assert channel.height == pytest.approx(2.0)
        assert channel.hydraulicHeight == pytest.approx(4.0)


@pytest.mark.unit
class TestSkinFriction:
    @pytest.mark.xfail(
        strict=True,
        reason="B33: both channelFlow.skin_friction (line 174) and "
               "couetteFlow.skin_friction (line 298) contain `numpy.log()` with no "
               "argument, so every call raises TypeError. Neither method can ever "
               "have run. See the consolidated findings issue.",
    )
    def test_channel_skin_friction_is_computable(self, channel):
        value = channel.skin_friction(0.3 * ureg.m / ureg.s, 5 * ureg.m / ureg.s)
        assert value > 0

    @pytest.mark.xfail(
        strict=True,
        reason="B33: the same missing argument in couetteFlow.skin_friction. "
               "See the consolidated findings issue.",
    )
    def test_couette_skin_friction_is_computable(self):
        flow = couetteFlow(
            Ra=1e-4 * ureg.m, nu=AIR_VISCOSITY, channelHeight=2 * ureg.m
        )
        assert flow.skin_friction(0.3 * ureg.m / ureg.s, 5 * ureg.m / ureg.s) > 0


@pytest.mark.unit
class TestMeanVelocityFromFriction:
    def test_it_returns_a_velocity(self, channel):
        result = channel.get_Umean_from_Ustar(0.3 * ureg.m / ureg.s)
        assert result.check("[length]/[time]")

    def test_the_mean_velocity_exceeds_the_friction_velocity(self, channel):
        """In a turbulent channel u_m is roughly 20x u*, never below it."""
        ustar = 0.3
        mean = channel.get_Umean_from_Ustar(ustar * ureg.m / ureg.s)
        assert mean.to(ureg.m / ureg.s).magnitude > ustar

    def test_it_increases_with_friction_velocity(self, channel):
        slow = channel.get_Umean_from_Ustar(0.1 * ureg.m / ureg.s)
        fast = channel.get_Umean_from_Ustar(0.5 * ureg.m / ureg.s)
        assert fast > slow


@pytest.mark.unit
class TestFunctionG:
    """G solves an implicit relation; the tests check it actually converges."""

    def test_it_returns_a_finite_root(self):
        value = functionG().solve(1.0, 0.5)
        assert np.isfinite(value)

    def test_the_root_satisfies_the_implicit_equation(self):
        """A root finder is only correct if the residual is zero at the answer."""
        solver = functionG()
        root = solver.solve(1.0, 0.5)
        assert solver._implicitG(root, 1.0, 0.5) == pytest.approx(0.0, abs=1e-6)

    @pytest.mark.parametrize("lam, d", [(0.5, 0.1), (1.0, 0.5), (2.0, 1.0)])
    def test_it_converges_across_the_parameter_range(self, lam, d):
        solver = functionG()
        root = solver.solve(lam, d)
        assert np.isfinite(root)
        assert solver._implicitG(root, lam, d) == pytest.approx(0.0, abs=1e-6)


@pytest.mark.unit
class TestConstruction:
    def test_the_viscosity_is_stored_with_units(self, channel):
        assert channel.kinematicViscosity.check("[length]**2/[time]")

    def test_plain_numbers_are_given_the_default_units(self):
        """The docstrings promise mks defaults for bare floats."""
        flow = channelFlow(Ra=1e-4, nu=1.5e-5, channelHeight=2)
        assert flow.kinematicViscosity.to(
            ureg.m**2 / ureg.s
        ).magnitude == pytest.approx(1.5e-5)

    def test_couette_is_a_near_wall_flow(self):
        flow = couetteFlow(Ra=1e-4 * ureg.m, nu=AIR_VISCOSITY, channelHeight=2 * ureg.m)
        assert isinstance(flow, nearWallFlow)

    def test_channel_is_a_near_wall_flow(self, channel):
        assert isinstance(channel, nearWallFlow)

    def test_the_two_share_the_roughness_treatment(self, channel):
        """Both inherit technical_roughness, so they must agree on it."""
        couette = couetteFlow(
            Ra=1e-4 * ureg.m, nu=AIR_VISCOSITY, channelHeight=2 * ureg.m
        )
        assert couette.technical_roughness == channel.technical_roughness
