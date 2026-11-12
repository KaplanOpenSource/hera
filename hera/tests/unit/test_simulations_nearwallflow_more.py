"""Near-wall flow, the members `test_simulations_nearwallflow.py` left over.

That file covers the roughness treatment (`technical_roughness`,
`technical_roughness_plus`, `Cplus`), `channelFlow.Re_tau`,
`channelFlow.get_Umean_from_Ustar`, `functionG` and the broken
`skin_friction`/`ReynoldsUm` (B33-B35).  The three velocity relations it
never reached are covered here:

* ``channelFlow.get_Ucenter_from_Ustar``  (Schlichting Eqn. 17.53)
* ``couetteFlow.get_Um_from_Ustar``      (Schlichting Eqn. 17.21)
* ``couetteFlow.Re_tau``

The assertions come from the log law itself rather than from the code:

    u+ = (1/kappa) ln(Re_tau) + C+ + <additive constant>

so (a) every profile must be a velocity, (b) the centreline speed must
exceed the bulk mean, (c) the gap between the two is fixed by the two
documented additive constants of Eqn. 17.91 (C_bar = 0.94 and
C_bar + C_bar_bar = -1.7) and is therefore independent of Reynolds
number, and (d) doubling Re_tau must raise u/u* by exactly ln(2)/kappa
-- the defining slope of the logarithmic layer.

Two defects surfaced while probing:

* B243: ``Cplus``'s fully rough branch (k+ >= 70) returns the bare
  constant 8.  That constant belongs to the log law written in wall
  units of the roughness height, u+ = ln(y/k)/kappa + 8.5; the law here
  is written in viscous units (``1/kappa*log(Re_tau) + Cplus``), so the
  fully rough value has to be 8 - ln(k+)/kappa -- exactly what the
  transitional formula on the line above tends to for large k+.  With
  the term missing, C+ jumps from -2.47 to +8 across k+ = 70, and a
  fully rough wall is predicted to carry a *faster* flow than a smooth
  one at the same friction velocity.  Both velocity relations covered
  here inherit it.  Note the pre-existing
  ``test_the_fully_rough_limit_is_eight`` records that plateau as
  intended; the two monotonicity tests below say why it cannot be.
* B12 (already pinned in ``test_utils_units.py``) reaches this module:
  every method documents ``ustar: float`` with mks defaults, but a bare
  float goes through ``tounit``, which builds a Quantity in pint's
  application registry instead of hera's, so it cannot be multiplied by
  ``_channelHeight``.

``skin_friction`` is not touched here -- it is dead on arrival (B33)
and already pinned.
"""
import numpy
import pytest

from hera.simulations.hydrodynamics.nearWallFlow import channelFlow, couetteFlow
from hera.utils.unitHandler import ureg

KAPPA = 0.41
AIR_VISCOSITY = 1.5e-5 * ureg.m**2 / ureg.s
HEIGHT = 2 * ureg.m

CPLUS_REASON = (
    "B243: Cplus's fully rough branch (k+ >= 70) returns the bare constant 8. "
    "That is the Nikuradse constant of the log law written in roughness units "
    "(u+ = ln(y/k)/kappa + 8.5); this module writes the law in viscous units "
    "(1/kappa*log(Re_tau) + Cplus), so the fully rough value must be "
    "8 - ln(k+)/kappa -- the large-k+ limit of the transitional formula on the "
    "line above. With the roughness logarithm missing, C+ rises from 5 (smooth) "
    "to 8 (fully rough), so roughening a wall ADDS momentum to the flow."
)


@pytest.fixture()
def smoothChannel():
    """Hydraulically smooth, so C+ is pinned at 5 and cannot drift with u*."""
    return channelFlow(Ra=1e-9 * ureg.m, nu=AIR_VISCOSITY, channelHeight=HEIGHT)


@pytest.fixture()
def smoothCouette():
    return couetteFlow(Ra=1e-9 * ureg.m, nu=AIR_VISCOSITY, channelHeight=HEIGHT)


def _speed(quantity):
    return quantity.to(ureg.m / ureg.s).magnitude


@pytest.mark.unit
class TestCouetteFrictionReynoldsNumber:
    def test_it_matches_its_definition(self, smoothCouette):
        """Re_tau = u* (H/2) / nu, with H the full channel height."""
        assert smoothCouette.Re_tau(0.3 * ureg.m / ureg.s) == pytest.approx(
            0.3 * (2.0 / 2) / 1.5e-5
        )

    def test_it_is_a_bare_number(self, smoothCouette):
        """A Reynolds number is dimensionless by definition."""
        assert isinstance(smoothCouette.Re_tau(0.3 * ureg.m / ureg.s), float)

    def test_it_scales_linearly_with_friction_velocity(self, smoothCouette):
        one = smoothCouette.Re_tau(0.2 * ureg.m / ureg.s)
        two = smoothCouette.Re_tau(0.4 * ureg.m / ureg.s)
        assert two / one == pytest.approx(2.0)

    @pytest.mark.xfail(
        strict=True,
        reason="B12 reaching a new call site: Re_tau documents `ustar: float, "
               "unum -- default units [m/s]`, but tounit builds a bare "
               "pint Quantity for a plain number, so it lands in pint's "
               "application registry while _channelHeight (built from a ureg "
               "Quantity) stays in hera's. Multiplying them raises "
               "ValueError('Cannot operate with Quantity and Quantity of "
               "different registries'). Every documented-float entry point of "
               "this module is affected the same way. "
               "See the consolidated findings issue.",
    )
    def test_a_plain_float_is_read_as_metres_per_second(self, smoothCouette):
        """The docstring promises mks defaults for bare numbers."""
        assert smoothCouette.Re_tau(0.3) == pytest.approx(
            smoothCouette.Re_tau(0.3 * ureg.m / ureg.s)
        )

    def test_a_plain_float_currently_raises_a_registry_error(self, smoothCouette):
        """Characterisation of B12 at this call site."""
        with pytest.raises(ValueError, match="different registries"):
            smoothCouette.Re_tau(0.3)

    def test_the_two_flows_agree_on_it(self, smoothChannel, smoothCouette):
        """channelFlow and couetteFlow define Re_tau identically."""
        ustar = 0.3 * ureg.m / ureg.s
        assert smoothCouette.Re_tau(ustar) == pytest.approx(smoothChannel.Re_tau(ustar))


@pytest.mark.unit
class TestChannelCentrelineVelocity:
    def test_it_returns_a_velocity(self, smoothChannel):
        assert smoothChannel.get_Ucenter_from_Ustar(0.3 * ureg.m / ureg.s).check(
            "[length]/[time]"
        )

    def test_the_centreline_speed_exceeds_the_friction_velocity(self, smoothChannel):
        """u_c/u* is order 20 in a turbulent channel, never below one."""
        assert _speed(smoothChannel.get_Ucenter_from_Ustar(0.3 * ureg.m / ureg.s)) > 0.3

    def test_it_increases_with_friction_velocity(self, smoothChannel):
        slow = smoothChannel.get_Ucenter_from_Ustar(0.1 * ureg.m / ureg.s)
        fast = smoothChannel.get_Ucenter_from_Ustar(0.5 * ureg.m / ureg.s)
        assert fast > slow

    def test_the_centreline_exceeds_the_bulk_mean(self, smoothChannel):
        """A channel profile peaks in the middle, so u_c > u_m always."""
        ustar = 0.3 * ureg.m / ureg.s
        assert _speed(smoothChannel.get_Ucenter_from_Ustar(ustar)) > _speed(
            smoothChannel.get_Umean_from_Ustar(ustar)
        )

    @pytest.mark.parametrize("ustar", [0.05, 0.3, 1.0])
    def test_the_gap_to_the_mean_is_set_by_the_documented_constants(
        self, smoothChannel, ustar
    ):
        """Eqn. 17.91 supplies C_bar = 0.94 for the centreline and
        C_bar + C_bar_bar = -1.7 for the mean, and the two profiles are
        otherwise the same expression, so

            (u_c - u_m) / u* = 0.94 - (-1.7) = 2.64

        for every Reynolds number."""
        quantity = ustar * ureg.m / ureg.s
        difference = _speed(
            smoothChannel.get_Ucenter_from_Ustar(quantity)
        ) - _speed(smoothChannel.get_Umean_from_Ustar(quantity))
        assert difference / ustar == pytest.approx(2.64)

    def test_doubling_the_reynolds_number_follows_the_log_law_slope(
        self, smoothChannel
    ):
        """u+/u* rises by ln(2)/kappa per doubling of Re_tau. C+ is fixed at
        5 here (the wall is hydraulically smooth), so nothing else moves."""
        low = 0.2
        high = 0.4
        rise = _speed(
            smoothChannel.get_Ucenter_from_Ustar(high * ureg.m / ureg.s)
        ) / high - _speed(
            smoothChannel.get_Ucenter_from_Ustar(low * ureg.m / ureg.s)
        ) / low
        assert rise == pytest.approx(numpy.log(2.0) / KAPPA)

    @pytest.mark.xfail(
        strict=True,
        reason=CPLUS_REASON + " Consequence: at u* = 0.3 m/s a wall 10^6 times "
               "rougher carries a centreline speed 0.9 m/s HIGHER than a smooth "
               "one (9.93 vs 9.03 m/s). See the consolidated findings issue.",
    )
    def test_a_rougher_wall_slows_the_centreline(self):
        """Roughness always retards the flow at a fixed friction velocity."""
        ustar = 0.3 * ureg.m / ureg.s
        smooth = channelFlow(Ra=1e-9 * ureg.m, nu=AIR_VISCOSITY, channelHeight=HEIGHT)
        rough = channelFlow(Ra=2e-3 * ureg.m, nu=AIR_VISCOSITY, channelHeight=HEIGHT)
        assert _speed(rough.get_Ucenter_from_Ustar(ustar)) < _speed(
            smooth.get_Ucenter_from_Ustar(ustar)
        )

    def test_a_rougher_wall_currently_speeds_the_centreline_up(self):
        """Characterisation of B243, in velocity terms."""
        ustar = 0.3 * ureg.m / ureg.s
        smooth = channelFlow(Ra=1e-9 * ureg.m, nu=AIR_VISCOSITY, channelHeight=HEIGHT)
        rough = channelFlow(Ra=2e-3 * ureg.m, nu=AIR_VISCOSITY, channelHeight=HEIGHT)
        difference = _speed(rough.get_Ucenter_from_Ustar(ustar)) - _speed(
            smooth.get_Ucenter_from_Ustar(ustar)
        )
        # (8 - 5) * u*, i.e. exactly the jump in the two C+ plateaus.
        assert difference == pytest.approx((8 - 5) * 0.3)


@pytest.mark.unit
class TestCouetteMeanVelocity:
    def test_it_returns_a_velocity(self, smoothCouette):
        assert smoothCouette.get_Um_from_Ustar(0.3 * ureg.m / ureg.s).check(
            "[length]/[time]"
        )

    def test_the_mean_speed_exceeds_the_friction_velocity(self, smoothCouette):
        assert _speed(smoothCouette.get_Um_from_Ustar(0.3 * ureg.m / ureg.s)) > 0.3

    def test_it_increases_with_friction_velocity(self, smoothCouette):
        slow = smoothCouette.get_Um_from_Ustar(0.1 * ureg.m / ureg.s)
        fast = smoothCouette.get_Um_from_Ustar(0.5 * ureg.m / ureg.s)
        assert fast > slow

    def test_doubling_the_reynolds_number_follows_the_log_law_slope(
        self, smoothCouette
    ):
        low, high = 0.2, 0.4
        rise = _speed(
            smoothCouette.get_Um_from_Ustar(high * ureg.m / ureg.s)
        ) / high - _speed(
            smoothCouette.get_Um_from_Ustar(low * ureg.m / ureg.s)
        ) / low
        assert rise == pytest.approx(numpy.log(2.0) / KAPPA)

    def test_it_carries_no_additive_constant_beyond_cplus(
        self, smoothChannel, smoothCouette
    ):
        """Eqn. 17.21 has no C_bar term at all, while the channel
        centreline (17.53) adds C_bar = 0.94, so the two differ by exactly
        0.94 u* for the same wall and the same u*."""
        ustar = 0.3
        quantity = ustar * ureg.m / ureg.s
        difference = _speed(
            smoothChannel.get_Ucenter_from_Ustar(quantity)
        ) - _speed(smoothCouette.get_Um_from_Ustar(quantity))
        assert difference / ustar == pytest.approx(0.94)

    @pytest.mark.xfail(
        strict=True,
        reason=CPLUS_REASON + " The Couette profile inherits the same Cplus, so "
               "a fully rough wall speeds it up too. "
               "See the consolidated findings issue.",
    )
    def test_a_rougher_wall_slows_the_flow(self):
        ustar = 0.3 * ureg.m / ureg.s
        smooth = couetteFlow(Ra=1e-9 * ureg.m, nu=AIR_VISCOSITY, channelHeight=HEIGHT)
        rough = couetteFlow(Ra=2e-3 * ureg.m, nu=AIR_VISCOSITY, channelHeight=HEIGHT)
        assert _speed(rough.get_Um_from_Ustar(ustar)) < _speed(
            smooth.get_Um_from_Ustar(ustar)
        )


@pytest.mark.unit
class TestFullyRoughPlateau:
    """B243: the mechanism behind the two xfails above."""

    @staticmethod
    def _flowWithKplus(kplus, ustar=0.3):
        """Pick Ra so that k+ = u* (3.5 Ra) / nu takes the wanted value."""
        Ra = kplus * 1.5e-5 / (3.5 * ustar)
        return channelFlow(
            Ra=Ra * ureg.m, nu=AIR_VISCOSITY, channelHeight=HEIGHT
        )

    @pytest.mark.xfail(
        strict=True,
        reason=CPLUS_REASON + " C+ is therefore discontinuous at k+ = 70: it "
               "reaches -2.47 just below and jumps to +8 just above. "
               "See the consolidated findings issue.",
    )
    def test_cplus_is_continuous_across_the_fully_rough_threshold(self):
        """The transitional formula and the fully rough one describe the same
        curve, so they must agree where they meet."""
        ustar = 0.3 * ureg.m / ureg.s
        below = self._flowWithKplus(69.9).Cplus(ustar)
        above = self._flowWithKplus(70.1).Cplus(ustar)
        assert above == pytest.approx(below, abs=0.1)

    @pytest.mark.xfail(
        strict=True,
        reason=CPLUS_REASON + " See the consolidated findings issue.",
    )
    def test_cplus_decreases_as_the_wall_roughens(self):
        """C+ is a monotonically decreasing function of k+ over the whole
        range; roughness cannot add momentum to the flow."""
        ustar = 0.3 * ureg.m / ureg.s
        values = [
            self._flowWithKplus(kplus).Cplus(ustar)
            for kplus in (1.0, 10.0, 50.0, 100.0, 1000.0)
        ]
        assert values == sorted(values, reverse=True)

    def test_cplus_currently_jumps_to_the_nikuradse_constant(self):
        """Characterisation of B243: the value above k+ = 70 is the bare
        constant 8, ten and a half units above the transitional curve it is
        supposed to continue."""
        ustar = 0.3 * ureg.m / ureg.s
        below = self._flowWithKplus(69.9).Cplus(ustar)
        above = self._flowWithKplus(70.1).Cplus(ustar)
        assert above == 8
        assert below == pytest.approx(8 - numpy.log(3.4 + 69.9) / KAPPA)
        assert above - below == pytest.approx(10.47, abs=0.01)

    def test_the_missing_term_is_the_roughness_logarithm(self):
        """Characterisation of B243's fix: because the surrounding log law
        is written in viscous units (ln Re_tau), the fully rough constant has
        to be 8 - ln(k+)/kappa, which is what the transitional formula tends
        to for large k+ -- the two differ by only 0.1 at k+ = 1000."""
        ustar = 0.3 * ureg.m / ureg.s
        kplus = 1000.0
        transitional = 8 - numpy.log(3.4 + kplus) / KAPPA
        expected = 8 - numpy.log(kplus) / KAPPA
        assert transitional == pytest.approx(expected, abs=0.1)
        assert self._flowWithKplus(kplus).Cplus(ustar) == 8
