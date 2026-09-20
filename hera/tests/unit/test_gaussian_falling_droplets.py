"""simulations/gaussian/FallingNonEvaporatingDroplets.py: a gaussian model
for droplets falling at terminal velocity.

B96: the constructor is completely broken --
``self._meteorology = MeteorologyFactory().getMeteorology(meteorologyName, **met_kwargs)``
calls a method that does not exist. ``MeteorologyFactory`` only defines
``getMeteorologyFromU10`` and ``getMeteorologyFromURefHeight``, each with
a different (non-string-name) parameter list. Every single instantiation
of ``FallingNonEvaporatingDroplets(...)`` raises ``AttributeError`` before
anything else runs -- the class is entirely unusable through its public
constructor.

Because __init__ never completes, the tests below construct instances via
``FallingNonEvaporatingDroplets.__new__(...)`` and set only the attributes
each method actually touches, to independently verify the physics that
would otherwise work correctly once the constructor is fixed. The
``buildDroplet`` fixture goes one step further and also injects a real
``MeteorologyFactory`` profile (no database, no network), which makes the
settling physics reachable: getTerminalVelocity, _VTFunc,
correctionCloud_Plume/_Puff, solveToTime and _fallingParticle.

Three further findings came out of exercising those:

* B119 (terminal velocity): ``_VTFunc`` computes the Reynolds number as
  ``(Vt * d / mu).magnitude`` -- the air density is missing, so the
  expression is not dimensionless and ``.magnitude`` silently yields a
  number 10/rho_air ~ 8.3 times the true Reynolds number.
  ``_fallingParticle`` does the same calculation correctly
  (``rho_air * Uabs * d / nu_air``), which shows the intent. The error
  cancels exactly in the Stokes regime (where Cd*Re = 24 is constant) but
  not above it: a 3 mm droplet is reported to fall at 2.83 m/s instead of
  the ~8.2 m/s a 3 mm drop actually falls at.
* B120 (cloud correction): ``correctionCloud_Plume`` / ``_Puff`` average a
  *python list of pint Quantities* with ``numpy.mean``, which always
  raises ValueError. Neither function can ever return a value, and
  "Plume" is the default, so ``solveToTime`` is unusable too.
* B12 (already documented) cascades into ``_fallingParticle``: ``tounit``
  builds its Quantity from the bare ``pint`` import, i.e. in the default
  application registry, so ``u_p - U`` mixes two registries and raises.
  ``solveToTime`` therefore dies inside ``solve_ivp``.

Deliberately left uncovered: lines 5-6 (the ImportError fallback for a
scipy without ``solve_ivp``, unreachable with a working scipy) and the
body of ``__init__`` (unreachable while B96 stands, and already pinned
below).
"""
import numpy
import pandas
import pytest
from pint import DimensionalityError

from hera.simulations.gaussian.FallingNonEvaporatingDroplets import (
    FallingNonEvaporatingDroplets,
    hit_ground,
)
from hera.simulations.gaussian.Meteorology import MeteorologyFactory
from hera.utils import ureg

STANDARD_GRAVITY = 9.80665


@pytest.fixture()
def droplet():
    return FallingNonEvaporatingDroplets.__new__(FallingNonEvaporatingDroplets)


@pytest.fixture()
def buildDroplet():
    """Build a *usable* droplet without running the B96-broken constructor.

    Every attribute __init__ would have set is set here by hand -- the
    established technique in this suite for a class whose constructor
    cannot run.  The meteorology is a real object built by the real
    MeteorologyFactory (pure arithmetic: no database, no network), so the
    settling physics is exercised against exactly the profiles __init__
    would have handed it.
    """
    def _build(diameter=100 * ureg.um,
               height=0 * ureg.m,
               u10=5 * ureg.m / ureg.s,
               drag="Haugen",
               correction="Plume",
               profile="uniformWind",
               stability="D",
               rho_p=0.9 * ureg.g / ureg.cm ** 3,
               Q=1 * ureg.kg):
        built = FallingNonEvaporatingDroplets.__new__(FallingNonEvaporatingDroplets)
        built._meteorology = MeteorologyFactory().getMeteorologyFromU10(
            u10=u10,
            inversion=1000 * ureg.m,
            verticalProfileType=profile,
            stability=stability,
        )
        built._dragfunc = getattr(built, "_DragCoefficient_%s" % drag)
        built._correctionfunc = getattr(built, "correctionCloud_%s" % correction)
        built.rho_p = rho_p
        built.particleDiameter = diameter
        built.Q = Q
        built.cloudSigma = (0 * ureg.m, 0 * ureg.m, 0 * ureg.m)
        built.position = (0 * ureg.m, 0 * ureg.m, height)
        return built

    return _build


def airProperties(droplet, height=0 * ureg.m):
    """rho_air [kg/m3] and mu [Pa s], read straight off the meteorology.

    Both are independently checked in test_gaussian_meteorology.py against
    the textbook values for dry air at 20 C (1.204 kg/m3, ~1.81e-5 Pa s),
    so using them here keeps the expected settling velocities free of any
    number taken from FallingNonEvaporatingDroplets itself.
    """
    rho_air = droplet.meteorology.getAirDensity(height).m_as(ureg.kg / ureg.m ** 3)
    mu = droplet.meteorology.getAirDynamicViscosity(height).m_as(ureg.Pa * ureg.s)
    return rho_air, mu


@pytest.mark.unit
class TestConstructorIsBroken:
    def test_meteorology_factory_has_no_getmeteorology_method(self):
        """Characterisation of B96, at its exact source."""
        assert not hasattr(MeteorologyFactory(), "getMeteorology")

    @pytest.mark.xfail(
        strict=True,
        reason="B96: __init__ calls MeteorologyFactory().getMeteorology(...), "
               "a method that does not exist on MeteorologyFactory (only "
               "getMeteorologyFromU10/getMeteorologyFromURefHeight are "
               "defined). Every construction fails. "
               "See the consolidated findings issue.",
    )
    def test_it_should_be_constructible(self):
        FallingNonEvaporatingDroplets(particleDiameter=100 * ureg.um)

    def test_it_currently_raises_on_every_construction(self):
        """Characterisation of B96."""
        with pytest.raises(AttributeError, match="getMeteorology"):
            FallingNonEvaporatingDroplets(particleDiameter=100 * ureg.um)


@pytest.mark.unit
class TestPositionAndSigmaProperties:
    def test_position_exposes_x_y_z(self, droplet):
        droplet.position = (1 * ureg.m, 2 * ureg.m, 3 * ureg.m)
        assert droplet.x == 1 * ureg.m
        assert droplet.y == 2 * ureg.m
        assert droplet.z == 3 * ureg.m

    def test_position_rejects_a_non_3_tuple(self, droplet):
        with pytest.raises(ValueError, match="Must be a 3 tuple"):
            droplet.position = (1 * ureg.m, 2 * ureg.m)

    def test_cloud_sigma_exposes_x_y_z(self, droplet):
        droplet.cloudSigma = (1 * ureg.m, 2 * ureg.m, 3 * ureg.m)
        assert droplet.cloudSigmaX == 1 * ureg.m
        assert droplet.cloudSigmaY == 2 * ureg.m
        assert droplet.cloudSigmaZ == 3 * ureg.m

    def test_cloud_sigma_rejects_a_non_3_tuple(self, droplet):
        with pytest.raises(ValueError, match="Must be a 3 tuple"):
            droplet.cloudSigma = (1 * ureg.m,)


@pytest.mark.unit
class TestParticleMassDerivation:
    def test_particle_mass_is_none_until_both_diameter_and_density_are_set(self, droplet):
        droplet.rho_p = 0.9 * ureg.g / ureg.cm ** 3
        assert droplet.particleMass is None

    def test_setting_both_computes_the_sphere_mass(self, droplet):
        droplet.rho_p = 1000 * ureg.kg / ureg.m ** 3
        droplet.particleDiameter = 1000 * ureg.um  # 1mm sphere, water density
        # V = pi/6 * d^3 ; d=1mm -> V ~ 5.236e-10 m^3 ; m = rho*V
        assert droplet.particleMass.m_as(ureg.kg) == pytest.approx(5.236e-10 * 1000, rel=1e-3)


@pytest.mark.unit
class TestSimpleProperties:
    def test_beta_is_a_fixed_constant(self, droplet):
        assert droplet.beta == 5.0

    def test_spread_factor_is_a_fixed_constant(self, droplet):
        assert droplet.SpreadFactor == 4.5

    def test_g_is_standard_gravity_with_units(self, droplet):
        assert droplet.g.m_as(ureg.m / ureg.s ** 2) == pytest.approx(9.80665)

    def test_n_is_total_mass_over_particle_mass(self, droplet):
        droplet.rho_p = 1000 * ureg.kg / ureg.m ** 3
        droplet.particleDiameter = 1000 * ureg.um
        droplet._cloudQ = droplet.particleMass * 10  # exactly 10 particles' worth
        assert droplet.N == pytest.approx(10)


@pytest.mark.unit
class TestDragCoefficientIk:
    def test_very_low_reynolds_uses_stokes_law(self, droplet):
        assert droplet._DragCoefficient_Ik(0.02) == pytest.approx(24.0 / 0.02)

    def test_it_is_continuous_and_positive_across_all_branches(self, droplet):
        for re in [0.02, 1.0, 100, 2000, 3000]:
            assert droplet._DragCoefficient_Ik(re) > 0

    def test_the_two_supercritical_branches_join_smoothly_at_2905(self, droplet):
        """A piecewise fit of a physical curve must not jump at a breakpoint.

        The two linear pieces are 0.21667354 + 1.44e-4 Re and
        0.11174 + 1.80124e-4 Re; they were fitted to agree at Re = 2905, so
        a mismatch there would mean one of the constants had been mistyped.
        """
        below = droplet._DragCoefficient_Ik(2904.999)
        above = droplet._DragCoefficient_Ik(2905.001)
        assert below == pytest.approx(above, rel=1e-4)

    def test_the_correlation_joins_the_linear_branch_at_2500(self, droplet):
        """Same argument across the other breakpoint, Re = 2500."""
        below = droplet._DragCoefficient_Ik(2499.999)
        above = droplet._DragCoefficient_Ik(2500.001)
        assert below == pytest.approx(above, rel=1e-4)

    def test_the_supercritical_drag_is_of_order_one_half(self, droplet):
        """The drag crisis has not been reached yet at Re ~ 3e3.

        A sphere sits on the Cd ~ 0.4-0.6 plateau from Re ~ 1e3 to ~ 2e5.
        """
        assert 0.4 < droplet._DragCoefficient_Ik(2600) < 0.7


@pytest.mark.unit
class TestDragCoefficientKelbaliyev15:
    def test_very_low_reynolds_uses_stokes_law(self, droplet):
        assert droplet._DragCoefficient_Kelbaliyev15(0.005) == pytest.approx(24.0 / 0.005)

    def test_high_reynolds_saturates_at_0_44(self, droplet):
        assert droplet._DragCoefficient_Kelbaliyev15(20000) == pytest.approx(0.44)

    def test_drag_decreases_as_reynolds_increases(self, droplet):
        values = [droplet._DragCoefficient_Kelbaliyev15(re) for re in [5, 100, 1000, 5000]]
        assert values == sorted(values, reverse=True)


@pytest.mark.unit
class TestDragCoefficientFan:
    def test_stokes_regime_below_2(self, droplet):
        assert droplet._DragCoefficient_Fan(1) == pytest.approx(24.0)

    def test_high_reynolds_saturates_at_0_44(self, droplet):
        assert droplet._DragCoefficient_Fan(600) == pytest.approx(0.44)

    def test_intermediate_regime_uses_the_power_law(self, droplet):
        assert droplet._DragCoefficient_Fan(10) == pytest.approx(18.5 / 10 ** 0.6)


@pytest.mark.unit
class TestDragCoefficientHaugen:
    def test_high_reynolds_saturates_at_0_44(self, droplet):
        assert droplet._DragCoefficient_Haugen(2000) == pytest.approx(0.44)

    def test_low_reynolds_uses_the_correction_formula(self, droplet):
        re = 10
        assert droplet._DragCoefficient_Haugen(re) == pytest.approx(24 / re * (1 + 0.15 * re ** 0.687))


@pytest.mark.unit
class TestCorrectionCloudNone:
    def test_it_is_a_no_op_factor(self, droplet):
        assert droplet.correctionCloud_None() == 1


@pytest.mark.unit
class TestCloudMassAndFootprint:
    def test_q_is_converted_to_kilograms(self, droplet):
        droplet.Q = 500 * ureg.g
        assert droplet.Q.m_as(ureg.kg) == pytest.approx(0.5)

    def test_the_meteorology_property_returns_the_injected_profile(self, buildDroplet):
        built = buildDroplet()
        assert built.meteorology is built._meteorology

    def test_position_and_cloud_sigma_are_readable_as_triples(self, buildDroplet):
        built = buildDroplet(height=40 * ureg.m)
        assert [q.m_as(ureg.m) for q in built.position] == [0.0, 0.0, 40.0]
        assert [q.m_as(ureg.m) for q in built.cloudSigma] == [0.0, 0.0, 0.0]

    def test_the_footprint_matches_the_documented_closed_form(self, buildDroplet):
        r"""The docstring reduces N * pi (d/2)^2 * S to 3 Q S / (2 rho d).

        1 kg of oil (900 kg/m3) as 1 mm droplets, spread factor 4.5:
        3 * 1 * 4.5 / (2 * 900 * 1e-3) = 7.5 m2.
        """
        built = buildDroplet(diameter=1000 * ureg.um)
        expected = 3 * 1.0 * 4.5 / (2 * 900.0 * 1e-3)
        assert built.AreaOnSurface.m_as(ureg.m ** 2) == pytest.approx(expected)
        assert expected == pytest.approx(7.5)

    def test_the_footprint_grows_as_the_droplets_get_smaller(self, buildDroplet):
        """Splitting the same mass into finer droplets wets more ground: 1/d."""
        coarse = buildDroplet(diameter=1000 * ureg.um).AreaOnSurface.m_as(ureg.m ** 2)
        fine = buildDroplet(diameter=100 * ureg.um).AreaOnSurface.m_as(ureg.m ** 2)
        assert fine == pytest.approx(10 * coarse)

    def test_the_particle_count_is_the_mass_ratio(self, buildDroplet):
        """1 kg / (pi/6 * 900 * (1e-3)^3 kg) ~ 2.12e6 millimetric droplets."""
        built = buildDroplet(diameter=1000 * ureg.um)
        expected = 1.0 / (900.0 * numpy.pi / 6.0 * (1e-3) ** 3)
        assert built.N == pytest.approx(expected)


@pytest.mark.unit
class TestTerminalVelocity:
    def test_a_twenty_micron_droplet_settles_at_the_stokes_velocity(self, buildDroplet):
        r"""Stokes: Vt = (rho_p - rho_air) g d^2 / (18 mu).

        For 20 um oil (900 kg/m3) in air at 20 C that is about 1.06 cm/s --
        the familiar "a fine mist barely falls" number, and the one regime
        where the missing air density in _VTFunc (B119 below) does no harm,
        because Cd*Re = 24 there and the Reynolds number cancels out of the
        force balance entirely.  The Haugen correction 1 + 0.15 Re^0.687 is
        evaluated at a Reynolds number 8.3x too large, which is what leaves
        a residual 3% deficit here; hence the 5% tolerance.
        """
        built = buildDroplet(diameter=20 * ureg.um)
        rho_air, mu = airProperties(built)
        stokes = (900.0 - rho_air) * STANDARD_GRAVITY * (20e-6) ** 2 / (18 * mu)
        assert stokes == pytest.approx(0.0106, rel=0.02)
        assert built.getTerminalVelocity().m_as(ureg.m / ureg.s) == pytest.approx(
            stokes, rel=0.05
        )

    def test_terminal_velocity_grows_with_droplet_diameter(self, buildDroplet):
        """Monotone by construction: gravity scales as d^3, drag at most d^2."""
        velocities = [
            buildDroplet(diameter=d * ureg.um).getTerminalVelocity().m_as(
                ureg.m / ureg.s
            )
            for d in (10, 50, 100, 500, 1000, 3000)
        ]
        assert velocities == sorted(velocities)

    def test_terminal_velocity_grows_as_the_square_of_small_diameters(self, buildDroplet):
        """In the Stokes regime Vt ~ d^2, so 5 um -> 10 um is a factor 4.

        Deep in the Stokes regime only: the Haugen correction is evaluated
        at an inflated Reynolds number (B119), and by 50 um that alone has
        already bent the exponent noticeably away from 2.
        """
        small = buildDroplet(diameter=5 * ureg.um).getTerminalVelocity()
        larger = buildDroplet(diameter=10 * ureg.um).getTerminalVelocity()
        assert (larger / small).magnitude == pytest.approx(4.0, rel=0.02)

    def test_an_omitted_height_falls_back_to_the_cloud_height(self, buildDroplet):
        built = buildDroplet(diameter=20 * ureg.um, height=250 * ureg.m)
        assert built.getTerminalVelocity() == built.getTerminalVelocity(
            height=250 * ureg.m
        )

    def test_droplets_settle_faster_in_the_thinner_air_aloft(self, buildDroplet):
        """Air is cooler and hence less viscous with height, and Vt ~ 1/mu.

        The lapse rate is 6.5 K/km, so at 300 m the viscosity correlation
        gives about 0.6% less than at the ground -- small, but the sign is
        the physics being asserted.
        """
        built = buildDroplet(diameter=20 * ureg.um)
        assert built.getTerminalVelocity(height=300 * ureg.m) > built.getTerminalVelocity(
            height=0 * ureg.m
        )

    @pytest.mark.xfail(
        strict=True,
        reason="B119: _VTFunc computes Re as (Vt * d / nu_air).magnitude, "
               "omitting the air density -- the expression is not "
               "dimensionless and .magnitude yields 10/rho_air ~ 8.3 times "
               "the true Reynolds number, so the drag coefficient is read "
               "off the wrong part of the curve.  _fallingParticle computes "
               "the same Reynolds number correctly as rho_air*Uabs*d/nu_air. "
               "Above the Stokes regime the terminal velocity comes out "
               "sqrt(rho_air/10) ~ 0.347 of the physical value. "
               "See the consolidated findings issue.",
    )
    def test_a_three_millimetre_drop_matches_the_high_reynolds_force_balance(
        self, buildDroplet
    ):
        r"""At Re > 1000 Haugen fixes Cd = 0.44, so the balance closes:

            (1/2) Cd rho_air A Vt^2 = (rho_p - rho_air) V g
            =>  Vt = sqrt(4 g d (rho_p - rho_air) / (3 Cd rho_air))

        which for 3 mm gives 8.2 m/s -- and a 3 mm raindrop is measured to
        fall at about 8 m/s, so the formula is independently anchored.
        """
        built = buildDroplet(diameter=3000 * ureg.um)
        rho_air, _ = airProperties(built)
        expected = numpy.sqrt(
            4 * STANDARD_GRAVITY * 3e-3 * (900.0 - rho_air) / (3 * 0.44 * rho_air)
        )
        assert expected == pytest.approx(8.2, rel=0.02)
        assert built.getTerminalVelocity().m_as(ureg.m / ureg.s) == pytest.approx(
            expected, rel=0.05
        )

    def test_it_currently_returns_the_physical_value_scaled_down(self, buildDroplet):
        """Characterisation of B119.

        Because the only defect is the argument handed to Cd(Re), and Haugen
        is constant above Re = 1000, the error collapses to one clean
        factor: Re is inflated by k = 10/rho_air, the balance divides by
        Cd*Re, so the returned velocity is the physical one over sqrt(k).
        """
        built = buildDroplet(diameter=3000 * ureg.um)
        rho_air, _ = airProperties(built)
        physical = numpy.sqrt(
            4 * STANDARD_GRAVITY * 3e-3 * (900.0 - rho_air) / (3 * 0.44 * rho_air)
        )
        inflation = 10.0 / rho_air
        assert built.getTerminalVelocity().m_as(ureg.m / ureg.s) == pytest.approx(
            physical / numpy.sqrt(inflation), rel=1e-3
        )

    def test_the_reynolds_expression_without_density_is_not_dimensionless(
        self, buildDroplet
    ):
        """Characterisation of B119, at its source.

        Re = rho V d / mu is dimensionless; V d / mu is not.  pint refuses
        the conversion, and .magnitude is what hides that from the caller.
        """
        built = buildDroplet()
        rho_air_q = built.meteorology.getAirDensity(0 * ureg.m)
        mu_q = built.meteorology.getAirDynamicViscosity(0 * ureg.m)
        speed = 1 * ureg.m / ureg.s
        diameter = 1000 * ureg.um

        asWritten = speed * diameter / mu_q
        with pytest.raises(DimensionalityError):
            asWritten.m_as(ureg.dimensionless)

        asIntended = rho_air_q * speed * diameter / mu_q
        assert asIntended.m_as(ureg.dimensionless) == pytest.approx(64.99, rel=1e-3)


@pytest.mark.unit
class TestTerminalVelocityResidual:
    def test_the_residual_vanishes_at_the_solved_velocity(self, buildDroplet):
        """getTerminalVelocity returns the root of _VTFunc, by definition."""
        built = buildDroplet(diameter=100 * ureg.um)
        solved = built.getTerminalVelocity().m_as(ureg.m / ureg.s)
        assert built._VTFunc([solved]) == pytest.approx(0.0, abs=1e-9)

    def test_the_residual_is_positive_below_and_negative_above_the_root(
        self, buildDroplet
    ):
        """The balance velocity falls as the trial velocity rises, so the
        residual GuessVt - Vt crosses zero once, downwards."""
        built = buildDroplet(diameter=100 * ureg.um)
        assert built._VTFunc([1e-4]) > 0
        assert built._VTFunc([50.0]) < 0

    def test_the_residual_uses_the_requested_height(self, buildDroplet):
        built = buildDroplet(diameter=100 * ureg.um, height=0 * ureg.m)
        assert built._VTFunc([0.5], height=300 * ureg.m) != built._VTFunc([0.5])


@pytest.mark.unit
class TestCloudSigmaCorrection:
    """Csanady's correction for a settling plume: sigma is reduced by
    (1 + (beta Vt / Ubar)^2)^-1/4 for a plume and ^-1/2 for a puff, so the
    puff factor is the square of the plume factor and both lie in (0, 1].
    """

    @pytest.mark.xfail(
        strict=True,
        reason="B120: correctionCloud_Plume/_Puff call "
               "numpy.mean([...]) on a python LIST of pint Quantities "
               "(one per metre of fall height).  numpy cannot build an "
               "array from Quantity objects, so both raise "
               "ValueError('setting an array element with a sequence') "
               "for every input -- the correction can never be computed, "
               "and since 'Plume' is the default correctionCloudFunc, "
               "solveToTime is unusable as well.  numpy.mean over a single "
               "Quantity array (ureg.Quantity([...], 'm/s')) works, which "
               "is what was meant. See the consolidated findings issue.",
    )
    def test_the_puff_correction_is_the_square_of_the_plume_correction(
        self, buildDroplet
    ):
        built = buildDroplet(diameter=1000 * ureg.um, height=100 * ureg.m)
        plume = built.correctionCloud_Plume()
        puff = built.correctionCloud_Puff()
        assert 0 < plume < 1
        assert puff == pytest.approx(plume ** 2)

    def test_the_plume_correction_currently_always_raises(self, buildDroplet):
        """Characterisation of B120."""
        built = buildDroplet(diameter=1000 * ureg.um, height=100 * ureg.m)
        with pytest.raises(ValueError, match="setting an array element with a sequence"):
            built.correctionCloud_Plume()

    def test_the_puff_correction_currently_always_raises(self, buildDroplet):
        """Characterisation of B120."""
        built = buildDroplet(diameter=1000 * ureg.um, height=100 * ureg.m)
        with pytest.raises(ValueError, match="setting an array element with a sequence"):
            built.correctionCloud_Puff()

    def test_numpy_cannot_average_a_list_of_quantities(self):
        """Characterisation of B120, isolated from the class.

        The same average over one Quantity array is fine, which is the
        difference between the code as written and as intended.
        """
        with pytest.raises(ValueError):
            numpy.mean([1 * ureg.m / ureg.s, 2 * ureg.m / ureg.s])

        assert numpy.mean(ureg.Quantity([1.0, 2.0], "m/s")).m_as(
            ureg.m / ureg.s
        ) == pytest.approx(1.5)


@pytest.mark.unit
class TestFallingParticleDerivative:
    @pytest.mark.xfail(
        strict=True,
        reason="B12 cascade: tounit() builds its Quantity from the bare "
               "`from pint import Quantity` import, i.e. in pint's default "
               "application registry, while the meteorology's wind speed "
               "lives in hera's own ureg.  _fallingParticle then evaluates "
               "u_p - U across the two registries and pint raises "
               "ValueError('Cannot operate with Quantity and Quantity of "
               "different registries').  The right-hand side of the ODE "
               "can therefore never be evaluated. "
               "See the consolidated findings issue.",
    )
    def test_a_droplet_released_from_rest_starts_in_free_fall(self, buildDroplet):
        """With no initial velocity, drag vanishes on the vertical axis, so
        dw/dt = -g exactly, dz/dt = 0, and the wind starts pushing the
        droplet downwind (du_p/dt > 0)."""
        built = buildDroplet(diameter=1000 * ureg.um, height=50 * ureg.m)
        derivative = built._fallingParticle(
            0.0, numpy.array([0.0, 0.0, 50.0, 0.0, 0.0, 0.0])
        )
        assert derivative[2] == pytest.approx(0.0)
        assert derivative[5] == pytest.approx(-STANDARD_GRAVITY)
        assert derivative[4] > 0

    def test_the_derivative_currently_raises_on_mixed_registries(self, buildDroplet):
        """Characterisation of the B12 cascade."""
        built = buildDroplet(diameter=1000 * ureg.um, height=50 * ureg.m)
        with pytest.raises(ValueError, match="different registries"):
            built._fallingParticle(0.0, numpy.array([0.0, 0.0, 50.0, 0.0, 0.0, 0.0]))

    def test_numpy_max_with_a_scalar_second_argument_does_not_clip(self):
        """_fallingParticle line 402 reads `z = numpy.max(z, 0)`.

        That second positional argument is `axis`, not a floor, so the
        intended "never go below the ground" clamp is a no-op.  Recorded
        here because the line is unreachable while the B12 cascade above
        stands; numpy.clip / max(z, 0) is what it wanted.
        """
        assert numpy.max(-5.0, 0) == -5.0


@pytest.mark.unit
class TestSolveToTime:
    @pytest.mark.xfail(
        strict=True,
        reason="B12 cascade: solveToTime hands _fallingParticle to "
               "scipy.solve_ivp, and the first right-hand-side evaluation "
               "raises ValueError('Cannot operate with Quantity and "
               "Quantity of different registries') -- see "
               "TestFallingParticleDerivative.  Independently, the default "
               "correctionCloudFunc 'Plume' would raise as well (B120). "
               "See the consolidated findings issue.",
    )
    def test_it_returns_the_documented_trajectory_table(self, buildDroplet):
        """The docstring promises times plus position, Q and sigma; the
        droplet must also actually descend."""
        built = buildDroplet(
            diameter=1000 * ureg.um, height=50 * ureg.m, correction="None"
        )
        trajectory = built.solveToTime(200)
        assert isinstance(trajectory, pandas.DataFrame)
        for column in ("t", "x", "z", "distance", "u_p", "w_p", "Q", "diameter"):
            assert column in trajectory.columns
        assert trajectory["z"].iloc[-1] < trajectory["z"].iloc[0]

    def test_it_currently_raises_from_inside_the_integrator(self, buildDroplet):
        """Characterisation of the B12 cascade."""
        built = buildDroplet(
            diameter=1000 * ureg.um, height=50 * ureg.m, correction="None"
        )
        with pytest.raises(ValueError, match="different registries"):
            built.solveToTime(200)


@pytest.mark.unit
class TestHitGroundEvent:
    def test_it_is_terminal(self):
        assert hit_ground.terminal is True

    def test_it_is_zero_essentially_at_the_ground(self):
        assert hit_ground()(0, [0, 0, 1e-3, 0, 0, 0]) == pytest.approx(0, abs=1e-9)

    def test_it_is_positive_above_ground(self):
        assert hit_ground()(0, [0, 0, 5, 0, 0, 0]) > 0

    def test_it_is_negative_below_ground(self):
        assert hit_ground()(0, [0, 0, -1, 0, 0, 0]) < 0
