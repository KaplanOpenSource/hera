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
each method actually touches, to independently verify the pure-math parts
(the four drag-coefficient formulas, the property setters/getters) that
would otherwise work correctly once the constructor is fixed. Methods
that need ``self.meteorology`` (getTerminalVelocity, _VTFunc,
correctionCloud_Plume/_Puff, solveToTime, _fallingParticle) are left
uncovered -- they cannot be reached at all while B96 stands.
"""
import pytest

from hera.simulations.gaussian.FallingNonEvaporatingDroplets import (
    FallingNonEvaporatingDroplets,
    hit_ground,
)
from hera.simulations.gaussian.Meteorology import MeteorologyFactory
from hera.utils import ureg


@pytest.fixture()
def droplet():
    return FallingNonEvaporatingDroplets.__new__(FallingNonEvaporatingDroplets)


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
class TestHitGroundEvent:
    def test_it_is_terminal(self):
        assert hit_ground.terminal is True

    def test_it_is_zero_essentially_at_the_ground(self):
        assert hit_ground()(0, [0, 0, 1e-3, 0, 0, 0]) == pytest.approx(0, abs=1e-9)

    def test_it_is_positive_above_ground(self):
        assert hit_ground()(0, [0, 0, 5, 0, 0, 0]) > 0

    def test_it_is_negative_below_ground(self):
        assert hit_ground()(0, [0, 0, -1, 0, 0, 0]) < 0
