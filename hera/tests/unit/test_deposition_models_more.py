"""Dry deposition: the two paths the earlier deposition files leave alone.

``test_deposition_models.py`` covers construction with a dict surface, the
tangled ``ustar``/``diameter``/``heatFlux`` properties (B46-B49) and the
dispatch; ``test_deposition_models_canopy.py`` covers the
``obstacleShape`` branch and its nested PhiM/PsiM.  What neither reaches
is:

* the constructor's *other* surface branch -- a surface given by name
  instead of by dict, which is looked up in the project's Cache
  collection;
* the bare-ground Brownian arithmetic itself: the Cunningham slip factor,
  the Stokes-Einstein diffusivity, the Schmidt number and the ground
  deposition efficiency ``egb``, together with the in-canopy resistance
  ``Qsol`` and the nested ``PhiH`` that only the *non*-obstacle path
  exercises.

Because of B46 the surviving ``ustar`` property returns ``_diameter``, so
``diameter`` is what the model uses as friction velocity; the tests pass
``diameter=0.3`` where a realistic u* is wanted, exactly as the canopy
file does, and read ``model.ustar`` when they need the value the
calculation actually saw.

Three defects surfaced:

* B254: ``surface`` given as a name runs
  ``p.getCacheDocuments(type="surface", surface=surface)[0]``.  Nothing
  checks the list, so an unregistered surface name raises
  ``IndexError: list index out of range`` -- which names neither the
  surface nor the collection it was looked for in.
* B255: ``depositionRate_Petroff`` starts with
  ``dpm = 0.000001*dp``, a micron-to-metre conversion, but ``dp`` comes
  from ``self.diameter``, which the constructor has already normalised to
  metres (``tonumber(tounit(diameter, ureg.m), ureg.m)``).  Every other
  constant in the block is SI (mean free path 6.7e-8 m, viscosity
  1.89e-5 Pa s), so the model evaluates a 1 micron particle as a 1
  picometre one: the ground deposition efficiency ``vds/u*`` comes out at
  8.5e4 instead of order 1e-5, i.e. the deposition velocity exceeds the
  friction velocity by five orders of magnitude.
* B257: the Boltzmann constant in the same block is written
  ``kb = 1.83E-23``.  Its value is 1.380649e-23 J/K -- the digits are
  transposed -- so the Brownian diffusivity is 33% too large wherever
  the model is used.
"""
import re

import pytest

from hera.simulations.deposition.models import depositionModels
from hera.utils.unitHandler import ureg

DESERT = {"name": "Desert", "zrough": 0.04}
# A canopy WITHOUT obstacleShape: takes the bare-ground formula, but with a
# non-zero hcanop, so lm and Qsol stop being zero and PhiH is exercised.
GRASSLAND = {
    "name": "Grassland", "zrough": 0.1, "hcanop": 10.0, "displacementHeight": 6.0,
}


@pytest.fixture()
def make(unit_project):
    """depositionModels opens a real Project, so it must be pointed at the
    mongomock-backed test project rather than its own default."""

    def _make(**kwargs):
        kwargs.setdefault("surface", DESERT)
        kwargs.setdefault("diameter", 0.3)
        return depositionModels(projectName=unit_project.projectName, **kwargs)

    return _make


def _registerSurface(project, name, **fields):
    project.addCacheDocument(
        resource="", dataFormat="string", type="surface",
        desc=dict(surface=name, name=name, **fields),
    )


@pytest.mark.unit
class TestSurfaceLookedUpByName:
    def test_a_registered_surface_name_is_resolved(self, unit_project):
        """The docstring's alternative to a dict: the name is matched against
        the Cache collection's `surface` field."""
        _registerSurface(unit_project, "Sand", zrough=0.002)
        built = depositionModels(
            projectName=unit_project.projectName, surface="Sand", diameter=0.3
        )
        assert built.surface["zrough"] == pytest.approx(0.002)

    def test_the_whole_descriptor_becomes_the_surface(self, unit_project):
        """asDict()["desc"] is handed over as-is, so extra fields survive."""
        _registerSurface(unit_project, "Sand", zrough=0.002, hcanop=0.5)
        built = depositionModels(
            projectName=unit_project.projectName, surface="Sand", diameter=0.3
        )
        assert built.surface["hcanop"] == pytest.approx(0.5)
        assert built.surface["surface"] == "Sand"

    def test_a_looked_up_surface_reaches_the_calculation(self, unit_project):
        """A canopy registered by name must change the answer the same way a
        canopy passed as a dict does."""
        _registerSurface(
            unit_project, "Meadow", zrough=0.1, hcanop=10.0, displacementHeight=6.0
        )
        byName = depositionModels(
            projectName=unit_project.projectName, surface="Meadow", diameter=0.3
        ).depositionRate()
        byDict = depositionModels(
            projectName=unit_project.projectName, surface=GRASSLAND, diameter=0.3
        ).depositionRate()
        assert byName == pytest.approx(byDict)

    def test_a_dict_surface_never_touches_the_database(self, make):
        """The `type(surface)==dict` test short-circuits the lookup."""
        assert make(surface=DESERT).surface is DESERT

    @pytest.mark.xfail(
        strict=True,
        reason="B254: an unregistered surface name indexes the query result "
               "with [0] and no emptiness check, so it raises IndexError('list "
               "index out of range') instead of an error naming the surface that "
               "could not be found. See the consolidated findings issue.",
    )
    def test_an_unknown_surface_name_is_reported_clearly(self, unit_project):
        with pytest.raises(ValueError, match="NoSuchSurface"):
            depositionModels(
                projectName=unit_project.projectName, surface="NoSuchSurface"
            )

    def test_an_unknown_surface_name_currently_raises_indexerror(self, unit_project):
        """Characterisation of B254."""
        with pytest.raises(IndexError, match="list index out of range"):
            depositionModels(
                projectName=unit_project.projectName, surface="NoSuchSurface"
            )

    def test_a_surface_registered_under_another_type_is_not_found(self, unit_project):
        """Characterisation of B254: the query is narrowed by type, so a
        near-miss fails the same silent way."""
        unit_project.addCacheDocument(
            resource="", dataFormat="string", type="NotASurface",
            desc=dict(surface="Sand", zrough=0.002),
        )
        with pytest.raises(IndexError):
            depositionModels(projectName=unit_project.projectName, surface="Sand")


@pytest.mark.unit
class TestGroundDepositionEfficiency:
    """B255: see the module docstring."""

    def test_the_bare_ground_rate_is_positive_and_finite(self, make):
        rate = make(diameter=1e-6).depositionRate()
        assert rate > 0
        assert rate == rate

    def test_with_no_canopy_the_rate_is_exactly_ustar_times_the_efficiency(
        self, make
    ):
        """hcanop defaults to 0, so lm is 0, the `if lm != 0 else 0` guard sets
        Qsol to 0 and vds collapses to u* egb."""
        model = make(diameter=1e-6)
        efficiency = model.depositionRate() / model.ustar
        assert efficiency == pytest.approx(85142.1, rel=1e-4)

    @pytest.mark.xfail(
        strict=True,
        reason="B255: depositionRate_Petroff computes dpm = 0.000001*dp, a "
               "micron-to-metre conversion, but dp is self.diameter, which the "
               "constructor already normalised to metres. The surrounding "
               "constants are SI (mean free path 6.7e-8 m, dynamic viscosity "
               "1.89e-5 Pa s), so a 1 micron particle is evaluated as a 1 "
               "picometre one, deep in the free-molecule regime. The Brownian "
               "ground efficiency vds/u* comes out at 8.5e4; being an efficiency "
               "it cannot exceed 1, and for a 1 micron particle it should be of "
               "order 1e-5. See the consolidated findings issue.",
    )
    def test_a_deposition_velocity_cannot_exceed_the_friction_velocity(self, make):
        """vds = u* egb/(1 + Qsol) with egb a dimensionless ground transfer
        efficiency, so vds <= u* for every particle."""
        model = make(diameter=1e-6)
        assert model.depositionRate() <= model.ustar

    def test_the_efficiency_currently_exceeds_one_by_five_orders(self, make):
        """Characterisation of B255."""
        model = make(diameter=1e-6)
        assert model.depositionRate() / model.ustar > 1e4

    def test_the_conversion_is_applied_to_an_already_converted_diameter(self, make):
        """Characterisation of B255's mechanism, read off both sites."""
        import inspect

        constructor = inspect.getsource(depositionModels.__init__)
        assert "self._diameter = tonumber(tounit(diameter, ureg.m), ureg.m)" in (
            constructor
        )
        rate = inspect.getsource(depositionModels.depositionRate_Petroff)
        assert "dpm = 0.000001*dp" in rate
        assert make(diameter=1 * ureg.micrometer)._diameter == pytest.approx(1e-6)

    def test_brownian_deposition_weakens_as_the_particle_grows(self, make):
        """Brownian diffusivity falls with size, so the diffusive ground sink
        must weaken monotonically -- true whatever units dpm is in."""
        efficiencies = []
        for diameter in (1e-7, 1e-6, 1e-5):
            model = make(diameter=diameter)
            efficiencies.append(model.depositionRate() / model.ustar)
        assert efficiencies == sorted(efficiencies, reverse=True)


@pytest.mark.unit
class TestBoltzmannConstant:
    """B257: see the module docstring."""

    @staticmethod
    def _constantInSource():
        """Returned in units of 1e-23 J/K: pytest.approx's default absolute
        tolerance of 1e-12 swamps any comparison at the raw magnitude."""
        import inspect

        source = inspect.getsource(depositionModels.depositionRate_Petroff)
        match = re.search(r"kb\s*=\s*([0-9.]+[Ee][-+]?[0-9]+)", source)
        assert match is not None, "the kb assignment moved"
        return float(match.group(1)) * 1e23

    @pytest.mark.xfail(
        strict=True,
        reason="B257: the Brownian diffusivity uses kb = 1.83E-23 where the "
               "Boltzmann constant is 1.380649e-23 J/K -- the first two digits "
               "are transposed. Db = cu kb T / (3 pi muA dpm) is therefore 33% "
               "too large, and the Schmidt number 33% too small, everywhere the "
               "Petroff model is used. See the consolidated findings issue.",
    )
    def test_it_is_the_boltzmann_constant(self):
        assert self._constantInSource() == pytest.approx(1.380649, rel=1e-3)

    def test_it_is_currently_a_transposition_of_it(self):
        """Characterisation of B257, with the size of the error."""
        used = self._constantInSource()
        assert used == pytest.approx(1.83)
        assert used / 1.380649 == pytest.approx(1.3255, rel=1e-3)

    def test_the_other_constants_in_the_block_are_right(self):
        """So B257 is a typo rather than a different unit system: every
        other constant is its standard SI value for air at 293 K."""
        import inspect

        source = inspect.getsource(depositionModels.depositionRate_Petroff)
        for assignment in (
            "nuA = 0.0000157",     # kinematic viscosity, m2/s
            "muA = 0.0000189",     # dynamic viscosity, Pa s
            "g = 9.81",            # m/s2
            "lpm = 0.000000067",   # mean free path, m
            "rhoA = 1.2",          # kg/m3
            "cpA = 1000",          # J/(kg K)
        ):
            assert assignment in source


@pytest.mark.unit
class TestInCanopyResistanceWithoutAnObstacleShape:
    """The `obstacleShape is None` branch with a non-zero canopy height, the
    only path on which PhiH's return value matters."""

    def test_a_canopy_height_adds_resistance(self, make):
        """With hcanop > 0, Qsol = hcanop/lm * egb > 0 and
        vds = u* egb/(1 + Qsol) must fall below the bare-ground u* egb."""
        bare = make(surface=DESERT).depositionRate()
        canopy = make(surface=GRASSLAND).depositionRate()
        assert 0 < canopy < bare

    @staticmethod
    def _qsol(make, heatFlux, surface=GRASSLAND):
        """Recover Qsol from the two runs: vds = u* egb/(1 + Qsol), and the
        bare-ground run at the same particle size measures egb (which depends
        on the particle and the air, not on the surface)."""
        bare = make(surface=DESERT)
        canopy = make(surface=surface, heatFlux=heatFlux)
        egb = bare.depositionRate() / bare.ustar
        return (egb * canopy.ustar) / canopy.depositionRate() - 1

    @staticmethod
    def _expectedQsol(make, heatFlux, hcanop=10.0, hdepl=6.0):
        """Qsol = hcanop/lm * egb with lm = kappa (hcanop - hdepl)/PhiH(xi),
        xi = (hcanop - hdepl)/L and L = -u*^3/kappa * T/g * rho cp / H --
        the documented Monin-Obukhov form, evaluated here rather than read
        from the code (kappa = 0.4, T = 293 K, rho = 1.2, cp = 1000)."""
        bare = make(surface=DESERT)
        egb = bare.depositionRate() / bare.ustar
        kappa = 0.4
        obukhov = (
            -(bare.ustar**3) / kappa * 293.0 / 9.81 * 1.2 * 1000 / heatFlux
        )
        xi = (hcanop - hdepl) / obukhov
        phiH = (1 - 16 * xi) ** -0.5 if xi < 0 else 1 + 5 * xi
        lm = kappa * (hcanop - hdepl) / phiH
        return hcanop / lm * egb

    def test_the_resistance_matches_the_documented_form(self, make):
        """Qsol recovered from the deposition velocities must equal the
        Monin-Obukhov expression it is defined by."""
        assert self._qsol(make, 0.1) == pytest.approx(
            self._expectedQsol(make, 0.1), rel=1e-6
        )

    def test_the_resistance_matches_it_under_strong_convection_too(self, make):
        """Same check where the stability correction is no longer a small
        perturbation: PhiH is 0.524 at H = 100 W/m2 rather than 0.999."""
        assert self._qsol(make, 100.0) == pytest.approx(
            self._expectedQsol(make, 100.0), rel=1e-6
        )

    def test_the_stability_correction_scales_the_resistance(self, make):
        """Qsol is proportional to PhiH, so the stable-to-unstable ratio at
        +-100 W/m2 must be (1 + 5 xi)/(1 - 16 xi)^-0.5 = 1.827/0.524 = 3.487
        with xi = 0.1653."""
        ratio = self._qsol(make, -100.0) / self._qsol(make, 100.0)
        assert ratio == pytest.approx(3.487, rel=1e-3)

    def test_a_stable_surface_uses_the_linear_stability_function(self, make):
        """A downward heat flux makes the Obukhov length positive, so
        xi = (hcanop - hdepl)/L > 0 and PhiH takes its 1 + 5 xi branch."""
        stable = make(surface=GRASSLAND, heatFlux=-100.0)
        assert stable.depositionRate() > 0

    def test_a_convective_surface_deposits_faster_than_a_stable_one(self, make):
        """PhiH is below 1 when the flux is upward and above 1 when it is
        downward, so the mixing length is longer and the in-canopy resistance
        Qsol smaller under convection -- which is the whole point of carrying
        a stability correction."""
        unstable = make(surface=GRASSLAND, heatFlux=100.0).depositionRate()
        stable = make(surface=GRASSLAND, heatFlux=-100.0).depositionRate()
        assert stable < unstable

    def test_a_strongly_convective_surface_is_out_of_range(self, make):
        """PhiH accepts only -2 <= xi <= 1; a large upward heat flux shortens
        the Obukhov length until xi drops below -2."""
        with pytest.raises(KeyError, match="PhiH is out of stability range"):
            make(surface=GRASSLAND, heatFlux=5000.0).depositionRate()

    def test_a_strongly_stable_surface_is_out_of_range_too(self, make):
        """And a large downward flux pushes xi above 1."""
        with pytest.raises(KeyError, match="PhiH is out of stability range"):
            make(surface=GRASSLAND, heatFlux=-5000.0).depositionRate()

    def test_a_canopy_with_no_height_skips_the_resistance_entirely(self, make):
        """hcanop absent means lm == 0, and the `if lm != 0 else 0` guard is
        what keeps that from dividing by zero."""
        model = make(surface={"name": "Flat", "zrough": 0.01})
        assert model.depositionRate() == pytest.approx(
            model.ustar * (model.depositionRate() / model.ustar)
        )
        assert model.depositionRate() > 0
