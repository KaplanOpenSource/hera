"""Dry-deposition models.

Three of this class's properties are tangled together by two copy-paste
errors, and one of them corrupts the friction velocity that the deposition
rate is computed from.  The tests below establish what each property does
today and assert what it should do.
"""
import pytest

from hera.simulations.deposition.models import depositionModels
from hera.utils.unitHandler import ureg


@pytest.fixture()
def model(unit_project):
    """depositionModels constructs a Project, so it needs the in-memory store."""
    return depositionModels(projectName=unit_project.projectName)


@pytest.mark.unit
class TestConstruction:
    def test_the_documented_defaults_are_applied(self, model):
        assert model.temperature == pytest.approx(293.0)
        assert model.density == 1500
        assert model.surface["name"] == "Desert"
        assert model.surface["zrough"] == pytest.approx(0.04)

    def test_the_default_diameter_is_one_micron(self, model):
        assert model._diameter == pytest.approx(1e-6)

    def test_the_default_heat_flux(self, model):
        assert model._heatFlux == pytest.approx(0.1)

    def test_temperature_is_converted_to_kelvin(self, unit_project):
        """Offset units need Quantity(); `20 * ureg.degC` is ambiguous in pint."""
        built = depositionModels(
            projectName=unit_project.projectName,
            temperature=ureg.Quantity(20, ureg.degC),
        )
        assert built.temperature == pytest.approx(293.15)

    def test_a_diameter_in_microns_is_converted_to_metres(self, unit_project):
        built = depositionModels(
            projectName=unit_project.projectName, diameter=5 * ureg.micrometer
        )
        assert built._diameter == pytest.approx(5e-6)

    def test_a_surface_dict_is_taken_as_given(self, unit_project):
        built = depositionModels(
            projectName=unit_project.projectName,
            surface={"name": "Grass", "zrough": 0.01},
        )
        assert built.surface["zrough"] == pytest.approx(0.01)

    def test_the_model_name_is_recorded(self, model):
        assert model._depositionModel == "Petroff"


@pytest.mark.unit
class TestTangledProperties:
    r"""Two copy-paste errors leave three properties reading each other.

        @heatFlux.setter
        def heatFlux(self, newheatFlux):
            self._ustar = newheatFlux     # writes ustar, not heatFlux

        @diameter.setter
        def ustar(self, newdiameter):     # named ustar, so it REPLACES the
            self._diameter = newdiameter  # ustar property defined above
    """

    def test_ustar_currently_reads_the_diameter(self, model):
        """Characterisation: the diameter setter's name replaced the property."""
        model._ustar = 111
        model._diameter = 999
        assert model.ustar == 999

    def test_diameter_has_no_setter_at_all(self, model):
        """Because @diameter.setter produced a property named ustar instead."""
        with pytest.raises(AttributeError):
            model.diameter = 7

    def test_assigning_heat_flux_currently_writes_ustar(self, model):
        model._ustar = 111
        model.heatFlux = 42
        assert model._ustar == 42

    @pytest.mark.xfail(
        strict=True,
        reason="B46: the diameter setter is named `ustar`, so @diameter.setter "
               "builds a new property under that name and replaces the real ustar "
               "property. Reading obj.ustar returns _diameter. ustar is the "
               "central parameter of the deposition calculation. "
               "See the consolidated findings issue.",
    )
    def test_ustar_reads_the_friction_velocity(self, model):
        model._ustar = 111
        model._diameter = 999
        assert model.ustar == 111

    @pytest.mark.xfail(
        strict=True,
        reason="B46: assigning ustar writes _diameter. "
               "See the consolidated findings issue.",
    )
    def test_assigning_ustar_writes_the_friction_velocity(self, model):
        model.ustar = 5
        assert model._ustar == 5

    @pytest.mark.xfail(
        strict=True,
        reason="B47: the heatFlux setter assigns self._ustar, so heatFlux never "
               "changes and the friction velocity is silently overwritten. "
               "See the consolidated findings issue.",
    )
    def test_assigning_heat_flux_writes_the_heat_flux(self, model):
        model.heatFlux = 42
        assert model._heatFlux == 42

    @pytest.mark.xfail(
        strict=True,
        reason="B47: and it leaves ustar alone. "
               "See the consolidated findings issue.",
    )
    def test_assigning_heat_flux_leaves_ustar_alone(self, model):
        model._ustar = 111
        model.heatFlux = 42
        assert model._ustar == 111

    def test_the_untangled_properties_behave_normally(self, model):
        """density, temperature and surface are unaffected, for contrast."""
        model.density = 2000
        model.temperature = 300
        model.surface = {"name": "X", "zrough": 0.5}
        assert model.density == 2000
        assert model.temperature == 300
        assert model.surface["name"] == "X"


@pytest.mark.unit
class TestDepositionRateDispatch:
    def test_the_model_name_selects_the_method(self, model):
        """depositionRate defers to depositionRate_<name>."""
        assert hasattr(model, "depositionRate_Petroff")

    def test_an_unknown_model_name_raises(self, unit_project):
        built = depositionModels(
            projectName=unit_project.projectName, depositionModel="NoSuchModel"
        )
        with pytest.raises(AttributeError):
            built.depositionRate()

    def test_petroff_returns_a_finite_positive_rate(self, model):
        """A deposition velocity is positive and finite for ordinary inputs."""
        rate = model.depositionRate()
        assert rate > 0
        assert rate == rate  # not NaN

    @staticmethod
    def _rateAt(project, effectiveUstar=0.3, **kwargs):
        """Rate at a realistic friction velocity, working around B46.

        obj.ustar reads _diameter, so the only way to feed the calculation a
        plausible u* is to write that field.  Doing it here keeps the
        sensitivity tests below about roughness and density rather than about
        B46.
        """
        built = depositionModels(projectName=project.projectName, **kwargs)
        built._diameter = effectiveUstar
        return built.depositionRate()

    @pytest.mark.xfail(
        strict=True,
        reason="B48: the deposition rate is completely insensitive to surface "
               "roughness. z0 = 0.001 and z0 = 0.5 give a bit-identical result "
               "(0.0000079953) at a realistic u*, so the `surface` parameter -- "
               "the reason the class takes one -- has no effect on the answer. "
               "See the consolidated findings issue.",
    )
    def test_the_rate_responds_to_the_surface_roughness(self, unit_project):
        smooth = self._rateAt(unit_project, surface={"name": "S", "zrough": 0.001})
        rough = self._rateAt(unit_project, surface={"name": "R", "zrough": 0.5})
        assert smooth != rough

    @pytest.mark.xfail(
        strict=True,
        reason="B49: the deposition rate is completely insensitive to particle "
               "density. 500 and 5000 kg/m3 give a bit-identical result, though "
               "gravitational settling scales with density. "
               "See the consolidated findings issue.",
    )
    def test_the_rate_responds_to_particle_density(self, unit_project):
        light = self._rateAt(unit_project, density=500)
        heavy = self._rateAt(unit_project, density=5000)
        assert heavy != light

    def test_roughness_and_density_are_bit_identical_today(self, unit_project):
        """Characterisation, so B48 and B49 are numbers rather than adjectives."""
        smooth = self._rateAt(unit_project, surface={"name": "S", "zrough": 0.001})
        rough = self._rateAt(unit_project, surface={"name": "R", "zrough": 0.5})
        light = self._rateAt(unit_project, density=500)
        heavy = self._rateAt(unit_project, density=5000)
        assert smooth == rough == light == heavy

    def test_the_friction_velocity_path_itself_works(self, unit_project):
        """Feeding u* through the backing field DOES change the rate.

        Which localises B46: the calculation uses u* correctly, it is only
        handed the wrong value by the broken property.
        """
        calm = self._rateAt(unit_project, effectiveUstar=0.1)
        windy = self._rateAt(unit_project, effectiveUstar=2.0)
        assert windy > calm

    @pytest.mark.xfail(
        strict=True,
        reason="B46 in effect: depositionRate_Petroff reads `ustar = self.ustar`, "
               "which returns the particle DIAMETER. Passing a real u* to the "
               "constructor therefore changes nothing -- a length is substituted "
               "for a velocity. See the consolidated findings issue.",
    )
    def test_the_constructor_ustar_reaches_the_calculation(self, unit_project):
        calm = depositionModels(
            projectName=unit_project.projectName, ustar=0.1
        ).depositionRate()
        windy = depositionModels(
            projectName=unit_project.projectName, ustar=2.0
        ).depositionRate()
        assert calm != windy
