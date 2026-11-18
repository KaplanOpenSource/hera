"""Evaporation models: the agent-facing members.

``test_evaporation_models.py`` covers the two diffusion correlations, the
air-viscosity power law, the dimensionless groups, the flux dispatch and
the constructor (B50/B51).  What it never touches is the agent plumbing
those correlations read their molecular data from:

* ``evaporationModels.agent`` (property and setter)
* ``evaporationModels.Magent``  -- molecular weight, g/mol
* ``evaporationModels.Vagent``  -- molar volume, cm3/mol

The assertions are the published forms of the two correlations, which fix
exactly how the two agent numbers must enter:

    FSG   D ~ sqrt(1/Mair + 1/Magent) / (Vair^1/3 + Vagent^1/3)^2
    EPA   D ~ sqrt(1/Mair + 1/Magent) / Magent^1/3

so a second agent with a different molecular weight pins the dependence
without any reference to the implementation.

One defect surfaced:

* B256, fixed: the ``agent`` setter called ``RiskToolkit.getAgent(newAgent)``
  on the class, but ``getAgent`` is an instance method
  ``(self, nameOrDesc, version=None)``, so the new agent was bound to
  ``self``, ``nameOrDesc`` was left missing and every assignment raised
  TypeError.  Same shape as B50 in ``__init__``, a different code site.
  It now goes through the toolkit the constructor builds -- so an object
  assembled by hand, as the helper below does, has to be given one.  The
  class also declared ``agent`` twice; the first, setter-less getter was
  dead code and is gone.
"""
import numpy as np
import pytest

from hera.simulations.evaporation.models import evaporationModels
from hera.utils.unitHandler import ureg


def _agentDescriptor(molecularWeight="100 g/mol"):
    """The minimum an Agent accepts; note the capital P in physicalProperties."""
    return {
        "name": "test_agent",
        "effectParameters": {"tenbergeCoefficient": 1},
        "physicalProperties": {
            "molecularWeight": molecularWeight,
            "molecularVolume": 80.0,
        },
        "effects": {
            "inhalation": {
                "type": "Lognormal10",
                "calculator": {"TenBerge": {}},
                "units": "mg/m**3",
                "parameters": {
                    "type": "Lognormal10DoseResponse",
                    "levels": ["Severe"],
                    "parameters": {"Severe": {"TL_50": 400, "sigma": 0.5}},
                },
            }
        },
    }


def _build(molecularWeight="100 g/mol", molarVolume=80.0, diffusionModel="FSG"):
    """Assemble the object field by field, which keeps these tests off the
    database and off the constructor's molecularVolume fallback."""
    from hera.riskassessment.agents.Agents import Agent

    built = evaporationModels.__new__(evaporationModels)
    built._agent = Agent(_agentDescriptor(molecularWeight))
    built._Mair = 28.967
    built._Vair = 20.1
    built._Vagent = molarVolume
    built._molecularDiffusionModel = diffusionModel
    built._dinamicViscocityModel = "powerLaw"
    built._evaporationModel = "US"
    return built


@pytest.fixture()
def models():
    return _build()


@pytest.mark.unit
class TestAgentProperty:
    def test_the_agent_is_exposed(self, models):
        assert models.agent is models._agent

    def test_the_agent_carries_its_physical_properties(self, models):
        assert models.agent.physicalproperties is not None

    def test_the_agent_property_has_a_getter_and_a_setter(self):
        agentProperty = evaporationModels.__dict__["agent"]
        assert agentProperty.fget is not None
        assert agentProperty.fset is not None

    def test_the_agent_property_is_declared_once(self):
        """B256: the class body used to define `agent` twice -- a bare getter
        near the top that the later getter/setter pair shadowed."""
        import inspect

        source = inspect.getsource(evaporationModels)
        assert source.count("def agent(self):") == 1

    def test_the_agent_can_be_reassigned(self, models, unit_project):
        """A hand-assembled object has no toolkit, so it is given one."""
        from hera.riskassessment import RiskToolkit

        models._riskToolkit = RiskToolkit(projectName=unit_project.projectName)
        models.agent = _agentDescriptor("50 g/mol")
        assert models.Magent == pytest.approx(50.0)

    def test_the_setter_goes_through_the_toolkit_instance(self):
        """B256's mechanism: getAgent needs an instance, not the class."""
        import inspect

        from hera.riskassessment import RiskToolkit

        parameters = list(inspect.signature(RiskToolkit.getAgent).parameters)
        assert parameters[:2] == ["self", "nameOrDesc"]
        setterSource = inspect.getsource(evaporationModels)
        assert "self._agent = self._riskToolkit.getAgent(newAgent)" in setterSource


@pytest.mark.unit
class TestMolecularWeight:
    def test_it_is_reported_in_grams_per_mole(self, models):
        """The descriptor says 100 g/mol, and Magent is documented as g/mol."""
        assert models.Magent == pytest.approx(100.0)

    def test_it_is_a_bare_number(self, models):
        """m_as strips the units, so the correlations can use it directly."""
        assert not hasattr(models.Magent, "units")

    def test_a_weight_given_in_kilograms_per_mole_is_converted(self):
        """0.25 kg/mol is 250 g/mol."""
        assert _build(molecularWeight="0.25 kg/mol").Magent == pytest.approx(250.0)

    def test_a_heavier_agent_diffuses_more_slowly_under_epa(self):
        """EPA: D ~ sqrt(1/Mair + 1/Mb) / Mb^(1/3), decreasing in Mb."""
        light = _build(molecularWeight="50 g/mol", diffusionModel="EPA")
        heavy = _build(molecularWeight="200 g/mol", diffusionModel="EPA")
        assert heavy.molecularDiffusion(300.0 * ureg.K) < light.molecularDiffusion(
            300.0 * ureg.K
        )

    def test_the_epa_ratio_follows_the_published_dependence(self):
        """The temperature factor cancels in a ratio, leaving only the
        molecular-weight dependence to check."""
        light = _build(molecularWeight="50 g/mol", diffusionModel="EPA")
        heavy = _build(molecularWeight="200 g/mol", diffusionModel="EPA")
        expected = (
            np.sqrt(1 / 28.967 + 1 / 200.0)
            / np.cbrt(200.0)
            / (np.sqrt(1 / 28.967 + 1 / 50.0) / np.cbrt(50.0))
        )
        ratio = heavy.molecularDiffusion_EPA(300.0) / light.molecularDiffusion_EPA(300.0)
        assert ratio == pytest.approx(expected)

    def test_a_heavier_agent_diffuses_more_slowly_under_fsg_too(self):
        """FSG: only through sqrt(1/Mair + 1/Mb), so the effect is weaker but
        has the same sign."""
        light = _build(molecularWeight="50 g/mol")
        heavy = _build(molecularWeight="200 g/mol")
        assert heavy.molecularDiffusion_FSG(300.0) < light.molecularDiffusion_FSG(300.0)

    def test_the_weight_dependence_is_weaker_under_fsg_than_under_epa(self):
        """FSG has the molar-volume term carrying the size dependence, so its
        mass ratio must be the milder of the two."""
        light, heavy = _build("50 g/mol"), _build("200 g/mol")
        fsgRatio = heavy.molecularDiffusion_FSG(300.0) / light.molecularDiffusion_FSG(
            300.0
        )
        epaRatio = heavy.molecularDiffusion_EPA(300.0) / light.molecularDiffusion_EPA(
            300.0
        )
        assert epaRatio < fsgRatio < 1


@pytest.mark.unit
class TestMolarVolume:
    def test_it_reports_the_stored_molar_volume(self, models):
        assert models.Vagent == pytest.approx(80.0)

    def test_a_bulkier_molecule_diffuses_more_slowly(self):
        """FSG: D ~ 1/(Vair^1/3 + Vb^1/3)^2, decreasing in Vb."""
        small = _build(molarVolume=20.0)
        large = _build(molarVolume=200.0)
        assert large.molecularDiffusion_FSG(300.0) < small.molecularDiffusion_FSG(300.0)

    def test_the_ratio_follows_the_inverse_square_of_the_summed_cube_roots(self):
        small, large = _build(molarVolume=20.0), _build(molarVolume=200.0)
        expected = (
            (np.cbrt(20.1) + np.cbrt(20.0)) / (np.cbrt(20.1) + np.cbrt(200.0))
        ) ** 2
        ratio = large.molecularDiffusion_FSG(300.0) / small.molecularDiffusion_FSG(300.0)
        assert ratio == pytest.approx(expected)

    def test_the_molar_volume_does_not_enter_the_epa_form(self):
        """EPA exists precisely for agents whose molar volume is unknown."""
        small = _build(molarVolume=20.0, diffusionModel="EPA")
        large = _build(molarVolume=200.0, diffusionModel="EPA")
        assert small.molecularDiffusion(300.0 * ureg.K) == pytest.approx(
            large.molecularDiffusion(300.0 * ureg.K)
        )

    def test_a_missing_molar_volume_breaks_the_fsg_form(self):
        """Which is why the constructor is meant to fall back to EPA: with
        Vagent None the cube root raises rather than returning nonsense."""
        broken = _build(molarVolume=None)
        with pytest.raises(TypeError):
            broken.molecularDiffusion_FSG(300.0)


@pytest.mark.unit
class TestDinamicViscocityModelProperty:
    """The last uncovered member of this class: a plain round-tripping
    property, spelled `dinamicViscocityModel` (sic -- the misspelling is
    part of the public surface)."""

    def test_it_round_trips(self, models):
        models.dinamicViscocityModel = "somethingElse"
        assert models.dinamicViscocityModel == "somethingElse"

    def test_the_constructor_value_is_visible_through_the_getter(self, models):
        assert models.dinamicViscocityModel == models._dinamicViscocityModel
