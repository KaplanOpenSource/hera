"""gaussianToolkit: the sigma registry and the meteorology pass-throughs.

The toolkit inherits from Project, so it needs a database to construct -- the
mongomock seam supplies one.  The behaviour under test is pure dispatch: which
sigma class a name resolves to, and that the meteorology helpers agree with the
factory they delegate to.
"""
import pytest

from hera import toolkitHome
from hera.simulations.gaussian.Sigma import AbstractSigma, BriggsRural
from hera.utils.unitHandler import ureg


@pytest.fixture()
def toolkit(unit_toolkit_factory):
    return unit_toolkit_factory(toolkitHome.GAUSSIANDISPERSION)


@pytest.mark.unit
class TestSigmaRegistry:
    def test_briggs_rural_is_registered(self, toolkit):
        assert "briggsRural" in toolkit.listSigmaTypes()

    def test_the_name_resolves_to_an_instance_not_the_class(self, toolkit):
        """getSigmaType calls the class, so callers get a usable object."""
        resolved = toolkit.getSigmaType("briggsRural")
        assert isinstance(resolved, BriggsRural)
        assert resolved is not BriggsRural

    def test_the_resolved_object_is_a_sigma_model(self, toolkit):
        assert isinstance(toolkit.getSigmaType("briggsRural"), AbstractSigma)

    def test_the_resolved_object_computes(self, toolkit):
        """Resolution is only useful if the result actually works."""
        sigma = toolkit.getSigmaType("briggsRural")
        published = 0.08 * 1000.0 * (1 + 1e-4 * 1000.0) ** -0.5
        assert sigma.getSigma(1000.0, "D", units=False)["sigmaY"][0] == pytest.approx(
            published
        )

    def test_two_resolutions_give_independent_instances(self, toolkit):
        assert toolkit.getSigmaType("briggsRural") is not toolkit.getSigmaType(
            "briggsRural"
        )

    def test_an_unknown_name_lists_the_valid_ones(self, toolkit):
        """The error has to tell the caller what they could have written."""
        with pytest.raises(ValueError) as raised:
            toolkit.getSigmaType("pasquillGifford")

        message = str(raised.value)
        assert "pasquillGifford" in message
        assert "briggsRural" in message

    def test_the_error_is_a_valueerror_not_a_keyerror(self, toolkit):
        """The KeyError is translated deliberately; a bare KeyError would leak
        the internal dict as the API surface."""
        with pytest.raises(ValueError):
            toolkit.getSigmaType("nope")

    def test_listing_matches_what_resolves(self, toolkit):
        for name in toolkit.listSigmaTypes():
            assert toolkit.getSigmaType(name) is not None


@pytest.mark.unit
class TestMeteorologyPassThrough:
    def test_from_u10_agrees_with_the_factory(self, toolkit):
        from hera.simulations.gaussian.Meteorology import MeteorologyFactory

        through_toolkit = toolkit.getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s, inversion=1000 * ureg.m
        )
        through_factory = MeteorologyFactory().getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s,
            inversion=1000 * ureg.m,
            verticalProfileType="powerLaw",
        )
        assert type(through_toolkit) is type(through_factory)
        assert through_toolkit.getWindVelocity(50 * ureg.m) == through_factory.getWindVelocity(
            50 * ureg.m
        )

    def test_the_toolkit_defaults_to_the_power_law(self, toolkit):
        """Note the difference from MeteorologyFactory, which defaults to log."""
        from hera.simulations.gaussian.Meteorology import (
            StandardMeteorolgyConstant_powerLaw,
        )

        built = toolkit.getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s, inversion=1000 * ureg.m
        )
        assert type(built) is StandardMeteorolgyConstant_powerLaw

    def test_the_two_defaults_differ_between_toolkit_and_factory(self, toolkit):
        """Pinned because it is a trap: the same call through two entry points
        gives two different wind profiles."""
        from hera.simulations.gaussian.Meteorology import (
            MeteorologyFactory,
            StandardMeteorolgyConstant_log,
            StandardMeteorolgyConstant_powerLaw,
        )

        viaToolkit = toolkit.getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s, inversion=1000 * ureg.m
        )
        viaFactory = MeteorologyFactory().getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s, inversion=1000 * ureg.m
        )
        assert type(viaToolkit) is StandardMeteorolgyConstant_powerLaw
        assert type(viaFactory) is StandardMeteorolgyConstant_log

    def test_from_reference_height_reproduces_the_reference_wind(self, toolkit):
        built = toolkit.getMeteorologyFromURefHeight(
            u=8 * ureg.m / ureg.s, refHeight=50 * ureg.m, inversion=1000 * ureg.m
        )
        assert built.getWindVelocity(50 * ureg.m).to(
            ureg.m / ureg.s
        ).magnitude == pytest.approx(8.0)

    @pytest.mark.parametrize("profileType", ["powerLaw", "log", "uniformWind"])
    def test_every_profile_type_is_reachable(self, toolkit, profileType):
        built = toolkit.getMeteorologyFromU10(
            u10=5 * ureg.m / ureg.s,
            inversion=1000 * ureg.m,
            verticalProfileType=profileType,
        )
        assert built.getWindVelocity(20 * ureg.m) is not None


@pytest.mark.unit
class TestToolkitStructure:
    def test_the_toolkit_name_is_the_documented_one(self, toolkit):
        assert toolkit.toolkitName == "gaussianToolkit"

    def test_a_presentation_layer_is_attached(self, toolkit):
        assert toolkit.presentation is not None

    def test_presentation_is_composition_not_inheritance(self, toolkit):
        """CLAUDE.md: analysis and presentation must not subclass the toolkit."""
        from hera.toolkit import abstractToolkit

        assert not isinstance(toolkit.presentation, abstractToolkit)

    def test_the_toolkit_is_a_project(self, toolkit):
        from hera.datalayer import Project

        assert isinstance(toolkit, Project)

    def test_datasource_management_is_inherited(self, toolkit):
        """It comes from abstractToolkit, so it must work here unchanged."""
        toolkit.addDataSource("GAUSS_SRC", "/data/case", "string", version=(0, 0, 1))
        assert toolkit.getDataSourceDocument("GAUSS_SRC").resource == "/data/case"

    def test_its_datasources_are_tagged_with_its_own_name(self, toolkit):
        toolkit.addDataSource("GAUSS_SRC", "/data/case", "string")
        assert (
            toolkit.getDataSourceDocument("GAUSS_SRC").desc["toolkit"]
            == "gaussianToolkit"
        )
