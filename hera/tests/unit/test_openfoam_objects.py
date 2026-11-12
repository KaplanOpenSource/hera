"""OpenFOAM preprocessing objects.

The heaviest area in hera and the least tested -- 257 public callables, of
which 2 were exercised before this file.  It imports PyFoam and evtk at module
level, so everything here runs on the Phase 0 stub layer.

The parts covered are the ones that are genuine string and dictionary
manipulation: the dimension vector, the field-definition registry, and the
zero values for scalar, vector and tensor fields.  OpenFOAM's dimension vector
is [kg m s K mol A cd], which makes those assertions checkable against the
format rather than against current output.
"""
import pandas as pd
import pytest

from hera.simulations.openFoam.preprocessOFObjects.OFObjectHome import OFObjectHome


@pytest.fixture()
def home():
    return OFObjectHome()


@pytest.mark.unit
class TestDimensionVector:
    """OpenFOAM writes dimensions as [kg m s K mol A cd]."""

    def test_a_dimensionless_quantity_is_all_zeros(self, home):
        assert home.getDimensions() == "[0 0 0 0 0 0 0]"

    def test_length_sits_in_the_second_slot(self, home):
        assert home.getDimensions(m=1) == "[0 1 0 0 0 0 0]"

    def test_mass_sits_in_the_first_slot(self, home):
        assert home.getDimensions(kg=1) == "[1 0 0 0 0 0 0]"

    def test_time_sits_in_the_third_slot(self, home):
        assert home.getDimensions(s=1) == "[0 0 1 0 0 0 0]"

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            (dict(m=1, s=-1), "[0 1 -1 0 0 0 0]"),                  # velocity
            (dict(m=2, s=-2), "[0 2 -2 0 0 0 0]"),                  # TKE
            (dict(kg=1, m=-1, s=-2), "[1 -1 -2 0 0 0 0]"),          # pressure
            (dict(kg=1, m=-3), "[1 -3 0 0 0 0 0]"),                 # density
            (dict(K=1), "[0 0 0 1 0 0 0]"),                         # temperature
        ],
    )
    def test_physical_quantities_get_their_known_vectors(self, home, kwargs, expected):
        assert home.getDimensions(**kwargs) == expected

    def test_all_seven_slots_are_addressable(self, home):
        assert home.getDimensions(kg=1, m=2, s=3, K=4, mol=5, A=6, cd=7) == (
            "[1 2 3 4 5 6 7]"
        )

    def test_it_is_callable_on_the_class(self, home):
        """Declared @staticmethod, so both call forms must agree."""
        assert OFObjectHome.getDimensions(m=1) == home.getDimensions(m=1)


@pytest.mark.unit
class TestFieldDefinitionRegistry:
    def test_definitions_are_loaded_at_construction(self, home):
        assert isinstance(home.fieldDefinitions, dict)
        assert len(home.fieldDefinitions) > 0

    def test_the_common_fields_are_known(self, home):
        """U, p and T are the fields any OpenFOAM case has."""
        known = set(home.fieldDefinitions)
        assert {"U", "p"} <= known

    def test_each_definition_carries_a_type_and_dimensions(self, home):
        for name, definition in home.fieldDefinitions.items():
            assert "fieldType" in definition, f"{name} has no fieldType"
            assert "dimensions" in definition, f"{name} has no dimensions"

    def test_velocity_is_a_vector(self, home):
        assert home.fieldDefinitions["U"]["fieldType"].lower() == "vector"

    def test_pressure_is_a_scalar(self, home):
        assert home.fieldDefinitions["p"]["fieldType"].lower() == "scalar"

    def test_a_new_definition_can_be_added(self, home):
        home.addFieldDefinitions(
            fieldName="myTracer",
            dimensions=home.getDimensions(kg=1, m=-3),
            fieldType="scalar",
        )
        assert "myTracer" in home.fieldDefinitions
        assert home.fieldDefinitions["myTracer"]["fieldType"] == "scalar"

    def test_adding_an_existing_name_without_overwrite_is_refused(self, home):
        with pytest.raises(Exception):
            home.addFieldDefinitions(
                fieldName="U",
                dimensions=home.getDimensions(m=1),
                fieldType="scalar",
            )

    def test_overwrite_replaces_the_definition(self, home):
        home.addFieldDefinitions(
            fieldName="U",
            dimensions=home.getDimensions(m=1),
            fieldType="scalar",
            overwrite=True,
        )
        assert home.fieldDefinitions["U"]["fieldType"] == "scalar"

    def test_the_registry_is_per_instance(self):
        """Two homes must not share mutable state."""
        first, second = OFObjectHome(), OFObjectHome()
        first.addFieldDefinitions(
            fieldName="onlyInFirst",
            dimensions=first.getDimensions(),
            fieldType="scalar",
        )
        assert "onlyInFirst" not in second.fieldDefinitions


@pytest.mark.unit
class TestEmptyFieldConstruction:
    def test_a_scalar_field_zero_is_a_plain_zero(self, home):
        field = home.getEmptyField("p", flowType="incompressible")
        assert field.getZeroValue() == 0

    def test_a_vector_field_zero_comes_from_pyfoam_vector(self, home):
        """The value itself is a PyFoam Vector, which is stubbed here.

        So the assertion is about the dispatch -- a vector field must not get
        the scalar zero -- rather than about the three components, which only a
        real PyFoam install could show.
        """
        vectorZero = home.getEmptyField("U", flowType="incompressible").getZeroValue()
        scalarZero = home.getEmptyField("p", flowType="incompressible").getZeroValue()
        assert scalarZero == 0
        assert vectorZero != 0
        assert type(vectorZero).__name__ != "int"

    def test_the_field_remembers_its_name(self, home):
        assert home.getEmptyField("p", flowType="incompressible").name == "p"

    def test_the_field_carries_its_dimensions(self, home):
        field = home.getEmptyField("U", flowType="incompressible")
        assert field.dimensions is not None

    def test_an_unknown_field_name_is_refused(self, home):
        with pytest.raises(Exception):
            home.getEmptyField("noSuchFieldAnywhere", flowType="incompressible")

    def test_a_uniform_internal_value_can_be_set(self, home):
        field = home.getEmptyField("p", flowType="incompressible")
        field.setInternalUniformFieldValue(42)
        assert "42" in str(field.data) or field.data is not None

    def test_a_boundary_field_can_be_added(self, home):
        field = home.getEmptyField("p", flowType="incompressible")
        field.addBoundaryField("inlet", type="fixedValue", value="uniform 0")
        assert field.data is not None


@pytest.mark.unit
class TestPandasToFoamFormat:
    """B52: this method cannot be called successfully in any way.

        @staticmethod
        def pandasToFoamFormat(self, data):
            ...
            D = data if self.componentNames is None else data[self.componentNames]

    It is decorated @staticmethod but declared with `self`, so the ordinary
    call passes the DataFrame as `self` and leaves `data` missing.  Supplying
    an instance explicitly gets further and then fails too, because
    componentNames is an OFField attribute -- OFObjectHome does not have one.
    """

    FRAME = pd.DataFrame({"a": [1.0, 2.0, 3.0]})

    def test_the_decorator_and_the_signature_disagree(self):
        import inspect

        assert isinstance(
            inspect.getattr_static(OFObjectHome, "pandasToFoamFormat"), staticmethod
        )
        parameters = list(
            inspect.signature(OFObjectHome.pandasToFoamFormat).parameters
        )
        assert parameters[0] == "self"

    def test_the_attribute_it_needs_is_on_a_different_class(self, home):
        assert not hasattr(home, "componentNames")

    def test_calling_it_normally_raises(self, home):
        with pytest.raises(TypeError, match="missing 1 required positional argument"):
            home.pandasToFoamFormat(self.FRAME)

    def test_passing_an_instance_explicitly_also_raises(self, home):
        with pytest.raises(AttributeError, match="componentNames"):
            OFObjectHome.pandasToFoamFormat(home, self.FRAME)

    @pytest.mark.xfail(
        strict=True,
        reason="B52: pandasToFoamFormat is decorated @staticmethod but declared "
               "with `self`, and its body reads self.componentNames -- an OFField "
               "attribute that OFObjectHome does not have. Neither call form "
               "works, so the method has never run. "
               "See the consolidated findings issue.",
    )
    def test_it_produces_the_openfoam_list_format(self, home):
        """Documented output:  <count>\\n(\\n<values>\\n);\\n"""
        formatted = home.pandasToFoamFormat(self.FRAME)
        assert formatted.startswith("3\n(\n")
        assert formatted.rstrip().endswith(");")
