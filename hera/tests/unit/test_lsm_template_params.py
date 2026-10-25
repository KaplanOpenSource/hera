"""LSM template parameter preparation.

prepareParams merges a template's defaults with a caller's overrides, converts
every unit-bearing value into the units the template declares, and strips the
units so the model sees plain numbers.  Each of those steps is asserted here
against a template whose declared units make the expected numbers obvious.
"""
import pytest

from hera.simulations.LSM.template import LSMTemplate
from hera.utils.unitHandler import ureg


@pytest.mark.unit
class TestUnitConversion:
    def test_a_declared_unit_drives_the_conversion(self):
        """The template asks for metres; 1 km must arrive as 1000."""
        template = {"params": {}, "units": {"height": "m"}}
        prepared = LSMTemplate.prepareParams(template, {"height": 1 * ureg.km})
        assert prepared["height"] == pytest.approx(1000.0)

    def test_the_result_carries_no_units(self):
        """The model reads plain numbers, so units must be stripped."""
        template = {"params": {}, "units": {"height": "m"}}
        prepared = LSMTemplate.prepareParams(template, {"height": 1 * ureg.km})
        assert not hasattr(prepared["height"], "magnitude")

    def test_a_string_with_units_is_parsed(self):
        template = {"params": {}, "units": {"height": "m"}}
        prepared = LSMTemplate.prepareParams(template, {"height": "2 km"})
        assert prepared["height"] == pytest.approx(2000.0)

    def test_duration_is_special_cased_to_minutes(self):
        """A bare number for duration is read as minutes, not as base units."""
        template = {"params": {}, "units": {"duration": "min"}}
        prepared = LSMTemplate.prepareParams(template, {"duration": 30})
        assert prepared["duration"] == pytest.approx(30.0)

    def test_duration_is_exempt_from_standardisation(self):
        """Everything else goes to base units; duration stays in minutes.

        Without the exemption 30 minutes would come out as 1800 seconds and
        the model would run for thirty hours.
        """
        template = {"params": {}, "units": {"duration": "min", "height": "m"}}
        prepared = LSMTemplate.prepareParams(
            template, {"duration": 30, "height": 1 * ureg.km}
        )
        assert prepared["duration"] == pytest.approx(30.0)
        assert prepared["height"] == pytest.approx(1000.0)

    def test_a_unitless_value_where_units_are_declared_is_rejected(self):
        template = {"params": {}, "units": {"height": "m"}}
        with pytest.raises(ValueError, match="must use either pint or unum"):
            LSMTemplate.prepareParams(template, {"height": 5})

    def test_an_undeclared_parameter_passes_through(self):
        template = {"params": {}, "units": {}}
        prepared = LSMTemplate.prepareParams(template, {"label": "run1"})
        assert prepared["label"] == "run1"


@pytest.mark.unit
class TestDefaultsAndOverrides:
    def test_template_defaults_are_kept(self):
        template = {"params": {"a": 1}, "units": {}}
        assert LSMTemplate.prepareParams(template, {})["a"] == 1

    def test_the_caller_overrides_the_template(self):
        template = {"params": {"a": 1}, "units": {}}
        assert LSMTemplate.prepareParams(template, {"a": 99})["a"] == 99

    def test_both_sources_are_merged(self):
        template = {"params": {"a": 1}, "units": {}}
        prepared = LSMTemplate.prepareParams(template, {"b": 2})
        assert prepared == {"a": 1, "b": 2}

    def test_a_missing_template_is_tolerated(self):
        assert LSMTemplate.prepareParams(None, {"a": 1})["a"] == 1

    def test_a_template_without_units_is_tolerated(self):
        assert LSMTemplate.prepareParams({"params": {"a": 1}}, {})["a"] == 1


@pytest.mark.unit
class TestIntegerFields:
    @pytest.mark.parametrize("field", ["TopoXn", "TopoYn"])
    def test_a_float_grid_size_is_cast_to_int(self, field):
        """Grid counts index arrays, so a float would fail much later."""
        prepared = LSMTemplate.prepareParams({"params": {}, "units": {}}, {field: 10.7})
        assert prepared[field] == 10
        assert isinstance(prepared[field], int)

    @pytest.mark.parametrize("field", ["TopoXn", "TopoYn"])
    def test_an_integer_is_left_alone(self, field):
        prepared = LSMTemplate.prepareParams({"params": {}, "units": {}}, {field: 12})
        assert prepared[field] == 12

    def test_other_numeric_fields_are_not_cast(self):
        prepared = LSMTemplate.prepareParams(
            {"params": {}, "units": {}}, {"someOtherSize": 10.7}
        )
        assert prepared["someOtherSize"] == pytest.approx(10.7)


@pytest.mark.unit
class TestTheTemplateIsNotMutated:
    """A template is reused across simulations, so preparing must not touch it."""

    @pytest.mark.xfail(
        strict=True,
        reason="B45: even a single call writes the caller's overrides into "
               "template_desc['params']. See the consolidated findings issue.",
    )
    def test_the_declared_defaults_survive_one_call(self):
        template = {"params": {"a": 1}, "units": {}}
        LSMTemplate.prepareParams(template, {"b": 2})
        assert template["params"] == {"a": 1}

    def test_the_mutation_is_what_it_looks_like(self):
        """Characterisation of B45's mechanism, so the cause is unambiguous."""
        template = {"params": {"a": 1}, "units": {}}
        LSMTemplate.prepareParams(template, {"b": 2})
        assert template["params"] == {"a": 1, "b": 2}

    @pytest.mark.xfail(
        strict=True,
        reason="B45: prepareParams does params = template_desc.get('params', {}) "
               "and then params.update(...), which writes the caller's overrides "
               "into the template's own dict. A template reused for a second "
               "simulation carries the first run's parameters. "
               "See the consolidated findings issue.",
    )
    def test_a_reused_template_does_not_accumulate_earlier_runs(self):
        """The scenario, verified: one template, two simulations.

            prepareParams(template, {"b": 2})   # template gains b
            prepareParams(template, {"c": 3})   # second run also sees b
        """
        template = {"params": {"a": 1}, "units": {}}

        LSMTemplate.prepareParams(template, {"b": 2})
        second = LSMTemplate.prepareParams(template, {"c": 3})

        assert "b" not in second, "the second run inherited the first run's parameter"
        assert second == {"a": 1, "c": 3}
