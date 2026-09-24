"""JSON <-> configuration conversion, including the unit round trip.

The central contract, from the docstrings: ConfigurationToJSON turns a config
dict (which may carry units) into JSON-safe values, and JSONToConfiguration
reverses it.  Most assertions below are round trips, because that is the
property callers actually depend on.
"""
import json

import pytest

from hera.utils.jsonutils import (
    ConfigurationToJSON,
    JSONToConfiguration,
    loadJSON,
    setJSONPath,
    stripConfigurationUnits,
)
from hera.utils.unitHandler import ureg


@pytest.mark.unit
class TestConfigurationToJSON:
    def test_plain_values_pass_through(self):
        assert ConfigurationToJSON({"a": 1, "b": "x"}) == {"a": 1, "b": "x"}

    def test_quantity_becomes_a_parseable_string(self):
        assert ConfigurationToJSON({"v": 5 * ureg.m}) == {"v": "5 meter"}

    def test_standardize_converts_to_base_units(self):
        """MKS: 5 km must come out as 5000 m."""
        assert ConfigurationToJSON({"v": 5 * ureg.km}, standardize=True) == {
            "v": "5000.0 meter"
        }

    def test_nesting_is_preserved(self):
        assert ConfigurationToJSON({"a": {"b": 2 * ureg.s}}) == {"a": {"b": "2 second"}}

    def test_lists_are_converted_elementwise(self):
        assert ConfigurationToJSON({"l": [1 * ureg.m, 2]}) == {"l": ["1 meter", 2]}

    def test_empty_dict_stays_empty(self):
        assert ConfigurationToJSON({}) == {}

    def test_result_is_json_serialisable(self):
        """The point of the function -- the output has to survive json.dumps."""
        encoded = ConfigurationToJSON({"v": 5 * ureg.km, "n": 3, "s": "text"})
        assert json.loads(json.dumps(encoded)) == encoded


@pytest.mark.unit
class TestSplitUnits:
    """splitUnits stores magnitude and units in separate fields."""

    def test_produces_the_three_documented_fields(self):
        encoded = ConfigurationToJSON({"v": 5 * ureg.km}, splitUnits=True)["v"]
        assert set(encoded) == {"magnitude", "units", "value"}

    def test_magnitude_is_numeric_and_queryable(self):
        """The stated purpose: allow __lt / __gt queries on the magnitude."""
        encoded = ConfigurationToJSON({"v": 5 * ureg.km}, splitUnits=True)["v"]
        assert encoded["magnitude"] == pytest.approx(5)
        assert isinstance(encoded["magnitude"], (int, float))

    def test_standardize_puts_base_units_in_the_magnitude(self):
        encoded = ConfigurationToJSON(
            {"v": 5 * ureg.km}, splitUnits=True, standardize=True
        )["v"]
        assert encoded["magnitude"] == pytest.approx(5000.0)

    def test_value_always_keeps_the_original_units(self):
        encoded = ConfigurationToJSON(
            {"v": 5 * ureg.km}, splitUnits=True, standardize=True
        )["v"]
        assert encoded["value"] == "5 kilometer"

    def test_keep_original_units_false_drops_the_value_field(self):
        """Used when building query dicts, where the units are not needed."""
        encoded = ConfigurationToJSON(
            {"v": 5 * ureg.km}, splitUnits=True, keepOriginalUnits=False
        )["v"]
        assert set(encoded) == {"magnitude", "units"}


@pytest.mark.unit
class TestRoundTrip:
    def test_quantity_survives_a_round_trip(self):
        original = {"v": 5 * ureg.km}
        restored = JSONToConfiguration(ConfigurationToJSON(original))
        assert restored["v"].to(ureg.m).magnitude == pytest.approx(5000.0)

    def test_plain_values_survive_a_round_trip(self):
        original = {"n": 3, "s": "hello", "b": True, "f": 1.5}
        assert JSONToConfiguration(ConfigurationToJSON(original)) == original

    def test_nested_structures_survive(self):
        original = {"outer": {"inner": 2 * ureg.s, "count": 7}}
        restored = JSONToConfiguration(ConfigurationToJSON(original))
        assert restored["outer"]["count"] == 7
        assert restored["outer"]["inner"].to(ureg.s).magnitude == pytest.approx(2.0)

    def test_split_units_round_trip_recovers_original_units(self):
        encoded = ConfigurationToJSON({"v": 5 * ureg.km}, splitUnits=True)
        restored = JSONToConfiguration(encoded)
        assert restored["v"].to(ureg.km).magnitude == pytest.approx(5.0)

    def test_standardized_split_units_round_trip(self):
        encoded = ConfigurationToJSON(
            {"v": 5 * ureg.km}, splitUnits=True, standardize=True
        )
        restored = JSONToConfiguration(encoded, returnStandardize=True)
        assert str(restored["v"].units) == "meter"
        assert restored["v"].magnitude == pytest.approx(5000.0)

    def test_a_bare_dict_is_not_mistaken_for_a_quantity(self):
        """Only a 3-key magnitude/units/value dict is a quantity."""
        original = {"magnitude": 1, "units": "m"}
        restored = JSONToConfiguration(original)
        assert isinstance(restored, dict)
        assert "magnitude" in restored

    @pytest.mark.xfail(
        strict=True,
        reason="B13: returnStandardize is documented as 'return the units in MKS', "
               "but the decoder reads the 'units' field, which only holds base "
               "units when the ENCODER was given standardize=True. Asking the "
               "decoder for MKS on a non-standardized document silently returns "
               "the original units. See the consolidated findings issue.",
    )
    def test_decoder_can_standardize_on_its_own(self):
        """returnStandardize=True must mean MKS regardless of how it was encoded."""
        encoded = ConfigurationToJSON({"v": 5 * ureg.km}, splitUnits=True)
        restored = JSONToConfiguration(encoded, returnStandardize=True)
        assert str(restored["v"].units) == "meter"
        assert restored["v"].magnitude == pytest.approx(5000.0)


@pytest.mark.unit
class TestStripConfigurationUnits:
    def test_quantity_becomes_its_magnitude(self):
        assert stripConfigurationUnits({"v": 5 * ureg.km}) == {"v": 5}

    def test_standardize_strips_the_base_unit_magnitude(self):
        assert stripConfigurationUnits({"v": 5 * ureg.km}, returnStandardize=True) == {
            "v": 5000.0
        }

    def test_ignore_standardization_exempts_a_key(self):
        stripped = stripConfigurationUnits(
            {"v": 5 * ureg.km}, returnStandardize=True, ignoreStandardization=["v"]
        )
        assert stripped == {"v": 5}

    def test_unitless_values_are_untouched(self):
        assert stripConfigurationUnits({"n": 3, "s": "x"}) == {"n": 3, "s": "x"}

    def test_nested_and_list_values_are_handled(self):
        stripped = stripConfigurationUnits(
            {"a": {"b": 2 * ureg.m}, "l": [3 * ureg.s, 4]}
        )
        assert stripped == {"a": {"b": 2}, "l": [3, 4]}

    def test_output_has_no_quantity_left(self):
        stripped = stripConfigurationUnits({"a": {"b": 2 * ureg.m}})
        assert not hasattr(stripped["a"]["b"], "magnitude")


@pytest.mark.unit
class TestLoadJSON:
    """Accepts a dict, a list, a JSON string, a path, or a file object."""

    def test_dict_passes_through_unchanged(self):
        payload = {"a": 1}
        assert loadJSON(payload) is payload

    def test_list_passes_through(self):
        assert loadJSON([1, 2]) == [1, 2]

    def test_json_string(self):
        assert loadJSON('{"a": 1}') == {"a": 1}

    def test_python_repr_string_is_accepted(self):
        """Single quotes, True/False/None are translated to JSON spelling."""
        assert loadJSON("{'a': True, 'b': None, 'c': False}") == {
            "a": True,
            "b": None,
            "c": False,
        }

    def test_path_on_disk(self, tmp_path):
        path = tmp_path / "conf.json"
        path.write_text('{"a": 1}', encoding="utf-8")
        assert loadJSON(str(path)) == {"a": 1}

    def test_file_object(self, tmp_path):
        path = tmp_path / "conf.json"
        path.write_text('{"a": 1}', encoding="utf-8")
        with open(path, encoding="utf-8") as handle:
            assert loadJSON(handle) == {"a": 1}

    def test_unparseable_string_names_the_input(self):
        with pytest.raises(ValueError, match="either bad format of file does not exist"):
            loadJSON("this is not json")

    def test_unsupported_type_raises(self):
        with pytest.raises(ValueError, match="is unknown"):
            loadJSON(42)


@pytest.mark.unit
class TestSetJSONPath:
    BASE = {"a": {"b": 1}, "c": {"d": 2}}

    def test_sets_a_nested_value_by_dotted_path(self):
        assert setJSONPath(self.BASE, {"a.b": 99}) == {"a": {"b": 99}, "c": {"d": 2}}

    def test_leaves_the_original_alone_by_default(self):
        base = {"a": {"b": 1}}
        setJSONPath(base, {"a.b": 99})
        assert base == {"a": {"b": 1}}, "the input was mutated without inPlace"

    def test_inplace_mutates_and_returns_the_same_object(self):
        base = {"a": {"b": 1}}
        returned = setJSONPath(base, {"a.b": 99}, inPlace=True)
        assert returned is base
        assert base["a"]["b"] == 99

    def test_explicit_jsonpath_prefix_is_accepted(self):
        assert setJSONPath(self.BASE, {"$.c.d": 7})["c"]["d"] == 7

    def test_bare_and_prefixed_paths_agree(self):
        assert setJSONPath(self.BASE, {"c.d": 7}) == setJSONPath(self.BASE, {"$.c.d": 7})

    def test_several_paths_at_once(self):
        updated = setJSONPath(self.BASE, {"a.b": 10, "c.d": 20})
        assert updated == {"a": {"b": 10}, "c": {"d": 20}}

    def test_untouched_keys_are_preserved(self):
        updated = setJSONPath(self.BASE, {"a.b": 99})
        assert updated["c"] == {"d": 2}

    def test_a_path_that_matches_nothing_raises(self):
        """Silently ignoring a typo'd path would corrupt a variation sweep."""
        with pytest.raises(IndexError):
            setJSONPath(self.BASE, {"a.nope": 1})
