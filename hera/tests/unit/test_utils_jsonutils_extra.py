"""utils/jsonutils.py: JSON<->pandas flattening, variation generation, and
cross-dataset comparison."""
import pytest

from hera.utils.jsonutils import (
    JSONVariations,
    compareJSONS,
    convertJSONtoPandas,
    processJSONToPandas,
)


@pytest.mark.unit
class TestConvertJSONToPandas:
    def test_nested_scalars_become_dotted_parameter_names(self):
        result = convertJSONtoPandas({"a": {"b": 1}})
        assert list(result["parameterNameFullPath"]) == ["a.b"]
        assert list(result["value"]) == [1]

    def test_a_list_is_flattened_with_index_suffixes(self):
        result = convertJSONtoPandas({"a": {"c": [1, 2, 3]}})
        assert list(result["parameterNameFullPath"]) == ["a.c_0", "a.c_1", "a.c_2"]


@pytest.mark.unit
class TestProcessJSONToPandas:
    def test_it_uses_the_requested_column_names(self):
        result = processJSONToPandas({"nodes": {"a": {"x": 1}}}, nameColumn="p", valueColumn="v")
        assert list(result.columns) == ["p", "v"]
        assert result["p"].iloc[0] == "nodes.a.x"


@pytest.mark.unit
class TestJSONVariations:
    def test_it_yields_the_cartesian_product_of_variation_groups(self):
        results = list(JSONVariations({"a": 1, "b": 2}, [{"a": [10, 20]}, {"b": [100, 200]}]))
        assert len(results) == 4
        assert {"a": 10, "b": 100} in results
        assert {"a": 20, "b": 200} in results

    def test_parameters_in_the_same_group_change_together(self):
        results = list(JSONVariations({"a": 1, "b": 2}, [{"a": [10, 20], "b": [100, 200]}]))
        assert len(results) == 2
        assert {"a": 10, "b": 100} in results
        assert {"a": 20, "b": 200} in results

    def test_unmentioned_base_keys_are_preserved(self):
        results = list(JSONVariations({"a": 1, "untouched": "x"}, [{"a": [10]}]))
        assert results[0]["untouched"] == "x"


@pytest.mark.unit
class TestCompareJSONS:
    def test_it_compares_common_fields_across_datasets(self):
        result = compareJSONS(one={"a": 1}, two={"a": 2})
        assert list(result.columns) == ["one", "two"]
        assert result.loc["a", "one"] == 1
        assert result.loc["a", "two"] == 2
