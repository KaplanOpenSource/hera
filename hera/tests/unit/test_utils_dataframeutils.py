"""compareDataframeConfigurations: which parameters differ between datasets.

The docstring specifies three accepted input shapes, two output shapes, and
that the result contains only the parameters that differ.  These tests assert
that contract.

Fixture: two configurations agreeing on x.a and disagreeing on x.b, so the
correct answer is always "x.b, and only x.b".
"""
import pandas as pd
import pytest

from hera.utils.dataframeutils import compareDataframeConfigurations


@pytest.fixture()
def first():
    return pd.DataFrame({"parameterName": ["x.a", "x.b"], "value": [1, 2]})


@pytest.fixture()
def second():
    return pd.DataFrame({"parameterName": ["x.a", "x.b"], "value": [1, 99]})


@pytest.mark.unit
class TestInputShapes:
    """Documented: list of tuples, list of dicts, or one long-format frame."""

    def test_list_of_tuples(self, first, second):
        result = compareDataframeConfigurations([("s1", first), ("s2", second)])
        assert result.index.tolist() == ["x.b"]
        assert result.loc["x.b"].to_dict() == {"s1": 2, "s2": 99}

    def test_list_of_dicts(self, first, second):
        result = compareDataframeConfigurations(
            [
                {"datasetName": "s1", "data": first},
                {"datasetName": "s2", "data": second},
            ]
        )
        assert result.index.tolist() == ["x.b"]

    def test_long_format_frame(self, first, second):
        combined = pd.concat(
            [first.assign(datasetName="s1"), second.assign(datasetName="s2")]
        )
        result = compareDataframeConfigurations(combined)
        assert result.index.tolist() == ["x.b"]

    def test_the_three_shapes_agree(self, first, second):
        from_tuples = compareDataframeConfigurations([("s1", first), ("s2", second)])
        from_dicts = compareDataframeConfigurations(
            [{"datasetName": "s1", "data": first}, {"datasetName": "s2", "data": second}]
        )
        assert from_tuples.equals(from_dicts)

    def test_custom_dataset_name_key(self, first, second):
        result = compareDataframeConfigurations(
            [{"setName": "s1", "data": first}, {"setName": "s2", "data": second}],
            datasetName="setName",
        )
        assert result.index.tolist() == ["x.b"]


@pytest.mark.unit
class TestWhatCountsAsADifference:
    def test_identical_configurations_yield_nothing(self, first):
        result = compareDataframeConfigurations([("s1", first), ("s2", first)])
        assert result.empty

    def test_a_parameter_missing_from_one_set_counts_as_differing(self, first):
        """Documented via the count < datasetCount branch: absence is a difference."""
        partial = pd.DataFrame({"parameterName": ["x.a"], "value": [1]})
        result = compareDataframeConfigurations([("s1", first), ("s2", partial)])
        assert result.index.tolist() == ["x.b"]
        assert result.columns.tolist() == ["s1"]

    def test_agreeing_parameters_are_excluded(self, first, second):
        result = compareDataframeConfigurations([("s1", first), ("s2", second)])
        assert "x.a" not in result.index.tolist()


@pytest.mark.unit
class TestOutputShapes:
    def test_wide_is_the_default(self, first, second):
        result = compareDataframeConfigurations([("s1", first), ("s2", second)])
        assert result.columns.tolist() == ["s1", "s2"]

    def test_long_format_keeps_one_row_per_dataset(self, first, second):
        result = compareDataframeConfigurations(
            [("s1", first), ("s2", second)], longFormat=True
        )
        assert len(result) == 2
        assert set(result["datasetName"]) == {"s1", "s2"}
        assert set(result["value"]) == {2, 99}

    def test_change_dot_to_underscore_rewrites_the_index(self, first, second):
        """So the result can be fed to DataFrame.query, which rejects dots."""
        result = compareDataframeConfigurations(
            [("s1", first), ("s2", second)], changeDotToUnderscore=True
        )
        assert result.index.tolist() == ["x_b"]


@pytest.mark.unit
class TestErrorHandling:
    def test_a_dict_without_the_name_key_raises_keyerror(self, first):
        with pytest.raises(KeyError, match="does not have the key"):
            compareDataframeConfigurations([{"data": first}, {"data": first}])

    def test_a_dict_without_data_raises_keyerror(self):
        with pytest.raises(KeyError, match="does not have key 'data'"):
            compareDataframeConfigurations([{"datasetName": "s1"}])

    def test_an_unsupported_list_element_raises_valueerror(self, first):
        with pytest.raises(ValueError, match="incorrect format"):
            compareDataframeConfigurations(["not a tuple or dict"])

    def test_an_unsupported_container_raises_keyerror(self):
        with pytest.raises(KeyError, match="must be a list of tuples"):
            compareDataframeConfigurations("neither a list nor a frame")

    @pytest.mark.xfail(
        strict=True,
        reason="B4: the empty-input branch calls pandas.Dataframe (lower-case f), "
               "so the documented 'returning {}' path raises AttributeError. "
               "See the consolidated findings issue.",
    )
    def test_empty_input_returns_an_empty_frame(self):
        """The code logs 'Empty input, returning {}' and then cannot deliver it."""
        result = compareDataframeConfigurations([])
        assert isinstance(result, pd.DataFrame)
        assert result.empty
