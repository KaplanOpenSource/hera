"""Query builders: pandas query strings and mongoengine-style key flattening.

Assertions come from the docstrings in hera/utils/query.py, which specify
both transforms precisely.
"""
import pytest

from hera.utils.query import andClause, dictToMongoQuery


@pytest.mark.unit
class TestAndClause:
    """Builds a pandas .query() string, one clause per keyword, joined by 'and'."""

    def test_string_value_is_quoted(self):
        assert andClause(name="A") == "name == 'A'"

    def test_numeric_value_is_bare(self):
        assert andClause(n=5) == "n == 5"

    def test_list_value_becomes_membership(self):
        assert andClause(n=[1, 2]) == "n in [1, 2]"

    def test_dict_value_supplies_its_own_operator(self):
        assert andClause(n={"operator": ">=", "value": 3}) == "n >= 3"

    def test_clauses_are_joined_with_and(self):
        assert andClause(a="x", b=2) == "a == 'x' and b == 2"

    def test_excluded_fields_are_dropped(self):
        assert andClause(excludeFields=["a"], a="x", b=2) == "b == 2"

    def test_no_keywords_yields_an_empty_string(self):
        """An empty clause, not 'True' or a stray 'and'."""
        assert andClause() == ""

    def test_excluding_every_field_yields_an_empty_string(self):
        assert andClause(excludeFields=["a", "b"], a="x", b=2) == ""

    def test_output_is_a_usable_pandas_query(self):
        """The point of the function: the string must actually run."""
        import pandas as pd

        frame = pd.DataFrame({"a": ["x", "y", "x"], "b": [2, 2, 9]})
        selected = frame.query(andClause(a="x", b=2))
        assert len(selected) == 1
        assert selected.index.tolist() == [0]


@pytest.mark.unit
class TestDictToMongoQuery:
    """Flattens a nested dict into mongoengine's double-underscore keys."""

    def test_flat_dict_is_unchanged(self):
        assert dictToMongoQuery({"a": 1}) == {"a": 1}

    def test_nested_dict_is_joined_with_double_underscore(self):
        assert dictToMongoQuery({"config": {"model": "LSM"}}) == {"config__model": "LSM"}

    def test_prefix_is_prepended(self):
        assert dictToMongoQuery({"a": 1}, prefix="desc") == {"desc__a": 1}

    def test_empty_dict_yields_empty_dict(self):
        assert dictToMongoQuery({}) == {}

    def test_deep_nesting(self):
        assert dictToMongoQuery({"a": {"b": {"c": 7}}}) == {"a__b__c": 7}

    def test_prefix_exclude_key_is_skipped(self):
        """{'desc': {...}} must not become 'desc__desc__...'."""
        assert dictToMongoQuery({"desc": {"f": 1}}) == {"f": 1}

    def test_prefix_exclude_is_configurable(self):
        assert dictToMongoQuery({"wrap": {"f": 1}}, prefixExclude="wrap") == {"f": 1}

    def test_list_values_are_indexed(self):
        """Documented behaviour: each element gets its own positional key.

        This is what lets a stored version tuple be queried element-wise.
        """
        assert dictToMongoQuery({"version": [0, 0, 1]}) == {
            "version__0": 0,
            "version__1": 0,
            "version__2": 1,
        }

    def test_list_of_dicts_indexes_then_recurses(self):
        assert dictToMongoQuery({"a": [{"b": 1}, {"b": 2}]}) == {
            "a__0__b": 1,
            "a__1__b": 2,
        }

    def test_trailing_type_key_is_escaped(self):
        """'...__type' collides with mongoengine's type discrimination.

        An extra '__' is appended to escape it.
        """
        assert dictToMongoQuery({"a": {"type": "X"}}) == {"a__type__": "X"}
        assert dictToMongoQuery({"type": "X"}, prefix="desc") == {"desc__type__": "X"}

    def test_bare_top_level_type_is_not_escaped(self):
        """'type' alone does not end with '__type', so no escape applies.

        Asymmetric with the prefixed case above; pinned so the asymmetry is a
        decision rather than an accident.
        """
        assert dictToMongoQuery({"type": "X"}) == {"type": "X"}

    @pytest.mark.xfail(
        strict=True,
        reason="B3: list flattening is applied to mongoengine operator suffixes too, "
               "so station__in=['A','C'] becomes desc__station__in__0/__1 and the "
               "query silently matches nothing. See the consolidated findings issue.",
    )
    def test_list_valued_operator_suffix_survives(self):
        """A list-valued query operator must stay an operator.

        'station__in': ['A', 'C'] has to reach mongo as a single __in clause.
        Indexing the list turns it into a lookup for desc.station.in.0, which
        matches no document and raises nothing -- the worst combination.
        """
        built = dictToMongoQuery({"station__in": ["A", "C"]}, prefix="desc")
        assert built == {"desc__station__in": ["A", "C"]}
