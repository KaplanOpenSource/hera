"""
Tests for the RFC 6902 comparisons of issue #565: `diffJSONS` for two
documents, and `compareJSONS`, which tabulates the patches from one base
document to each of the others.
"""
import math

import jsonpatch
import pytest

from hera.utils import compareJSONS, diffJSONS

pytestmark = pytest.mark.unit


def test_identical_documents_produce_an_empty_patch():
    """The CLIs read the result as `if not res`, so equal documents must be falsy."""
    assert diffJSONS({"a": 1, "b": {"c": [1, 2]}}, {"a": 1, "b": {"c": [1, 2]}}) == []


def test_changed_value_is_a_replace():
    assert diffJSONS({"a": 1}, {"a": 2}) == [{"op": "replace", "path": "/a", "value": 2}]


def test_added_and_removed_keys():
    patch = diffJSONS({"keep": 1, "gone": 2}, {"keep": 1, "new": 3})
    assert {"op": "add", "path": "/new", "value": 3} in patch
    assert {"op": "remove", "path": "/gone"} in patch


def test_nested_paths_use_json_pointer():
    patch = diffJSONS({"outer": {"inner": 1}}, {"outer": {"inner": 2}})
    assert patch == [{"op": "replace", "path": "/outer/inner", "value": 2}]


def test_patch_applied_to_source_yields_target():
    """The point of a patch: it is executable, not just a description."""
    jsonpatch = pytest.importorskip("jsonpatch")
    source = {"a": 1, "b": {"c": 2}, "d": [1, 2, 3]}
    target = {"a": 9, "b": {"c": 2, "e": 5}, "d": [1, 3]}
    assert jsonpatch.apply_patch(source, diffJSONS(source, target)) == target


def test_direction_is_source_to_target():
    """Argument order matters — the patch turns the first into the second."""
    forward = diffJSONS({"a": 1}, {"a": 2})
    backward = diffJSONS({"a": 2}, {"a": 1})
    assert forward == [{"op": "replace", "path": "/a", "value": 2}]
    assert backward == [{"op": "replace", "path": "/a", "value": 1}]


# --- compareJSONS: N documents, each as its patch from the base ---------------

def _patch(table, name):
    """The rows of `name` back as a JSON Patch (dropping the NaN that pandas
    puts where an operation has no value or from)."""
    rows = table[table["datasetName"] == name].drop(columns="datasetName").to_dict("records")
    return [{k: v for k, v in row.items() if not (isinstance(v, float) and math.isnan(v))}
            for row in rows]


def test_compareJSONS_takes_the_first_document_as_the_base_by_default():
    table = compareJSONS(run1={"dt": 0.1}, run2={"dt": 0.2}, run3={"dt": 0.1, "n": 5})
    assert table.to_dict("records") == [
        {"datasetName": "run2", "op": "replace", "path": "/dt", "value": 0.2},
        {"datasetName": "run3", "op": "add", "path": "/n", "value": 5},
    ]


def test_compareJSONS_compares_against_the_chosen_base():
    table = compareJSONS(base="run2", run1={"dt": 0.1}, run2={"dt": 0.2})
    assert table.to_dict("records") == [
        {"datasetName": "run1", "op": "replace", "path": "/dt", "value": 0.1},
    ]


def test_applying_a_documents_rows_to_the_base_rebuilds_that_document():
    docs = dict(a={"x": [1, 2], "y": {"z": 1}},
                b={"x": [2], "y": {}},
                c={"x": [1, 2, 3], "w": "new", "y": {"z": 1}})
    table = compareJSONS(**docs)
    for name in ("b", "c"):
        assert jsonpatch.apply_patch(docs["a"], _patch(table, name)) == docs[name]


def test_identical_documents_give_an_empty_table_with_the_patch_columns():
    table = compareJSONS(one={"a": 1}, two={"a": 1})
    assert table.empty
    assert list(table.columns) == ["datasetName", "op", "path", "value"]


def test_an_unknown_base_is_rejected():
    with pytest.raises(ValueError, match="'nope'"):
        compareJSONS(base="nope", one={"a": 1}, two={"a": 2})


def test_nothing_to_compare_is_an_error():
    """An empty table would read as "no differences"."""
    with pytest.raises(ValueError, match="at least one"):
        compareJSONS()


def test_an_argument_that_is_not_a_json_document_is_rejected():
    """The removed longFormat/changeDotToUnderscore flags would otherwise be
    compared as if they were documents."""
    with pytest.raises(TypeError, match="longFormat"):
        compareJSONS(longFormat=True, one={"a": 1}, two={"a": 2})
