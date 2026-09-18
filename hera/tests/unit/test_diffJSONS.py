"""
Tests for `hera.utils.diffJSONS` — the RFC 6902 two-document comparison added
for issue #565.

`compareJSONS` keeps its N-document table; this covers only the pairwise case
and the "no differences" signal the CLIs branch on.
"""
import pytest

from hera.utils import diffJSONS


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
