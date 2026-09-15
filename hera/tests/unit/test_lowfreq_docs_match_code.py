"""
Regression guard for issue #843 — the MeteoLowFreq documentation named methods
that do not exist (e.g. `plotSeasonalHourly`) and used wrong argument names.

Every `lf.<accessor>.<method>(...)` call in the doc's Python examples is checked
against the real classes, arguments included. Nothing is instantiated, so no
MongoDB is needed.
"""
import ast
import inspect
import pathlib
import re

import pytest

DOC = pathlib.Path(__file__).resolve().parents[3] / "docs/toolkits/measurements/meteorology/lowfreq.md"

# The object `lf.<accessor>` resolves to in the examples.
ACCESSORS = {
    "analysis": "hera.measurements.meteorology.lowfreqdata.analysis:analysis",
    "presentation.dailyPlots": "hera.measurements.meteorology.lowfreqdata.presentationLayer:DailyPlots",
    "presentation.seasonalPlots": "hera.measurements.meteorology.lowfreqdata.presentationLayer:SeasonalPlots",
}

CALL = re.compile(r"\blf\.((?:\w+\.)*\w+)\.(\w+)\(")


def _resolve(path):
    import importlib

    module, _, name = ACCESSORS[path].partition(":")
    return getattr(importlib.import_module(module), name)


def _doc_calls():
    text = DOC.read_text(encoding="utf-8")
    code = "\n".join(re.findall(r"```python\n(.*?)```", text, re.DOTALL))
    return sorted({(m.group(1), m.group(2)) for m in CALL.finditer(code)})


def test_doc_exists_and_has_examples():
    assert DOC.is_file(), f"{DOC} is missing"
    assert _doc_calls(), "no `lf.<accessor>.<method>(...)` examples found in the doc"


@pytest.mark.parametrize("accessor,method", _doc_calls(), ids=lambda v: str(v))
def test_documented_method_exists(accessor, method):
    assert accessor in ACCESSORS, (
        f"doc uses lf.{accessor}.{method}(), which this test does not know how to "
        f"resolve — add it to ACCESSORS"
    )
    owner = _resolve(accessor)
    assert hasattr(owner, method), f"doc calls lf.{accessor}.{method}() but {owner.__name__} has no such method"


def test_documented_keyword_arguments_exist():
    """Catches the `field=` / `density=` style errors the issue reported."""
    text = DOC.read_text(encoding="utf-8")
    code = "\n".join(re.findall(r"```python\n(.*?)```", text, re.DOTALL))
    bad = []
    for node in ast.walk(ast.parse(code)):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
            continue
        target = ast.unparse(node.func)
        if not target.startswith("lf."):
            continue
        accessor, _, method = target[len("lf."):].rpartition(".")
        if accessor not in ACCESSORS:
            continue
        owner = _resolve(accessor)
        if not hasattr(owner, method):
            continue  # reported by test_documented_method_exists
        accepted = set(inspect.signature(getattr(owner, method)).parameters)
        for keyword in node.keywords:
            if keyword.arg is not None and keyword.arg not in accepted:
                bad.append(f"lf.{accessor}.{method}(): unknown argument '{keyword.arg}'")
    assert not bad, "\n".join(bad)
