"""Validation of node parameter values (node_params + the endpoint).

The real Hermes executers need hera and MongoDB, so these tests stub the class
lookup and only check the wiring: which class is asked, what gets returned, and
when nothing is reported.
"""
import pytest

import node_params
import server
from api_models import ValidateNodeParamsPayload


class _WithChecks:
    @staticmethod
    def testParamValues(params):
        if params.get("Agent") == "H2S":
            return True, ""
        return False, "Agent 'x' doesn't exists, choose one of: H2S"


@pytest.fixture
def warmed(monkeypatch):
    monkeypatch.setattr(server.warmup, "_ready", True)


def test_bad_value_returns_the_node_message(monkeypatch):
    monkeypatch.setattr(node_params, "_executer_class", lambda node_type: _WithChecks)
    result = node_params.validate_node_params("RiskAssessment.calculateThresholds", {"Agent": "x"})
    assert result["ok"] is False
    assert "doesn't exists" in result["message"]


def test_good_value_is_ok_with_no_message(monkeypatch):
    monkeypatch.setattr(node_params, "_executer_class", lambda node_type: _WithChecks)
    result = node_params.validate_node_params("RiskAssessment.calculateThresholds", {"Agent": "H2S"})
    assert result == {"ok": True, "message": ""}


def test_class_that_is_not_a_node_executer_is_never_called(monkeypatch):
    # _executer_class only returns abstractExecuter subclasses, so anything else
    # reads as "no node" rather than as something to run.
    monkeypatch.setattr(node_params, "_executer_class", lambda node_type: None)
    assert node_params.validate_node_params("general.CopyFile", {}) == {"ok": True, "message": ""}


def test_unknown_type_reports_nothing():
    assert node_params.validate_node_params("no.such.node", {"a": 1}) == {"ok": True, "message": ""}


def test_type_that_is_not_a_dotted_name_is_never_imported():
    # The type comes from the client; a path or an expression must not be imported.
    assert node_params._executer_class("../../etc/passwd") is None
    assert node_params._executer_class("os.system('x')") is None


def test_real_executer_class_is_found_by_type():
    # No hera in the test env, so importing the module may fail; the point is that
    # the type maps to a directory that exists and carries an executer.py.
    node_dir = node_params._HERMES_ROOT / "hermes" / "Resources" / "RiskAssessment" / "calculateThresholds"
    assert (node_dir / "executer.py").is_file()


def test_endpoint_passes_type_and_params_through(monkeypatch, warmed):
    seen = {}

    def fake_validate(node_type, params):
        seen["type"] = node_type
        seen["params"] = params
        return {"ok": False, "message": "nope"}

    monkeypatch.setattr(server, "validate_node_params", fake_validate)
    resp = server.node_params_validate(ValidateNodeParamsPayload(type="general.CopyFile", params={"a": 1}))
    assert seen == {"type": "general.CopyFile", "params": {"a": 1}}
    assert resp.ok is False
    assert resp.message == "nope"


def test_endpoint_reports_nothing_before_warmup(monkeypatch):
    def fail(node_type, params):
        raise AssertionError("must not validate before hera is warm")

    monkeypatch.setattr(server, "validate_node_params", fail)
    resp = server.node_params_validate(ValidateNodeParamsPayload(type="general.CopyFile", params={}))
    assert resp.ok is True
    assert resp.message == ""
