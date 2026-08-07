import pytest
from pydantic import ValidationError

from api_models import (
    ExecPayload,
    ExecResponse,
    Problem,
    RunWorkflowPayload,
    RunWorkflowResponse,
)


def test_exec_response_defaults():
    resp = ExecResponse()
    assert resp.data is None
    assert resp.problem is None


def test_exec_response_with_problem():
    resp = ExecResponse(problem=Problem(error="e", traceback="tb"))
    assert resp.problem.error == "e"
    assert resp.problem.traceback == "tb"


def test_run_workflow_payload_requires_fields():
    with pytest.raises(ValidationError):
        RunWorkflowPayload(projectName="P")  # missing workflowName


def test_run_workflow_response_fields():
    resp = RunWorkflowResponse(dispatch_id="d1", output="log")
    assert resp.dispatch_id == "d1"
    assert resp.output == "log"


def test_exec_payload_requires_code():
    with pytest.raises(ValidationError):
        ExecPayload()
