import server
from api_models import ExecPayload, RunWorkflowPayload


def test_healthz():
    assert server.healthz() == {"status": "ok"}


def test_cors_info_shape():
    info = server.cors_info()
    assert "origins" in info


def test_exec_returns_result_value():
    resp = server.exec_code(ExecPayload(code="result = 1 + 2"))
    assert resp.problem is None
    assert resp.data == 3


def test_exec_reports_errors():
    resp = server.exec_code(ExecPayload(code="raise ValueError('boom')"))
    assert resp.data is None
    assert "ValueError" in resp.problem.error
    assert "boom" in resp.problem.error
    assert "Traceback" in resp.problem.traceback


def test_run_workflow_delegates_to_runner(monkeypatch):
    calls = {}

    def fake_run(project_name, workflow_name):
        calls["args"] = (project_name, workflow_name)
        return {"dispatch_id": "d1", "output": "some log"}

    monkeypatch.setattr(server.workflow_runner, "run", fake_run)

    resp = server.run_workflow(RunWorkflowPayload(projectName="P", workflowName="W"))

    assert calls["args"] == ("P", "W")
    assert resp.dispatch_id == "d1"
    assert resp.output == "some log"
