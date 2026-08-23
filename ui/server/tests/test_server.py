import os

import server
from api_models import ExecPayload, RunWorkflowPayload


def test_healthz():
    assert server.healthz() == {"status": "ok"}


def test_cors_info_shape():
    info = server.cors_info()
    assert "origins" in info


def test_exec_reports_warming_up_until_hera_is_imported():
    # No `warmed_up` fixture here: this is what the client gets while the
    # background warmup thread is still importing hera.
    resp = server.exec_code(ExecPayload(code="result = 1 + 2"))
    assert resp.data is None
    assert resp.problem.error == "WARMING_UP"


def test_exec_returns_result_value(warmed_up):
    resp = server.exec_code(ExecPayload(code="result = 1 + 2"))
    assert resp.problem is None
    assert resp.data == 3


def test_exec_without_result_is_none(warmed_up):
    resp = server.exec_code(ExecPayload(code="x = 5"))
    assert resp.problem is None
    assert resp.data is None


def test_exec_returns_jsonable_structure(warmed_up):
    resp = server.exec_code(ExecPayload(code="result = {'a': [1, 2], 'b': 'x'}"))
    assert resp.problem is None
    assert resp.data == {"a": [1, 2], "b": "x"}


def test_exec_reports_errors(warmed_up):
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


def test_run_workflow_endpoint_through_real_runner(install_fake_hera, tmp_path):
    # Exercises the actual WorkflowRunner (not mocked) against a fake hera, so the
    # route + runner + capture path all run together.
    def on_execute(workflow_name, scheduler):
        os.write(1, ("workflow %s done\n" % workflow_name).encode())
        return "dispatch-xyz"

    install_fake_hera(str(tmp_path), on_execute)

    resp = server.run_workflow(RunWorkflowPayload(projectName="P", workflowName="hello"))

    assert resp.dispatch_id == "dispatch-xyz"
    assert "workflow hello done" in resp.output


def test_truncate_for_log_keeps_short_values():
    assert server._truncate_for_log("hi") == "'hi'"


def test_truncate_for_log_cuts_long_values():
    truncated = server._truncate_for_log("x" * 1000)
    assert "truncated" in truncated
    assert len(truncated) < 1000
