import os
import time

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


def test_exec_without_result_is_none():
    resp = server.exec_code(ExecPayload(code="x = 5"))
    assert resp.problem is None
    assert resp.data is None


def test_exec_returns_jsonable_structure():
    resp = server.exec_code(ExecPayload(code="result = {'a': [1, 2], 'b': 'x'}"))
    assert resp.problem is None
    assert resp.data == {"a": [1, 2], "b": "x"}


def test_exec_reports_errors():
    resp = server.exec_code(ExecPayload(code="raise ValueError('boom')"))
    assert resp.data is None
    assert "ValueError" in resp.problem.error
    assert "boom" in resp.problem.error
    assert "Traceback" in resp.problem.traceback


def test_start_workflow_delegates_to_runner_start(monkeypatch):
    calls = {}

    def fake_start(project_name, workflow_name):
        calls["args"] = (project_name, workflow_name)
        return {"token": "t1"}

    monkeypatch.setattr(server.workflow_runner, "start", fake_start)

    resp = server.start_workflow(RunWorkflowPayload(projectName="P", workflowName="W"))

    assert calls["args"] == ("P", "W")
    assert resp.token == "t1"


def test_workflow_status_delegates_to_runner_poll(monkeypatch):
    monkeypatch.setattr(
        server.workflow_runner, "poll",
        lambda token: {"status": "done", "output": "log", "error": ""},
    )

    resp = server.workflow_status("t1")

    assert resp.status == "done"
    assert resp.output == "log"


def test_run_workflow_endpoint_through_real_runner(install_fake_hera, tmp_path):
    # Exercises the actual WorkflowRunner (not mocked) against a fake hera, so the
    # route + runner + capture path all run together: start, then poll until done.
    def on_execute(workflow_name, scheduler):
        os.write(1, ("workflow %s done\n" % workflow_name).encode())
        return "dispatch-xyz"

    install_fake_hera(str(tmp_path), on_execute)

    start = server.start_workflow(RunWorkflowPayload(projectName="P", workflowName="hello"))
    assert start.token

    deadline = time.time() + 10.0
    resp = server.workflow_status(start.token)
    while resp.status == "running" and time.time() < deadline:
        time.sleep(0.01)
        resp = server.workflow_status(start.token)

    assert resp.status == "done"
    assert "workflow hello done" in resp.output


def test_truncate_for_log_keeps_short_values():
    assert server._truncate_for_log("hi") == "'hi'"


def test_truncate_for_log_cuts_long_values():
    truncated = server._truncate_for_log("x" * 1000)
    assert "truncated" in truncated
    assert len(truncated) < 1000
