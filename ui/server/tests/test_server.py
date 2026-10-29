import os
import time

import args_settings
import dependencies
from api_models import RunWorkflowPayload
from routers import cors_info, health, workflow


def test_healthz():
    assert health.healthz() == {"status": "ok"}


def test_cors_info_shape():
    info = cors_info.cors_info()
    assert "origins" in info


def test_cors_info_answers_the_origins_the_server_started_with(monkeypatch):
    settings = args_settings.parse_args_settings(["--cors", "10.0.0.1", "--port", "9000"])
    settings.origins()
    monkeypatch.setattr(args_settings, "_settings", settings)

    assert cors_info.cors_info() == {"origins": ["http://10.0.0.1:9000"]}


def test_cors_info_has_no_origins_without_cors(monkeypatch):
    settings = args_settings.parse_args_settings([])
    settings.origins()
    monkeypatch.setattr(args_settings, "_settings", settings)

    assert cors_info.cors_info() == {"origins": None}


def test_start_workflow_delegates_to_runner_start(monkeypatch):
    calls = {}

    def fake_start(project_name, doc):
        calls["args"] = (project_name, doc)
        return {"token": "t1"}

    monkeypatch.setattr(dependencies.workflow_runner, "start", fake_start)

    doc = {"desc": {"workflow": {"solver": "s"}, "workflowName": "W"}, "resource": "/tmp/W.json"}
    resp = workflow.start_workflow(RunWorkflowPayload(projectName="P", doc=doc))

    assert calls["args"] == ("P", doc)
    assert resp.token == "t1"


def test_workflow_status_delegates_to_runner_poll(monkeypatch):
    monkeypatch.setattr(
        dependencies.workflow_runner, "poll",
        lambda token: {"status": "done", "error": "", "chunks": [{"name": "__between__", "text": "log"}]},
    )

    resp = workflow.workflow_status("t1")

    assert resp.status == "done"
    assert resp.chunks[0].text == "log"


def test_run_workflow_endpoint_through_real_runner(install_fake_hera, tmp_path):
    # Exercises the actual WorkflowRunner (not mocked) against a fake hera, so the
    # route + runner + capture path all run together: start, then poll until done.
    def on_execute(workflow_name):
        os.write(1, ("workflow %s done\n" % workflow_name).encode())
        return "dispatch-xyz"

    install_fake_hera(str(tmp_path), on_execute)

    hello_doc = {"desc": {"workflow": {"solver": "s"}, "workflowName": "hello"}, "resource": str(tmp_path / "hello.json")}
    start = workflow.start_workflow(RunWorkflowPayload(projectName="P", doc=hello_doc))
    assert start.token

    deadline = time.time() + 10.0
    resp = workflow.workflow_status(start.token)
    while resp.status == "running" and time.time() < deadline:
        time.sleep(0.01)
        resp = workflow.workflow_status(start.token)

    assert resp.status == "done"
    assert "workflow hello done" in "".join(chunk.text for chunk in resp.chunks)

