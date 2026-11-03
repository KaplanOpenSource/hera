import json
import os
import time

import pytest
from fastapi import HTTPException

import dependencies
from routers import workflow
from workflow_runner import WorkflowRunner, RunStatus
from workflow_runner.active_workflow_run_one import ActiveWorkflowRunOne
from workflow_runner.workflow_run_history import WorkflowRunHistory


def _doc(workflow_name):
    return {"_id": {"$oid": "id-" + workflow_name}, "desc": {"workflow": {}, "workflowName": workflow_name}, "resource": ""}


def _wait_done(runner, token, timeout=10.0):
    deadline = time.time() + timeout
    result = runner.poll(token)
    while result["status"] == RunStatus.RUNNING and time.time() < deadline:
        time.sleep(0.01)
        result = runner.poll(token)
    return result


def _write_run(tmp_path, workflow_name, run_id, start_time, status="done"):
    folder = tmp_path / "workflow_runs" / workflow_name
    folder.mkdir(parents=True, exist_ok=True)
    record = {"token": "t-" + run_id, "startTime": start_time, "endTime": start_time, "status": status,
              "error": "", "chunks": [{"name": "__between__", "text": "out " + run_id}]}
    (folder / f"{run_id}.json").write_text(json.dumps(record))


def test_run_folders_are_under_the_files_directory(install_fake_hera, tmp_path):
    install_fake_hera(str(tmp_path))
    run = ActiveWorkflowRunOne("PROJECT", _doc("WORKFLOW"))
    assert run.target_files_folder() == os.path.join(str(tmp_path), "WORKFLOW_targetFiles")
    assert run.run_history_folder() == os.path.join(str(tmp_path), "workflow_runs", "WORKFLOW")


def test_a_finished_run_is_saved_before_its_status_changes(install_fake_hera, tmp_path):
    def on_execute(workflow_name):
        os.write(1, b"saved line\n")

    install_fake_hera(str(tmp_path), on_execute)
    runner = WorkflowRunner()
    token = runner.start("PROJECT", _doc("WORKFLOW"))["token"]
    assert _wait_done(runner, token)["status"] == RunStatus.DONE

    history = WorkflowRunHistory("PROJECT", "WORKFLOW")
    summaries = history.summaries()
    assert len(summaries) == 1
    assert summaries[0]["token"] == token
    assert summaries[0]["status"] == "done"
    assert summaries[0]["startTime"] and summaries[0]["endTime"]
    record = history.read(summaries[0]["id"])
    assert "saved line" in "".join(chunk["text"] for chunk in record["chunks"])


def test_a_failed_run_is_saved_with_its_error(install_fake_hera, tmp_path):
    def on_execute(workflow_name):
        raise RuntimeError("workflow blew up")

    install_fake_hera(str(tmp_path), on_execute)
    runner = WorkflowRunner()
    token = runner.start("PROJECT", _doc("WORKFLOW"))["token"]
    assert _wait_done(runner, token)["status"] == RunStatus.ERROR

    history = WorkflowRunHistory("PROJECT", "WORKFLOW")
    record = history.read(history.summaries()[0]["id"])
    assert record["status"] == "error"
    assert "workflow blew up" in record["error"]


def test_summaries_are_newest_first_and_empty_without_a_folder(install_fake_hera, tmp_path):
    install_fake_hera(str(tmp_path))
    history = WorkflowRunHistory("PROJECT", "WORKFLOW")
    assert history.summaries() == []

    _write_run(tmp_path, "WORKFLOW", "old", "2026-01-01 10:00:00")
    _write_run(tmp_path, "WORKFLOW", "new", "2026-02-01 10:00:00")
    assert [summary["id"] for summary in history.summaries()] == ["new", "old"]


def test_read_refuses_a_run_id_outside_the_folder(install_fake_hera, tmp_path):
    install_fake_hera(str(tmp_path))
    _write_run(tmp_path, "WORKFLOW", "one", "2026-01-01 10:00:00")
    history = WorkflowRunHistory("PROJECT", "WORKFLOW")
    assert history.read("../WORKFLOW/one") is None
    assert history.read("missing") is None
    assert history.read("one")["chunks"][0]["text"] == "out one"


def test_routes_list_and_read_saved_runs(install_fake_hera, tmp_path, monkeypatch):
    monkeypatch.setattr(dependencies.warmup, "_ready", True)
    install_fake_hera(str(tmp_path))
    _write_run(tmp_path, "WORKFLOW", "one", "2026-01-01 10:00:00", status="stopped")

    runs = workflow.workflow_runs("PROJECT", "WORKFLOW").runs
    assert [(run.id, run.status) for run in runs] == [("one", "stopped")]

    record = workflow.workflow_run("PROJECT", "WORKFLOW", "one")
    assert record.chunks[0].text == "out one"

    with pytest.raises(HTTPException) as missing:
        workflow.workflow_run("PROJECT", "WORKFLOW", "missing")
    assert missing.value.status_code == 404


def test_routes_wait_for_warmup(monkeypatch):
    monkeypatch.setattr(dependencies.warmup, "_ready", False)
    with pytest.raises(HTTPException) as warming:
        workflow.workflow_runs("PROJECT", "WORKFLOW")
    assert warming.value.status_code == 503
