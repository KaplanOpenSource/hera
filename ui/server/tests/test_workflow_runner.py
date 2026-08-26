import os
import time

import pytest

from workflow_runner import WorkflowRunner, RunStatus


def _wait_done(runner, token, timeout=10.0):
    """Poll a token until the run leaves the running state; return the poll result."""
    deadline = time.time() + timeout
    result = runner.poll(token)
    while result["status"] == RunStatus.RUNNING and time.time() < deadline:
        time.sleep(0.01)
        result = runner.poll(token)
    return result


def test_run_returns_dispatch_id_and_captured_stdout(install_fake_hera, tmp_path):
    def on_execute(workflow_name, scheduler):
        # Real hera writes to the fds via a subprocess / os.system, so mimic that
        # (a plain print would go to pytest's replaced sys.stdout, not fd 1).
        os.write(1, ("ran %s on %s\n" % (workflow_name, scheduler)).encode())
        return "dispatch-123"

    install_fake_hera(str(tmp_path), on_execute)

    result = WorkflowRunner().run("PROJECT", "WORKFLOW")

    assert result["dispatch_id"] == "dispatch-123"
    assert "ran WORKFLOW on local" in result["output"]


def test_run_captures_stderr_too(install_fake_hera, tmp_path):
    def on_execute(workflow_name, scheduler):
        os.write(2, b"a warning on stderr\n")
        return "d"

    install_fake_hera(str(tmp_path), on_execute)

    result = WorkflowRunner().run("PROJECT", "WORKFLOW")

    assert "a warning on stderr" in result["output"]


def test_run_uses_the_workflows_toolkit_and_local_scheduler(install_fake_hera, tmp_path):
    record = install_fake_hera(str(tmp_path))

    WorkflowRunner().run("MY_PROJECT", "MY_WORKFLOW")

    assert record.get_toolkit_kwargs == {
        "toolkitName": "SIMULATIONS_WORKFLOWS",
        "projectName": "MY_PROJECT",
    }
    assert record.execute_calls == [("MY_WORKFLOW", "local")]


def test_pythonpath_contains_files_dir_during_the_run(install_fake_hera, tmp_path):
    seen = {}

    def on_execute(workflow_name, scheduler):
        seen["pythonpath"] = os.environ.get("PYTHONPATH")
        return "d"

    install_fake_hera(str(tmp_path), on_execute)

    WorkflowRunner().run("PROJECT", "WORKFLOW")

    assert seen["pythonpath"].split(os.pathsep)[0] == str(tmp_path)


def test_run_restores_pythonpath_when_absent(install_fake_hera, tmp_path, monkeypatch):
    monkeypatch.delenv("PYTHONPATH", raising=False)
    install_fake_hera(str(tmp_path))

    WorkflowRunner().run("PROJECT", "WORKFLOW")

    assert "PYTHONPATH" not in os.environ


def test_run_restores_previous_pythonpath(install_fake_hera, tmp_path, monkeypatch):
    monkeypatch.setenv("PYTHONPATH", "/original/path")
    install_fake_hera(str(tmp_path))

    WorkflowRunner().run("PROJECT", "WORKFLOW")

    assert os.environ["PYTHONPATH"] == "/original/path"


def test_run_restores_pythonpath_even_on_error(install_fake_hera, tmp_path, monkeypatch):
    monkeypatch.delenv("PYTHONPATH", raising=False)

    def on_execute(workflow_name, scheduler):
        raise RuntimeError("workflow blew up")

    install_fake_hera(str(tmp_path), on_execute)

    with pytest.raises(RuntimeError, match="workflow blew up"):
        WorkflowRunner().run("PROJECT", "WORKFLOW")

    assert "PYTHONPATH" not in os.environ


def test_start_then_poll_reports_done_with_output(install_fake_hera, tmp_path):
    def on_execute(workflow_name, scheduler):
        os.write(1, ("ran %s\n" % workflow_name).encode())
        return "dispatch-123"

    install_fake_hera(str(tmp_path), on_execute)
    runner = WorkflowRunner()

    start = runner.start("PROJECT", "WORKFLOW")
    assert start["token"]

    result = _wait_done(runner, start["token"])
    assert result["status"] == RunStatus.DONE
    assert "ran WORKFLOW" in result["output"]


def test_start_then_poll_reports_error(install_fake_hera, tmp_path):
    def on_execute(workflow_name, scheduler):
        raise RuntimeError("workflow blew up")

    install_fake_hera(str(tmp_path), on_execute)
    runner = WorkflowRunner()

    result = _wait_done(runner, runner.start("PROJECT", "WORKFLOW")["token"])
    assert result["status"] == RunStatus.ERROR
    assert "workflow blew up" in result["error"]


def test_poll_returns_partial_output_while_running(install_fake_hera, tmp_path):
    gate = tmp_path / "gate"

    def on_execute(workflow_name, scheduler):
        os.write(1, b"partial line\n")
        # Block until the test lets us finish, so it can observe running output.
        while not gate.exists():
            time.sleep(0.01)
        return "dispatch-123"

    install_fake_hera(str(tmp_path), on_execute)
    runner = WorkflowRunner()
    token = runner.start("PROJECT", "WORKFLOW")["token"]

    # Wait until the partial output shows up while the run is still going.
    deadline = time.time() + 10
    result = runner.poll(token)
    while "partial line" not in result["output"] and time.time() < deadline:
        assert result["status"] == RunStatus.RUNNING
        time.sleep(0.01)
        result = runner.poll(token)
    assert result["status"] == RunStatus.RUNNING
    assert "partial line" in result["output"]

    # Let the run finish; the final output still has the partial line.
    gate.write_text("go")
    final = _wait_done(runner, token)
    assert final["status"] == RunStatus.DONE
    assert "partial line" in final["output"]


def test_poll_unknown_token_is_not_found():
    assert WorkflowRunner().poll("nope") == {
        "status": RunStatus.NOT_FOUND, "output": "", "error": "",
    }


def test_start_reports_busy_while_a_run_is_in_progress():
    runner = WorkflowRunner()
    # Simulate a run in progress by leaving the runner in the running state.
    runner._token = "t"
    runner._status = RunStatus.RUNNING
    assert runner.start("PROJECT", "WORKFLOW") == {"status": RunStatus.BUSY}


def test_next_run_overwrites_the_finished_slot(install_fake_hera, tmp_path):
    install_fake_hera(str(tmp_path))
    runner = WorkflowRunner()

    first = runner.start("PROJECT", "WORKFLOW")["token"]
    assert _wait_done(runner, first)["status"] == RunStatus.DONE

    second = runner.start("PROJECT", "WORKFLOW")["token"]
    assert _wait_done(runner, second)["status"] == RunStatus.DONE
    assert second != first
    assert runner.poll(first)["status"] == RunStatus.NOT_FOUND
