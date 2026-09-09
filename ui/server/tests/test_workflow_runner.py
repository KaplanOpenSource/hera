import os
import time

import pytest

from workflow_runner import WorkflowRunner, RunStatus


def _doc(workflow_name, resource=""):
    """A minimal sent workflow document: desc.workflow, desc.workflowName, resource."""
    return {"desc": {"workflow": {}, "workflowName": workflow_name}, "resource": resource}


def _joined(chunks):
    """Join per-task chunks into the flat text the run produced."""
    return "".join(chunk["text"] for chunk in chunks)


def _wait_done(runner, token, timeout=10.0):
    """Poll a token until the run leaves the running state; return the poll result."""
    deadline = time.time() + timeout
    result = runner.poll(token)
    while result["status"] == RunStatus.RUNNING and time.time() < deadline:
        time.sleep(0.01)
        result = runner.poll(token)
    return result


def test_run_returns_dispatch_id_and_captured_stdout(install_fake_hera, tmp_path):
    def on_execute(workflow_name):
        # Real hera writes to the fds via a subprocess / os.system, so mimic that
        # (a plain print would go to pytest's replaced sys.stdout, not fd 1).
        os.write(1, ("ran %s\n" % workflow_name).encode())
        return "dispatch-123"

    install_fake_hera(str(tmp_path), on_execute)

    result = WorkflowRunner().run("PROJECT", _doc("WORKFLOW"))

    assert result.dispatch_id == "dispatch-123"
    assert "ran WORKFLOW" in _joined(result.chunks)


def test_run_captures_stderr_too(install_fake_hera, tmp_path):
    def on_execute(workflow_name):
        os.write(2, b"a warning on stderr\n")
        return "d"

    install_fake_hera(str(tmp_path), on_execute)

    result = WorkflowRunner().run("PROJECT", _doc("WORKFLOW"))

    assert "a warning on stderr" in _joined(result.chunks)


def test_run_executes_the_named_workflow_in_process(install_fake_hera, tmp_path):
    def on_execute(workflow_name):
        os.write(1, ("executing %s\n" % workflow_name).encode())
        return "d"

    install_fake_hera(str(tmp_path), on_execute)

    result = WorkflowRunner().run("MY_PROJECT", _doc("MY_WORKFLOW"))

    assert "executing MY_WORKFLOW" in _joined(result.chunks)


def test_pythonpath_contains_files_dir_during_the_run(install_fake_hera, tmp_path):
    def on_execute(workflow_name):
        first = os.environ.get("PYTHONPATH", "").split(os.pathsep)[0]
        os.write(1, ("pythonpath0=%s\n" % first).encode())
        return "d"

    install_fake_hera(str(tmp_path), on_execute)

    result = WorkflowRunner().run("PROJECT", _doc("WORKFLOW"))

    assert ("pythonpath0=%s" % tmp_path) in _joined(result.chunks)


def test_chunks_cover_the_whole_output(install_fake_hera, tmp_path):
    def on_execute(workflow_name):
        os.write(1, b"hello chunks\n")
        return "d"

    install_fake_hera(str(tmp_path), on_execute)

    result = WorkflowRunner().run("PROJECT", _doc("WORKFLOW"))

    assert "hello chunks" in _joined(result.chunks)


def test_run_raises_on_workflow_error(install_fake_hera, tmp_path):
    def on_execute(workflow_name):
        raise RuntimeError("workflow blew up")

    install_fake_hera(str(tmp_path), on_execute)

    with pytest.raises(RuntimeError, match="workflow blew up"):
        WorkflowRunner().run("PROJECT", _doc("WORKFLOW"))


def test_start_then_poll_reports_done_with_output(install_fake_hera, tmp_path):
    def on_execute(workflow_name):
        os.write(1, ("ran %s\n" % workflow_name).encode())
        return "dispatch-123"

    install_fake_hera(str(tmp_path), on_execute)
    runner = WorkflowRunner()

    start = runner.start("PROJECT", _doc("WORKFLOW"))
    assert start["token"]

    result = _wait_done(runner, start["token"])
    assert result["status"] == RunStatus.DONE
    assert "ran WORKFLOW" in _joined(result["chunks"])


def test_start_then_poll_reports_error(install_fake_hera, tmp_path):
    def on_execute(workflow_name):
        raise RuntimeError("workflow blew up")

    install_fake_hera(str(tmp_path), on_execute)
    runner = WorkflowRunner()

    result = _wait_done(runner, runner.start("PROJECT", _doc("WORKFLOW"))["token"])
    assert result["status"] == RunStatus.ERROR
    assert "workflow blew up" in result["error"]


def test_poll_returns_partial_output_while_running(install_fake_hera, tmp_path):
    gate = tmp_path / "gate"

    def on_execute(workflow_name):
        os.write(1, b"partial line\n")
        # Block until the test lets us finish, so it can observe running output.
        while not gate.exists():
            time.sleep(0.01)
        return "dispatch-123"

    install_fake_hera(str(tmp_path), on_execute)
    runner = WorkflowRunner()
    token = runner.start("PROJECT", _doc("WORKFLOW"))["token"]

    # Wait until the partial output shows up while the run is still going.
    deadline = time.time() + 10
    result = runner.poll(token)
    while "partial line" not in _joined(result["chunks"] or []) and time.time() < deadline:
        assert result["status"] == RunStatus.RUNNING
        time.sleep(0.01)
        result = runner.poll(token)
    assert result["status"] == RunStatus.RUNNING
    assert "partial line" in _joined(result["chunks"])

    # Let the run finish; the final output still has the partial line.
    gate.write_text("go")
    final = _wait_done(runner, token)
    assert final["status"] == RunStatus.DONE
    assert "partial line" in _joined(final["chunks"])


def test_poll_unknown_token_is_not_found():
    assert WorkflowRunner().poll("nope") == {
        "status": RunStatus.NOT_FOUND, "error": "", "chunks": None,
    }


def test_start_reports_busy_while_a_run_is_in_progress():
    runner = WorkflowRunner()
    # Simulate a run in progress by leaving the runner in the running state.
    runner._token = "t"
    runner._status = RunStatus.RUNNING
    assert runner.start("PROJECT", _doc("WORKFLOW")) == {"status": RunStatus.BUSY}


def test_next_run_overwrites_the_finished_slot(install_fake_hera, tmp_path):
    install_fake_hera(str(tmp_path))
    runner = WorkflowRunner()

    first = runner.start("PROJECT", _doc("WORKFLOW"))["token"]
    assert _wait_done(runner, first)["status"] == RunStatus.DONE

    second = runner.start("PROJECT", _doc("WORKFLOW"))["token"]
    assert _wait_done(runner, second)["status"] == RunStatus.DONE
    assert second != first
    assert runner.poll(first)["status"] == RunStatus.NOT_FOUND
