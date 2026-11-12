import os
import signal
import time

import pytest

from workflow_runner import WorkflowRunner, RunStatus
from workflow_runner.active_workflow_run_list import ActiveWorkflowRunList
from workflow_runner.active_workflow_run_one import ActiveWorkflowRunOne


def _doc(workflow_name, resource=""):
    """A minimal sent workflow document: _id, desc.workflow, desc.workflowName, resource."""
    return {
        "_id": {"$oid": "id-" + workflow_name},
        "desc": {"workflow": {}, "workflowName": workflow_name},
        "resource": resource,
    }


def _run(project_name, doc):
    """Run one workflow in the foreground and return its result."""
    return ActiveWorkflowRunOne(project_name, doc).run()


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


def test_run_returns_no_dispatch_id_and_captured_stdout(install_fake_hera, tmp_path):
    def on_execute(workflow_name):
        # Real hera writes to the fds via a subprocess / os.system, so mimic that
        # (a plain print would go to pytest's replaced sys.stdout, not fd 1).
        os.write(1, ("ran %s\n" % workflow_name).encode())

    install_fake_hera(str(tmp_path), on_execute)

    result = _run("PROJECT", _doc("WORKFLOW"))

    # The in-process run uses the local scheduler with the flat layout, no dispatch id.
    assert result.dispatch_id is None
    assert "ran WORKFLOW" in _joined(result.chunks)


def test_run_captures_stderr_too(install_fake_hera, tmp_path):
    def on_execute(workflow_name):
        os.write(2, b"a warning on stderr\n")
        return "d"

    install_fake_hera(str(tmp_path), on_execute)

    result = _run("PROJECT", _doc("WORKFLOW"))

    assert "a warning on stderr" in _joined(result.chunks)


def test_run_executes_the_named_workflow_in_process(install_fake_hera, tmp_path):
    def on_execute(workflow_name):
        os.write(1, ("executing %s\n" % workflow_name).encode())
        return "d"

    install_fake_hera(str(tmp_path), on_execute)

    result = _run("MY_PROJECT", _doc("MY_WORKFLOW"))

    assert "executing MY_WORKFLOW" in _joined(result.chunks)


def test_pythonpath_contains_files_dir_during_the_run(install_fake_hera, tmp_path):
    def on_execute(workflow_name):
        first = os.environ.get("PYTHONPATH", "").split(os.pathsep)[0]
        os.write(1, ("pythonpath0=%s\n" % first).encode())
        return "d"

    install_fake_hera(str(tmp_path), on_execute)

    result = _run("PROJECT", _doc("WORKFLOW"))

    assert ("pythonpath0=%s" % tmp_path) in _joined(result.chunks)


def test_chunks_cover_the_whole_output(install_fake_hera, tmp_path):
    def on_execute(workflow_name):
        os.write(1, b"hello chunks\n")
        return "d"

    install_fake_hera(str(tmp_path), on_execute)

    result = _run("PROJECT", _doc("WORKFLOW"))

    assert "hello chunks" in _joined(result.chunks)


def test_run_raises_on_workflow_error(install_fake_hera, tmp_path):
    def on_execute(workflow_name):
        raise RuntimeError("workflow blew up")

    install_fake_hera(str(tmp_path), on_execute)

    with pytest.raises(RuntimeError, match="workflow blew up"):
        _run("PROJECT", _doc("WORKFLOW"))


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
    # Simulate a run in progress: a run in the list that is still running, never started.
    runner._runs.add(ActiveWorkflowRunOne("PROJECT", _doc("OTHER")))
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


def _wait_for_output(runner, token, text, timeout=10.0):
    """Poll a running token until its output holds text."""
    deadline = time.time() + timeout
    result = runner.poll(token)
    while text not in _joined(result["chunks"] or []) and time.time() < deadline:
        time.sleep(0.01)
        result = runner.poll(token)
    return result


def test_stop_ends_a_running_run_and_frees_the_runner(install_fake_hera, tmp_path):
    def on_execute(workflow_name):
        os.write(1, b"started\n")
        # Run until stopped.
        while True:
            time.sleep(0.01)

    install_fake_hera(str(tmp_path), on_execute)
    runner = WorkflowRunner()
    token = runner.start("PROJECT", _doc("WORKFLOW"))["token"]
    assert "started" in _joined(_wait_for_output(runner, token, "started")["chunks"])

    runner.stop(token)
    result = _wait_done(runner, token)
    assert result["status"] == RunStatus.STOPPED
    # The output up to the stop is kept.
    assert "started" in _joined(result["chunks"])

    # The runner is free again: a new run can start. Stop it too, so it does not run on.
    second = runner.start("PROJECT", _doc("OTHER"))["token"]
    assert second
    runner.stop(second)
    assert _wait_done(runner, second)["status"] == RunStatus.STOPPED


def test_stop_unknown_token_is_not_found():
    assert WorkflowRunner().stop("nope")["status"] == RunStatus.NOT_FOUND


def test_run_keeps_the_doc_id_and_pid(install_fake_hera, tmp_path):
    install_fake_hera(str(tmp_path))
    run = ActiveWorkflowRunOne("PROJECT", _doc("WORKFLOW"))
    run.run()
    assert run.doc_id == "id-WORKFLOW"
    assert run.pid is not None


def test_a_child_that_dies_without_a_result_is_an_error(install_fake_hera, tmp_path):
    def on_execute(workflow_name):
        os._exit(3)

    install_fake_hera(str(tmp_path), on_execute)
    with pytest.raises(RuntimeError, match="exit code 3"):
        _run("PROJECT", _doc("WORKFLOW"))


def test_list_finds_a_run_by_token_and_by_project_and_doc_id():
    runs = ActiveWorkflowRunList()
    run = ActiveWorkflowRunOne("PROJECT", _doc("WORKFLOW"))
    runs.add(run)
    assert runs.get_by_token(run.token) is run
    assert runs.get_by_doc_id("PROJECT", "id-WORKFLOW") is run
    assert runs.get_by_doc_id("OTHER_PROJECT", "id-WORKFLOW") is None


def test_list_replaces_the_old_run_of_the_same_doc():
    runs = ActiveWorkflowRunList()
    first = ActiveWorkflowRunOne("PROJECT", _doc("WORKFLOW"))
    second = ActiveWorkflowRunOne("PROJECT", _doc("WORKFLOW"))
    runs.add(first)
    runs.add(second)
    assert runs.get_by_token(first.token) is None
    assert runs.get_by_doc_id("PROJECT", "id-WORKFLOW") is second


def _pid_alive(pid):
    """True while a process with this pid exists (a zombie counts as gone)."""
    try:
        with open(f"/proc/{pid}/stat") as handle:
            return handle.read().split(") ")[1][0] != "Z"
    except FileNotFoundError:
        return False


def test_stop_ends_a_slow_os_command_task(install_fake_hera, tmp_path, request):
    pid_file = tmp_path / "sleep.pid"
    # Like uvicorn: the server catches SIGTERM, and a forked child copies that handler.
    old_handler = signal.signal(signal.SIGTERM, lambda signum, frame: None)
    request.addfinalizer(lambda: signal.signal(signal.SIGTERM, old_handler))

    def on_execute(workflow_name):
        # Like the Slow Task node: RunOsCommand runs each command with os.system.
        os.system("echo 'starting slow task'")
        os.system(f"sh -c 'echo $$ > {pid_file}; exec sleep 30'")
        os.system("echo 'slow task done'")

    install_fake_hera(str(tmp_path), on_execute)
    runner = WorkflowRunner()
    token = runner.start("PROJECT", _doc("WORKFLOW"))["token"]

    _wait_for_output(runner, token, "starting slow task")
    deadline = time.time() + 10
    while not pid_file.exists() and time.time() < deadline:
        time.sleep(0.01)
    for _ in range(3):
        assert runner.poll(token)["status"] == RunStatus.RUNNING
        time.sleep(0.1)

    run = runner._runs.get_by_token(token)
    child_pid = run.pid
    sleep_pid = int(pid_file.read_text())
    assert _pid_alive(child_pid)
    assert _pid_alive(sleep_pid)

    runner.stop(token)
    result = _wait_done(runner, token)

    assert result["status"] == RunStatus.STOPPED
    assert "slow task done" not in _joined(result["chunks"])
    assert not _pid_alive(child_pid)
    assert not _pid_alive(sleep_pid)
