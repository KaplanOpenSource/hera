import os

import pytest

from workflow_runner import WorkflowRunner


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
