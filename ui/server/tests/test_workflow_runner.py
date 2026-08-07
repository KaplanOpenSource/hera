import os
import sys
import types

import pytest


def _install_fake_hera(monkeypatch, files_directory, on_execute):
    """Inject a fake `hera` module so WorkflowRunner.run doesn't need real hera."""
    toolkit = types.SimpleNamespace(
        FilesDirectory=files_directory,
        executeWorkflowFromDB=on_execute,
    )
    toolkit_home = types.SimpleNamespace(
        SIMULATIONS_WORKFLOWS="SIMULATIONS_WORKFLOWS",
        getToolkit=lambda toolkitName, projectName: toolkit,
    )
    fake_hera = types.ModuleType("hera")
    fake_hera.toolkitHome = toolkit_home
    monkeypatch.setitem(sys.modules, "hera", fake_hera)


def test_run_returns_dispatch_id_and_captured_output(monkeypatch, tmp_path):
    def on_execute(workflow_name, scheduler):
        # Real hera writes to the fds via a subprocess / os.system, so mimic that
        # (a plain print would go to pytest's replaced sys.stdout, not fd 1).
        os.write(1, ("ran %s on %s\n" % (workflow_name, scheduler)).encode())
        return "dispatch-123"

    _install_fake_hera(monkeypatch, str(tmp_path), on_execute)
    from workflow_runner import WorkflowRunner

    result = WorkflowRunner().run("PROJECT", "WORKFLOW")

    assert result["dispatch_id"] == "dispatch-123"
    assert "ran WORKFLOW on local" in result["output"]


def test_run_restores_pythonpath_when_absent(monkeypatch, tmp_path):
    monkeypatch.delenv("PYTHONPATH", raising=False)
    _install_fake_hera(monkeypatch, str(tmp_path), lambda name, scheduler: "d")
    from workflow_runner import WorkflowRunner

    WorkflowRunner().run("PROJECT", "WORKFLOW")

    assert "PYTHONPATH" not in os.environ


def test_run_restores_previous_pythonpath(monkeypatch, tmp_path):
    monkeypatch.setenv("PYTHONPATH", "/original/path")
    _install_fake_hera(monkeypatch, str(tmp_path), lambda name, scheduler: "d")
    from workflow_runner import WorkflowRunner

    WorkflowRunner().run("PROJECT", "WORKFLOW")

    assert os.environ["PYTHONPATH"] == "/original/path"


def test_run_restores_pythonpath_even_on_error(monkeypatch, tmp_path):
    monkeypatch.delenv("PYTHONPATH", raising=False)

    def on_execute(name, scheduler):
        raise RuntimeError("workflow blew up")

    _install_fake_hera(monkeypatch, str(tmp_path), on_execute)
    from workflow_runner import WorkflowRunner

    with pytest.raises(RuntimeError, match="workflow blew up"):
        WorkflowRunner().run("PROJECT", "WORKFLOW")

    assert "PYTHONPATH" not in os.environ
