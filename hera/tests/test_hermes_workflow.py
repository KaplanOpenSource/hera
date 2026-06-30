"""
Tests for running hermes workflows through hera's hermesWorkflowToolkit.

Covers:
  - Building and executing a minimal single-node hermes workflow with the
    local Luigi scheduler, verifying it completes and writes its targets.

These tests are skipped when hermes or luigi are not importable so the suite
remains green in environments without the full hermes stack.
"""

import os
import subprocess
import sys
import tempfile

import pytest


# ---------------------------------------------------------------------------
# Skip entire module if hermes / luigi unavailable
# ---------------------------------------------------------------------------

hermes = pytest.importorskip("hermes")
pytest.importorskip("luigi")


# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------

_HERMES_ROOT = os.path.dirname(hermes.__file__) + "/.."

_SIMPLE_WORKFLOW_JSON = {
    "workflow": {
        "nodes": {
            "RunPythonCode": {
                "Execution": {
                    "input_parameters": {
                        "ModulePath": "hello_mod",
                        "ClassName": "Hello",
                        "MethodName": "run",
                        "Parameters": {},
                    }
                },
                "type": "general.RunPythonCode",
            }
        }
    }
}

_HELLO_MODULE = (
    "class Hello:\n"
    "    def run(self):\n"
    "        print('hello from hermes')\n"
)


def _run_workflow(workdir: str, workflow_json: dict, module_name: str = "Workflow1"):
    """Build and execute a hermes workflow in *workdir* using the local scheduler.

    Returns (returncode, stdout, stderr).
    """
    from hermes import workflow

    wf = workflow(workflow_json, workdir, Resource_path=workdir)
    workflow_script = os.path.join(workdir, f"{module_name}.py")
    with open(workflow_script, "w") as fp:
        fp.write(wf.build(buildername=workflow.BUILDER_LUIGI))

    cmd = (
        f"{sys.executable} -m luigi"
        f" --module {module_name} finalnode_xx_0 --local-scheduler"
    )
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [_HERMES_ROOT, workdir, env.get("PYTHONPATH", "")]
    )
    result = subprocess.run(
        cmd, shell=True, cwd=workdir, env=env,
        capture_output=True, text=True, timeout=120,
    )
    return result.returncode, result.stdout, result.stderr


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestSimpleHermesWorkflow:
    """End-to-end tests running a hermes workflow with the local scheduler."""

    def test_single_node_workflow_completes(self, tmp_path):
        """A single RunPythonCode node should run to completion."""
        workdir = str(tmp_path)
        with open(os.path.join(workdir, "hello_mod.py"), "w") as fp:
            fp.write(_HELLO_MODULE)

        rc, stdout, stderr = _run_workflow(workdir, _SIMPLE_WORKFLOW_JSON)
        assert rc == 0, f"Workflow exited with code {rc}.\nstderr:\n{stderr}"

    def test_luigi_reports_success(self, tmp_path):
        """Luigi must report no failed tasks in its execution summary."""
        workdir = str(tmp_path)
        with open(os.path.join(workdir, "hello_mod.py"), "w") as fp:
            fp.write(_HELLO_MODULE)

        _, _, stderr = _run_workflow(workdir, _SIMPLE_WORKFLOW_JSON)
        assert "this progress looks :)" in stderr.lower(), (
            "Luigi did not report successful completion.\nstderr:\n" + stderr
        )

    def test_target_files_are_created(self, tmp_path):
        """Luigi target JSON files must be written for each node."""
        workdir = str(tmp_path)
        with open(os.path.join(workdir, "hello_mod.py"), "w") as fp:
            fp.write(_HELLO_MODULE)

        _run_workflow(workdir, _SIMPLE_WORKFLOW_JSON)

        target_dir = os.path.join(workdir, "Workflow1_targetFiles")
        assert os.path.isdir(target_dir), f"Target directory not found: {target_dir}"
        assert os.path.isfile(os.path.join(target_dir, "RunPythonCode_0.json"))
        assert os.path.isfile(os.path.join(target_dir, "finalnode_xx_0.json"))

    def test_workflow_is_idempotent(self, tmp_path):
        """Running the same workflow twice should succeed both times (Luigi re-uses targets)."""
        workdir = str(tmp_path)
        with open(os.path.join(workdir, "hello_mod.py"), "w") as fp:
            fp.write(_HELLO_MODULE)

        rc1, _, stderr1 = _run_workflow(workdir, _SIMPLE_WORKFLOW_JSON)
        rc2, _, stderr2 = _run_workflow(workdir, _SIMPLE_WORKFLOW_JSON)
        assert rc1 == 0, f"First run failed.\nstderr:\n{stderr1}"
        assert rc2 == 0, f"Second run failed.\nstderr:\n{stderr2}"
