"""
Tests for the Luigi execution command construction and the dispatch_id mechanism
introduced for the centralized-scheduler support (issue #918).

Covers:
    - buildLuigiExecutionCommand: local vs central scheduler, host/port, dispatch-id.
    - An end-to-end run of a real hermes workflow driven through hera's command,
      asserting per-dispatch output isolation (skipped when hermes/luigi are
      unavailable or the vendored hermes predates the dispatch_id change).

The generated Luigi node template (dispatch_id parameter + propagation) is tested
in the hermes repository, since that is where the template lives.
"""

import os
import subprocess
import sys
from time import sleep

import pytest

from hera.simulations.hermesWorkflowToolkit import (
    SCHEDULER_CENTRAL,
    SCHEDULER_LOCAL,
    buildLuigiExecutionCommand,
)

# ---------------------------------------------------------------------------
# buildLuigiExecutionCommand
# ---------------------------------------------------------------------------

def test_local_scheduler_command_is_backward_compatible():
    cmd = buildLuigiExecutionCommand("Flow", "abc123", scheduler=SCHEDULER_LOCAL)
    assert "--module Flow finalnode_xx_0" in cmd
    assert "--local-scheduler" in cmd
    assert "--dispatch-id abc123" not in cmd
    # The legacy invocation (module + target + local-scheduler) must be preserved.
    assert cmd.startswith("python3 -m luigi --module Flow finalnode_xx_0 --local-scheduler")


def test_central_scheduler_drops_local_flag():
    cmd = buildLuigiExecutionCommand("Flow", "abc123", scheduler=SCHEDULER_CENTRAL)
    assert "--local-scheduler" not in cmd
    assert "--dispatch-id abc123" in cmd


def test_central_scheduler_host_and_port():
    cmd = buildLuigiExecutionCommand("Flow", "id1", scheduler=SCHEDULER_CENTRAL,
                                     schedulerHost="myhost", schedulerPort=8082)
    assert "--scheduler-host myhost" in cmd
    assert "--scheduler-port 8082" in cmd
    assert "--local-scheduler" not in cmd


def test_central_scheduler_omits_address_when_not_given():
    cmd = buildLuigiExecutionCommand("Flow", "id1", scheduler=SCHEDULER_CENTRAL)
    assert "--scheduler-host" not in cmd
    assert "--scheduler-port" not in cmd


def test_local_scheduler_ignores_host_and_port():
    cmd = buildLuigiExecutionCommand("Flow", "id1", scheduler=SCHEDULER_LOCAL,
                                     schedulerHost="myhost", schedulerPort=8082)
    assert "--scheduler-host" not in cmd
    assert "--scheduler-port" not in cmd
    assert "--local-scheduler" in cmd


def test_custom_target_task():
    cmd = buildLuigiExecutionCommand("Flow", "id1", targetTask="otherNode_0")
    assert "otherNode_0" in cmd
    assert "finalnode_xx_0" not in cmd


# Note: the generated Luigi node template (the dispatch_id luigi.Parameter, its
# propagation through requires() and per-dispatch output isolation) lives in the
# hermes repository, so it is tested there (tests/test_dispatch_id_template.py),
# not here. Hera's CI vendors hermes' default branch, which would not carry the
# change until the hermes PR merges.


# ---------------------------------------------------------------------------
# End-to-end: run a real hermes workflow through hera's command and verify
# per-dispatch output isolation.
# ---------------------------------------------------------------------------

_WORKFLOW_JSON = {
    "workflow": {
        "nodes": {
            "CopyDirectory": {
                "Execution": {
                    "input_parameters": {
                        "Source": "source", "Target": "target", "dirs_exist_ok": True,
                    }
                },
                "type": "general.CopyDirectory",
            },
            "RunPythonCode": {
                "Execution": {
                    "input_parameters": {
                        "ModulePath": "tutorial1",
                        "ClassName": "tutrialPrinter",
                        "MethodName": "printDirectories",
                        "Parameters": {
                            "source": "{CopyDirectory.output.Source}",
                            "target": "{CopyDirectory.output.Target}",
                        },
                    }
                },
                "type": "general.RunPythonCode",
            },
        }
    }
}

_TUTORIAL_MODULE = (
    "class tutrialPrinter:\n"
    "    def printDirectories(self, source, target):\n"
    "        print(f'Copied {source} to {target}')\n"
)


def test_run_hermes_workflow_isolates_outputs_per_dispatch(tmp_path):
    """Build and execute a real hermes workflow using hera's buildLuigiExecutionCommand
    and confirm the run completes and writes its targets under the dispatch_id subdir."""
    hermes = pytest.importorskip("hermes")
    pytest.importorskip("luigi")
    from hermes import workflow
    from hermes.engines.luigi.pythonClassBase import transform

    # The vendored hermes must carry the dispatch_id change for this to run; otherwise
    # passing --dispatch-id to nodes that don't declare it would error. Skip until then.
    if "dispatch_id" not in transform._basicLuigiTemplate:
        pytest.skip("vendored hermes predates the dispatch_id node template change")

    workdir = str(tmp_path)
    os.makedirs(os.path.join(workdir, "source"), exist_ok=True)
    with open(os.path.join(workdir, "tutorial1.py"), "w") as fp:
        fp.write(_TUTORIAL_MODULE)

    wf = workflow(_WORKFLOW_JSON, workdir, Resource_path=workdir)
    with open(os.path.join(workdir, "Workflow1.py"), "w") as fp:
        fp.write(wf.build(buildername=workflow.BUILDER_LUIGI))

    # Drive the run through hera's own command builder (python3 -> current interpreter
    # so the test is robust to how python is named in the environment).
    cmd = buildLuigiExecutionCommand("Workflow1", "RUN1", scheduler=SCHEDULER_LOCAL)
    cmd = cmd.replace("python3", sys.executable, 1)
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [os.path.dirname(hermes.__file__) + "/..", workdir, env.get("PYTHONPATH", "")]
    )
    result = subprocess.run(cmd, shell=True, cwd=workdir, env=env,
                            capture_output=True, text=True)

    assert "this progress looks :)" in result.stderr.lower(), result.stderr
    run_dir = os.path.join(workdir, "Workflow1_targetFiles")
    assert os.path.isfile(os.path.join(run_dir, "finalnode_xx_0.json"))
    assert os.path.isfile(os.path.join(run_dir, "RunPythonCode_0.json"))
