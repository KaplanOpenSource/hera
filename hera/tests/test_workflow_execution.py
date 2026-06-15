"""
Tests for the Luigi execution command construction and the dispatch_id mechanism
introduced for the centralized-scheduler support (issue #918).

Covers:
    - buildLuigiExecutionCommand: local vs central scheduler, host/port, dispatch-id.
    - The generated Luigi node template carries a dispatch_id luigi.Parameter and
      propagates it through requires()/output() (guarded by hermes availability).
"""

import pytest

from hera.simulations.hermesWorkflowToolkit import (
    buildLuigiExecutionCommand,
    SCHEDULER_LOCAL,
    SCHEDULER_CENTRAL,
)


# ---------------------------------------------------------------------------
# buildLuigiExecutionCommand
# ---------------------------------------------------------------------------

def test_local_scheduler_command_is_backward_compatible():
    cmd = buildLuigiExecutionCommand("Flow", "abc123", scheduler=SCHEDULER_LOCAL)
    assert "--module Flow finalnode_xx_0" in cmd
    assert "--local-scheduler" in cmd
    assert "--dispatch-id abc123" in cmd
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


# ---------------------------------------------------------------------------
# Generated Luigi node template (requires hermes to be importable)
# ---------------------------------------------------------------------------

def test_generated_template_declares_and_propagates_dispatch_id():
    pytest.importorskip("hermes")
    from hermes.engines.luigi.pythonClassBase import transform

    template = transform._basicLuigiTemplate
    # The parameter must be declared on every generated node class...
    assert "dispatch_id = luigi.Parameter(default=\"\")" in template
    # ...propagated to required tasks (Luigi does not thread params automatically)...
    assert "dispatch_id=self.dispatch_id" in template
    # ...and used to isolate per-dispatch output targets.
    assert "self.dispatch_id or" in template
