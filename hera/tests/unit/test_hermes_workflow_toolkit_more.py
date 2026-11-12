"""hermesWorkflowToolkit.executeWorkflowFromDB -- the last uncovered method.

``test_hermes_workflow_toolkit_db.py`` covers the whole DB-backed lifecycle
and ``buildLuigiExecutionCommand``, and states that ``executeWorkflowFromDB``
is deliberately not covered because "it shells out to ``python3 -m luigi``
via subprocess and writes a generated module to disk".  Both halves of that
are seams rather than obstacles: the module's own ``subprocess`` is
monkeypatched, and the generated module lands under the toolkit's
per-test ``FilesDirectory`` (the ``unit_files_directory`` fixture), i.e.
inside ``tmp_path``.

Everything else in the method is real: the workflows are real
``hermes.workflow`` objects built from the five nodes
``workflow_Eulerian._requiredNodeList`` demands, stored in and read back out
of the in-memory datalayer, and really put through ``hermes.build()``.  What
is asserted is what crosses the subprocess seam, what is on disk while the
command runs, what survives afterwards, and what the method returns -- never
anything a stub invented.

Bug fixed here:

* B294 (fixed) ``executeWorkflowFromDB`` generated (or accepted) a
  ``dispatch_id``, logged it and returned it, but built the Luigi command
  line without it, so on a central scheduler every run carried the literal
  ``--dispatch-id None`` and concurrent runs could not be told apart.
"""
import os
import subprocess
import types

import pytest

import hera.simulations.hermesWorkflowToolkit as workflowToolkitModule
from hera.simulations.hermesWorkflowToolkit import SCHEDULER_CENTRAL, SCHEDULER_LOCAL

_REQUIRED_NODES = ["controlDict", "fvSolution", "fvSchemes", "fileWriter", "Parameters"]


def _workflowJSON(alpha=1.0, solver="simpleFoam"):
    """A minimal but real hermes workflow; only Parameters.alpha varies."""
    nodes = {nodeName: {"type": "general.Parameters", "Execution": {"input_parameters": {}}}
             for nodeName in _REQUIRED_NODES}
    nodes["Parameters"]["Execution"]["input_parameters"] = {"alpha": alpha}
    return {"workflow": {"nodeList": list(_REQUIRED_NODES), "nodes": nodes, "solver": solver}}


@pytest.fixture()
def tk(unit_toolkit_factory):
    from hera import toolkitHome

    return unit_toolkit_factory(toolkitHome.SIMULATIONS_WORKFLOWS)


class _LuigiSeam(list):
    """The recorded `subprocess.run` calls, plus a disk snapshot per call.

    A list subclass so the tests can index it directly while it still
    carries the directory to snapshot.
    """

    filesDirectory = None
    filesAtRunTime = None

    def record(self, command, stdout=None, stderr=None, check=None, **kwargs):
        self.append(dict(command=command, stdout=stdout, stderr=stderr, check=check))
        if self.filesDirectory is not None:
            self.filesAtRunTime = sorted(os.listdir(self.filesDirectory))
        return types.SimpleNamespace(returncode=0)


@pytest.fixture()
def luigiSeam(monkeypatch):
    """Record every `subprocess.run`, plus what was on disk at the time."""
    seam = _LuigiSeam()
    monkeypatch.setattr(workflowToolkitModule.subprocess, "run", seam.record)
    return seam


@pytest.fixture()
def oneWorkflow(tk, luigiSeam):
    """A single workflow in the DB; returns its name."""
    document = tk.addWorkflowToGroup(_workflowJSON(1.0), "flow")
    luigiSeam.filesDirectory = tk.FilesDirectory
    return document.desc["workflowName"]


# ---------------------------------------------------------------------------
# The Luigi command line
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestTheLuigiInvocation:
    def test_one_luigi_process_is_started_per_workflow(self, tk, oneWorkflow, luigiSeam):
        tk.executeWorkflowFromDB(oneWorkflow)
        assert len(luigiSeam) == 1

    def test_the_generated_module_is_the_one_luigi_is_pointed_at(self, tk, oneWorkflow,
                                                                 luigiSeam):
        tk.executeWorkflowFromDB(oneWorkflow)
        command = luigiSeam[0]["command"]
        assert command[:4] == ["python3", "-m", "luigi", "--module"]
        assert command[4] == oneWorkflow

    def test_the_terminal_task_triggers_the_whole_dag(self, tk, oneWorkflow, luigiSeam):
        tk.executeWorkflowFromDB(oneWorkflow)
        assert "finalnode_xx_0" in luigiSeam[0]["command"]

    def test_the_local_scheduler_is_the_default(self, tk, oneWorkflow, luigiSeam):
        tk.executeWorkflowFromDB(oneWorkflow)
        assert "--local-scheduler" in luigiSeam[0]["command"]

    def test_the_central_scheduler_drops_the_local_flag(self, tk, oneWorkflow, luigiSeam):
        tk.executeWorkflowFromDB(oneWorkflow, scheduler=SCHEDULER_CENTRAL)
        assert "--local-scheduler" not in luigiSeam[0]["command"]

    def test_a_central_host_and_port_reach_the_command_line(self, tk, oneWorkflow, luigiSeam):
        tk.executeWorkflowFromDB(oneWorkflow, scheduler=SCHEDULER_CENTRAL,
                                 schedulerHost="luigi.example", schedulerPort=9999)
        command = luigiSeam[0]["command"]
        assert command[command.index("--scheduler-host") + 1] == "luigi.example"
        assert command[command.index("--scheduler-port") + 1] == "9999"

    def test_the_command_is_split_into_a_list_so_no_shell_is_involved(self, tk, oneWorkflow,
                                                                     luigiSeam):
        assert isinstance(luigiSeam, list)
        tk.executeWorkflowFromDB(oneWorkflow)
        assert all(isinstance(token, str) for token in luigiSeam[0]["command"])

    def test_a_failing_luigi_run_is_not_swallowed(self, tk, oneWorkflow, monkeypatch):
        def failing(command, **kwargs):
            raise subprocess.CalledProcessError(returncode=1, cmd=command)

        monkeypatch.setattr(workflowToolkitModule.subprocess, "run", failing)
        with pytest.raises(subprocess.CalledProcessError):
            tk.executeWorkflowFromDB(oneWorkflow)

    def test_check_is_requested_so_a_nonzero_exit_raises(self, tk, oneWorkflow, luigiSeam):
        tk.executeWorkflowFromDB(oneWorkflow)
        assert luigiSeam[0]["check"] is True

    def test_the_streams_the_caller_gave_are_the_ones_luigi_writes_to(self, tk, oneWorkflow,
                                                                     luigiSeam, tmp_path):
        with open(tmp_path / "out.log", "w") as outputStream:
            tk.executeWorkflowFromDB(oneWorkflow, stdout=outputStream, stderr=outputStream)
        assert luigiSeam[0]["stdout"] is luigiSeam[0]["stderr"] is outputStream

    def test_by_default_the_streams_are_inherited(self, tk, oneWorkflow, luigiSeam):
        tk.executeWorkflowFromDB(oneWorkflow)
        assert luigiSeam[0]["stdout"] is None
        assert luigiSeam[0]["stderr"] is None


# ---------------------------------------------------------------------------
# dispatch_id -- B294
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestTheDispatchId:
    def test_an_omitted_dispatch_id_is_generated(self, tk, oneWorkflow, luigiSeam):
        dispatchId = tk.executeWorkflowFromDB(oneWorkflow)
        assert isinstance(dispatchId, str) and len(dispatchId) == 32
        int(dispatchId, 16)  # a uuid4 hex

    def test_two_runs_get_two_different_generated_ids(self, tk, oneWorkflow, luigiSeam):
        first = tk.executeWorkflowFromDB(oneWorkflow)
        second = tk.executeWorkflowFromDB(oneWorkflow)
        assert first != second

    def test_an_explicit_dispatch_id_is_the_one_returned(self, tk, oneWorkflow, luigiSeam):
        assert tk.executeWorkflowFromDB(oneWorkflow, dispatch_id="deadbeef") == "deadbeef"

    def test_the_local_scheduler_never_receives_a_dispatch_id(self, tk, oneWorkflow,
                                                              luigiSeam):
        """Correct today, and the reason the defect below is invisible locally."""
        tk.executeWorkflowFromDB(oneWorkflow, dispatch_id="deadbeef")
        assert "--dispatch-id" not in luigiSeam[0]["command"]

    def test_the_central_scheduler_should_receive_the_dispatch_id(self, tk, oneWorkflow,
                                                                  luigiSeam):
        tk.executeWorkflowFromDB(oneWorkflow, dispatch_id="deadbeef",
                                 scheduler=SCHEDULER_CENTRAL)
        command = luigiSeam[0]["command"]
        assert command[command.index("--dispatch-id") + 1] == "deadbeef"

    def test_two_concurrent_runs_are_told_apart_on_the_command_line(self, tk,
                                                                     oneWorkflow,
                                                                     luigiSeam):
        """B294: each run's own id reaches the central scheduler."""
        first = tk.executeWorkflowFromDB(oneWorkflow, scheduler=SCHEDULER_CENTRAL)
        second = tk.executeWorkflowFromDB(oneWorkflow, scheduler=SCHEDULER_CENTRAL)
        assert first != second
        assert luigiSeam[0]["command"] != luigiSeam[1]["command"]


# ---------------------------------------------------------------------------
# What lands on disk
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestWhatIsWrittenAndCleanedUp:
    def test_the_generated_python_module_exists_while_luigi_runs(self, tk, oneWorkflow,
                                                                 luigiSeam):
        tk.executeWorkflowFromDB(oneWorkflow)
        assert f"{oneWorkflow}.py" in luigiSeam.filesAtRunTime

    def test_the_generated_python_module_is_removed_afterwards(self, tk, oneWorkflow,
                                                               luigiSeam):
        tk.executeWorkflowFromDB(oneWorkflow)
        assert not os.path.exists(os.path.join(tk.FilesDirectory, f"{oneWorkflow}.py"))

    def test_the_workflow_json_is_written_and_kept(self, tk, oneWorkflow, luigiSeam):
        tk.executeWorkflowFromDB(oneWorkflow)
        assert os.path.isfile(os.path.join(tk.FilesDirectory, f"{oneWorkflow}.json"))

    def test_the_previous_target_files_are_removed_before_the_run(self, tk, oneWorkflow,
                                                                  luigiSeam):
        """Luigi treats a target file as "this task is done", so a re-run has
        to start from a clean directory."""
        stale = os.path.join(tk.FilesDirectory, f"{oneWorkflow}_targetFiles")
        os.makedirs(stale)
        with open(os.path.join(stale, "finalnode_xx_0"), "w") as staleTarget:
            staleTarget.write("done")
        tk.executeWorkflowFromDB(oneWorkflow)
        assert not os.path.exists(stale)

    def test_a_missing_target_file_directory_is_not_an_error(self, tk, oneWorkflow,
                                                             luigiSeam):
        assert not os.path.exists(os.path.join(tk.FilesDirectory,
                                               f"{oneWorkflow}_targetFiles"))
        tk.executeWorkflowFromDB(oneWorkflow)
        assert len(luigiSeam) == 1


# ---------------------------------------------------------------------------
# Selecting what to execute
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestWhichWorkflowsAreExecuted:
    def test_a_name_that_is_in_no_document_runs_nothing_and_still_returns_an_id(
            self, tk, luigiSeam):
        dispatchId = tk.executeWorkflowFromDB("noSuchWorkflow")
        assert luigiSeam == []
        assert isinstance(dispatchId, str)

    def test_a_list_of_names_executes_each_of_them(self, tk, luigiSeam):
        first = tk.addWorkflowToGroup(_workflowJSON(1.0), "flow").desc["workflowName"]
        second = tk.addWorkflowToGroup(_workflowJSON(2.0), "flow").desc["workflowName"]
        luigiSeam.filesDirectory = tk.FilesDirectory
        tk.executeWorkflowFromDB([first, second])
        modules = [run["command"][4] for run in luigiSeam]
        assert sorted(modules) == sorted([first, second])

    def test_every_workflow_of_a_list_shares_the_one_dispatch_id(self, tk, luigiSeam):
        first = tk.addWorkflowToGroup(_workflowJSON(1.0), "flow").desc["workflowName"]
        second = tk.addWorkflowToGroup(_workflowJSON(2.0), "flow").desc["workflowName"]
        luigiSeam.filesDirectory = tk.FilesDirectory
        dispatchId = tk.executeWorkflowFromDB([first, second], dispatch_id="shared")
        assert dispatchId == "shared"
        assert len(luigiSeam) == 2

    def test_the_workflow_can_be_identified_by_its_resource(self, tk, oneWorkflow,
                                                            luigiSeam):
        resource = tk.getWorkflowDocumentFromDB(oneWorkflow)[0].resource
        tk.executeWorkflowFromDB(resource)
        assert [run["command"][4] for run in luigiSeam] == [oneWorkflow]

    def test_the_scheduler_choice_applies_to_every_workflow_of_a_list(self, tk, luigiSeam):
        first = tk.addWorkflowToGroup(_workflowJSON(1.0), "flow").desc["workflowName"]
        second = tk.addWorkflowToGroup(_workflowJSON(2.0), "flow").desc["workflowName"]
        luigiSeam.filesDirectory = tk.FilesDirectory
        tk.executeWorkflowFromDB([first, second], scheduler=SCHEDULER_LOCAL)
        assert all("--local-scheduler" in run["command"] for run in luigiSeam)
