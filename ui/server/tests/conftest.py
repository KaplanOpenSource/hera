import os
import sys
import types

import pytest

# Make the server modules (server.py, api_models.py, ...) importable regardless of
# where pytest is run from.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# server.py calls argparse.parse_args() at import time; give it clean args so it
# doesn't try to parse pytest's argv.
sys.argv = ["server"]


@pytest.fixture
def install_fake_hera(monkeypatch):
    """Install fakes so the runner works without real hera or luigi.

    Returns a factory: ``install_fake_hera(files_directory, on_execute=None)``.

    It fakes two things the workflow child imports:
      - ``hera.toolkitHome`` -> a toolkit exposing ``FilesDirectory``.
      - ``execute_workflow_inprocess.executeWorkflowFromDB_inprocess`` -> a stub that
        calls ``on_execute(workflow_name)`` and returns a dispatch id. This avoids
        importing luigi and the real build/run.

    The child runs in a forked process, so anything ``on_execute`` records in memory
    (or writes to the recorded namespace) is NOT visible to the parent. Tests must
    assert on the run's observable result instead: the returned output / dispatch id
    / chunks, or what ``poll`` reports. ``on_execute`` can ``os.write(1/2, ...)`` to
    produce captured output, and return a dispatch id.
    """
    def _install(files_directory, on_execute=None):
        toolkit = types.SimpleNamespace(FilesDirectory=files_directory)

        def get_toolkit(toolkitName, projectName):
            return toolkit

        toolkit_home = types.SimpleNamespace(
            SIMULATIONS_WORKFLOWS="SIMULATIONS_WORKFLOWS",
            getToolkit=get_toolkit,
        )
        fake_hera = types.ModuleType("hera")
        fake_hera.toolkitHome = toolkit_home
        monkeypatch.setitem(sys.modules, "hera", fake_hera)

        def execute_inprocess(workflow_toolkit, workflow_name, workers=1):
            if on_execute is not None:
                return on_execute(workflow_name)
            return "dispatch-default"

        fake_exec = types.ModuleType("execute_workflow_inprocess")
        fake_exec.executeWorkflowFromDB_inprocess = execute_inprocess
        monkeypatch.setitem(sys.modules, "execute_workflow_inprocess", fake_exec)

    return _install
