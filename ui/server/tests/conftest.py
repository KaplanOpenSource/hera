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
    """Install a fake `hera` module so code under test doesn't need real hera.

    Returns a factory: ``install_fake_hera(files_directory, on_execute=None)``.
    The returned record captures the getToolkit kwargs and every
    executeWorkflowFromDB call so tests can assert on them.
    """
    def _install(files_directory, on_execute=None):
        record = types.SimpleNamespace(get_toolkit_kwargs=None, execute_calls=[])

        def execute_workflow_from_db(workflow_name, scheduler):
            record.execute_calls.append((workflow_name, scheduler))
            if on_execute is not None:
                return on_execute(workflow_name, scheduler)
            return "dispatch-default"

        toolkit = types.SimpleNamespace(
            FilesDirectory=files_directory,
            executeWorkflowFromDB=execute_workflow_from_db,
        )

        def get_toolkit(toolkitName, projectName):
            record.get_toolkit_kwargs = {"toolkitName": toolkitName, "projectName": projectName}
            return toolkit

        toolkit_home = types.SimpleNamespace(
            SIMULATIONS_WORKFLOWS="SIMULATIONS_WORKFLOWS",
            getToolkit=get_toolkit,
        )
        fake_hera = types.ModuleType("hera")
        fake_hera.toolkitHome = toolkit_home
        monkeypatch.setitem(sys.modules, "hera", fake_hera)
        return record

    return _install
