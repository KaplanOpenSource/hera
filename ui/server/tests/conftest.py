import json
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


class FakeHeraRecord:
    """What the fake hera saw, read back from a file.

    WorkflowRunner runs the workflow in a forked child (#1033), so anything the
    fake records in memory dies with the child. Each call is appended to a file
    instead, which the parent reads after the run.
    """

    def __init__(self, path):
        self._path = path

    def _events(self):
        if not os.path.exists(self._path):
            return []
        with open(self._path) as handle:
            return [json.loads(line) for line in handle if line.strip()]

    def _append(self, event):
        with open(self._path, "a") as handle:
            handle.write(json.dumps(event) + "\n")

    @property
    def get_toolkit_kwargs(self):
        for event in self._events():
            if event["event"] == "getToolkit":
                return {"toolkitName": event["toolkitName"], "projectName": event["projectName"]}
        return None

    @property
    def execute_calls(self):
        return [(e["workflowName"], e["scheduler"]) for e in self._events() if e["event"] == "execute"]

    @property
    def pythonpath_during_run(self):
        """PYTHONPATH as seen inside the child while the workflow was executing."""
        for event in self._events():
            if event["event"] == "execute":
                return event["pythonpath"]
        return None


@pytest.fixture
def install_fake_hera(monkeypatch, tmp_path):
    """Install a fake `hera` module so code under test doesn't need real hera.

    Returns a factory: ``install_fake_hera(files_directory, on_execute=None)``.
    The returned record captures the getToolkit kwargs and every
    executeWorkflowFromDB call so tests can assert on them.
    """
    def _install(files_directory, on_execute=None):
        record = FakeHeraRecord(str(tmp_path / "fake-hera-calls.jsonl"))

        def execute_workflow_from_db(workflow_name, scheduler):
            record._append({
                "event": "execute",
                "workflowName": workflow_name,
                "scheduler": scheduler,
                "pythonpath": os.environ.get("PYTHONPATH"),
            })
            if on_execute is not None:
                return on_execute(workflow_name, scheduler)
            return "dispatch-default"

        toolkit = types.SimpleNamespace(
            FilesDirectory=files_directory,
            executeWorkflowFromDB=execute_workflow_from_db,
        )

        def get_toolkit(toolkitName, projectName):
            record._append({"event": "getToolkit", "toolkitName": toolkitName, "projectName": projectName})
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


@pytest.fixture
def warmed_up(monkeypatch):
    """Pretend hera finished importing.

    /exec answers WARMING_UP until the background warmup thread imported hera
    (#1030); that thread only starts in the app's lifespan, which unit tests
    don't run.
    """
    import server

    monkeypatch.setattr(server, "warmup", types.SimpleNamespace(ready=True))
