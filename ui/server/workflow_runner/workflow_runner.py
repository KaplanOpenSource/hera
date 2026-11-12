import threading
from typing import Any, Dict

from .active_workflow_run_list import ActiveWorkflowRunList
from .active_workflow_run_one import ActiveWorkflowRunOne
from .run_status import RunStatus
from time_utils import now_readable


class WorkflowRunner:
    """What the routes call to start, poll and stop workflow runs.

    ``start`` runs the workflow in the background and returns a token right away;
    the caller polls ``poll(token)`` until the run is done and then reads the whole
    output. One run at a time: ``start`` reports busy while a run is in progress.
    """

    def __init__(self):
        # Serialize starts: the local Luigi scheduler / DB access is not meant to run concurrently.
        self._lock = threading.Lock()
        self._runs = ActiveWorkflowRunList()

    def start(self, project_name: str, doc: Dict[str, Any]) -> dict:
        """Start a run in the background. Returns ``{"token"}`` or ``{"status": "busy"}``."""
        with self._lock:
            if self._runs.any_running():
                return {"status": RunStatus.BUSY}
            run = ActiveWorkflowRunOne(project_name, doc)
            self._runs.add(run)
            run.start()
            return {"token": run.token}

    def poll(self, token: str) -> dict:
        """Return ``{"status", "error", "chunks"}`` for a token, or not_found."""
        run = self._runs.get_by_token(token)
        if run is None:
            print(f"[{now_readable()}] poll {token}: not_found")
            return {"status": RunStatus.NOT_FOUND, "error": "", "chunks": None}
        print(f"[{now_readable()}] poll {token}: {run.status}")
        return run.poll()

    def stop(self, token: str) -> dict:
        """Stop a running run. Returns its poll result, or not_found."""
        run = self._runs.get_by_token(token)
        if run is None:
            return {"status": RunStatus.NOT_FOUND, "error": "", "chunks": None}
        if run.status == RunStatus.RUNNING:
            print(f"[{now_readable()}] stop {token}: pid {run.pid}")
            run.stop()
        return run.poll()
