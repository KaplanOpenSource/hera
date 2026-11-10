from typing import Dict, Optional

from .active_workflow_run_one import ActiveWorkflowRunOne
from .run_status import RunStatus


class ActiveWorkflowRunList:
    """The workflow runs the server knows about, keyed by token. Kept in memory only.

    A new run of a doc replaces its old run, so the old token is no longer found.
    """

    def __init__(self):
        self._runs: Dict[str, ActiveWorkflowRunOne] = {}

    def add(self, run: ActiveWorkflowRunOne) -> None:
        old = self.get_by_doc_id(run.project_name, run.doc_id)
        if old is not None:
            del self._runs[old.token]
        self._runs[run.token] = run

    def get_by_token(self, token: str) -> Optional[ActiveWorkflowRunOne]:
        return self._runs.get(token)

    # For a client that forgot the token, e.g. after a browser reload.
    def get_by_doc_id(self, project_name: str, doc_id: str) -> Optional[ActiveWorkflowRunOne]:
        for run in self._runs.values():
            if run.project_name == project_name and run.doc_id == doc_id:
                return run
        return None

    def any_running(self) -> bool:
        for run in self._runs.values():
            if run.status == RunStatus.RUNNING:
                return True
        return False
