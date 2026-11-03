import json
import os
from typing import Any, Dict, List, Optional

# The fields of a saved run that the run list shows; the rest is only read with the run.
SUMMARY_FIELDS = ("id", "token", "startTime", "endTime", "status")


class WorkflowRunHistory:
    """The saved runs of one workflow: one JSON file per finished run. Never deleted."""

    def __init__(self, project_name: str, workflow_name: str):
        self.project_name = project_name
        self.workflow_name = workflow_name

    @staticmethod
    def files_directory(project_name: str) -> str:
        """The project's files directory, read the same way the workflow run reads it."""
        from hera import toolkitHome

        workflow_toolkit = toolkitHome.getToolkit(
            toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS,
            projectName=project_name,
        )
        return workflow_toolkit.FilesDirectory

    def folder(self) -> str:
        return os.path.join(self.files_directory(self.project_name), "workflow_runs", self.workflow_name)

    def save(self, run_id: str, record: Dict[str, Any]) -> None:
        """Write one run as ``<run_id>.json``."""
        os.makedirs(self.folder(), exist_ok=True)
        with open(os.path.join(self.folder(), f"{run_id}.json"), "w") as handle:
            json.dump(record, handle)

    def summaries(self) -> List[Dict[str, Any]]:
        """Every saved run, newest first, without its output."""
        folder = self.folder()
        if not os.path.isdir(folder):
            return []
        result = []
        for name in os.listdir(folder):
            if not name.endswith(".json"):
                continue
            record = self.read(name[: -len(".json")])
            result.append({field: record.get(field, "") for field in SUMMARY_FIELDS})
        result.sort(key=lambda summary: summary["startTime"], reverse=True)
        return result

    def read(self, run_id: str) -> Optional[Dict[str, Any]]:
        """One saved run with its output, or None when there is no such run."""
        # A run id is a plain file name; refuse anything that walks out of the folder.
        if os.path.basename(run_id) != run_id:
            return None
        path = os.path.join(self.folder(), f"{run_id}.json")
        if not os.path.isfile(path):
            return None
        with open(path) as handle:
            record = json.load(handle)
        record["id"] = run_id
        return record
