"""The repository files under a folder."""
from __future__ import annotations

import glob
import os
from typing import Any, Dict, List

from repository_api.errors import NoSuchRepository
from repository_api.repo_json_file import RepoJsonFile

# An OpenFOAM case configuration is JSON in the same folders, and never a repository.
CASE_CONFIGURATION = "caseConfiguration.json"


class RepoFolder:
    """One folder, walked when asked.

    The folder is a user setting and may well point at nothing yet, so a folder
    that is not there is empty and not an error.
    """

    def __init__(self, folder: str):
        self.folder = os.path.expanduser(folder)

    def repositories(self) -> Dict[str, Any]:
        """Every repository file under the folder, by path, with its tree."""
        if not os.path.isdir(self.folder):
            return {}
        found = {}
        pattern = os.path.join(self.folder, "**", "*.json")
        for path in sorted(glob.glob(pattern, recursive=True)):
            if os.path.basename(path) == CASE_CONFIGURATION:
                continue
            # Picking the repositories out of a folder of JSON files is this
            # method's job, so the one refusal that means "not one" is the
            # only thing skipped here.
            try:
                found[path] = RepoJsonFile(path).tree()
            except NoSuchRepository:
                continue
        return found

    def files(self) -> List[str]:
        """The repository files under the folder, sorted."""
        return list(self.repositories().keys())
