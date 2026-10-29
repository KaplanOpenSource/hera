"""The files under one project's own directory."""
from __future__ import annotations

import os
import shutil
from typing import Any, Dict, List


class ProjectFiles:
    """One project's files directory, and what inside it may be removed.

    Lives for one request. An empty directory means the project has none, so
    nothing is ever ours to remove.
    """

    # Written at project creation, so no document points at it.
    CASE_CONFIG = "caseConfiguration.json"

    def __init__(self, directory: str = ""):
        self.directory = ""
        if directory:
            self.directory = os.path.abspath(directory)

    def holds(self, path: str) -> bool:
        """True when the path sits in this directory, so it is ours to remove."""
        if not self.directory or not path:
            return False
        target = os.path.abspath(path)
        return target == self.directory or target.startswith(self.directory + os.sep)

    def list(self, relPath: str = "") -> List[Dict[str, Any]]:
        """The files and folders right inside relPath, or none when it is not a folder in here."""
        if not self.directory or not os.path.isdir(self.directory):
            return []
        root = os.path.realpath(self.directory)
        target = os.path.realpath(os.path.join(root, relPath))
        if target != root and not target.startswith(root + os.sep):
            return []
        if not os.path.isdir(target):
            return []
        entries = [{"name": e.name, "isFolder": e.is_dir()} for e in os.scandir(target)]
        entries.sort(key=lambda e: (not e["isFolder"], e["name"].lower()))
        return entries

    def remove(self, path: str) -> bool:
        """Remove the path if this directory holds it. Answers whether it went."""
        if not self.holds(path):
            return False
        target = os.path.abspath(path)
        if os.path.isdir(target):
            shutil.rmtree(target)
            return True
        if os.path.exists(target):
            os.remove(target)
            return True
        return False

    def remove_if_empty(self) -> bool:
        """Drop the directory itself, once nothing but our own leftovers is in it."""
        if not self.directory or not os.path.isdir(self.directory):
            return False
        caseConfig = os.path.join(self.directory, self.CASE_CONFIG)
        if os.path.exists(caseConfig):
            os.remove(caseConfig)
        if os.listdir(self.directory):
            return False
        os.rmdir(self.directory)
        return True
