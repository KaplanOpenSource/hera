"""One project, as the named project routes see it."""
from __future__ import annotations

import glob
import os
from typing import Any, Dict, List

from project_api.project_files import ProjectFiles

# Shown first in the project list; the rest are sorted by name.
DEFAULT_PROJECT = "defaultProject"


class ProjectApi:
    """What the UI can ask about one project, and do to it.

    Holds nothing but the project's name: it is built for one request and
    dropped, so every call reads the database again.
    """

    def __init__(self, projectName: str):
        self.projectName = projectName

    @staticmethod
    def names() -> List[Dict[str, str]]:
        """Every project in the database: the default project first, the rest by name."""
        from hera.datalayer.project import getProjectList
        names = getProjectList()
        first = [{"name": n} for n in names if n == DEFAULT_PROJECT]
        rest = [{"name": n} for n in names if n != DEFAULT_PROJECT]
        rest.sort(key=lambda p: p["name"].lower())
        return first + rest

    def default_files_directory(self) -> str:
        """Where the project's files go when the user names no directory."""
        return os.path.join(os.getcwd(), "projects", self.projectName)

    def create(self, filesDirectory: str = "", loadRepositories: bool = True) -> str:
        """Create the project in the database, and answer its files directory."""
        from types import SimpleNamespace

        from hera.utils.data.CLI import project_create
        from hera.datalayer.project import Project

        directory = filesDirectory or self.default_files_directory()
        project_create(SimpleNamespace(
            projectName=self.projectName,
            directory=directory,
            loadRepositories=loadRepositories,
            overwrite=False))
        # Creates a config document in MongoDB so the project appears in getProjectList()
        Project(projectName=self.projectName, filesDirectory=directory)
        return directory

    def exists(self) -> bool:
        """True when the database holds a project by this name."""
        from hera.datalayer.project import getProjectList
        return self.projectName in getProjectList()

    def files_directory(self) -> str:
        """The directory the project keeps its own files in, or "" if it has none."""
        from hera.datalayer.project import Project
        return Project(projectName=self.projectName).filesDirectory or ""

    def documents(self) -> Dict[str, Any]:
        """The project's documents from all three collections, with their ids."""
        from hera.datalayer import All
        docs = All.getDocumentsAsDict(self.projectName, with_id=True)
        return {"name": self.projectName, "documents": docs["documents"]}

    def delete(self, deleteFiles: bool = False) -> None:
        """Remove every document of the project, and its own files when asked."""
        from hera.datalayer import All

        files = ProjectFiles()
        if deleteFiles:
            files = ProjectFiles(self.files_directory())

        docs = All.getDocumentsAsDict(self.projectName, with_id=True)["documents"]
        # The config document goes last: while it exists, the project still exists.
        docs = sorted(docs, key=lambda d: d["type"] == f"{self.projectName}__config__")
        for doc in docs:
            resource = doc.get("resource")
            if isinstance(resource, str):
                files.remove(resource)
            All.deleteDocumentByID(doc["_id"]["$oid"])

        files.remove_if_empty()

    def detect_notebooks(self, filesDirectory: str) -> int:
        """Register every notebook under the directory that is not tracked yet.

        Answers how many were added. Jupyter's .ipynb_checkpoints copies are skipped.
        """
        from hera.datalayer import Cache_Collection
        filesDir = os.path.expanduser(filesDirectory)
        existing = set()
        for doc in Cache_Collection().getDocuments(projectName=self.projectName, type="notebook"):
            if doc.resource:
                existing.add(os.path.abspath(os.path.expanduser(doc.resource)))
        found = glob.glob(os.path.join(filesDir, "**", "*.ipynb"), recursive=True)
        found = [f for f in found if ".ipynb_checkpoints" not in f]
        addedCount = 0
        for path in found:
            if os.path.abspath(path) in existing:
                continue
            Cache_Collection().addDocument(
                projectName=self.projectName,
                resource=path,
                dataFormat="JSON_dict",
                type="notebook",
                desc={"datasourceName": os.path.splitext(os.path.basename(path))[0]},
            )
            addedCount += 1
        return addedCount
