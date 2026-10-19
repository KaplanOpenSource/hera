"""One repository file on disk."""
from __future__ import annotations

import json
import os
from typing import Any, Dict

from repository_api.errors import NoSuchRepository
from repository_api.repo_json import RepositoryJson


class RepoJsonFile:
    """One path, read when asked.

    Reading is the only place the shape rule is applied, so reading one file and
    scanning a folder cannot disagree about what a repository is. Anything that
    goes wrong keeps the reason the operating system or the JSON parser gave.
    """

    def __init__(self, path: str):
        self.path = os.path.expanduser(path)

    def directory(self) -> str:
        """The folder the file is in: what its relative paths are against."""
        return os.path.dirname(os.path.abspath(self.path))

    def tree(self) -> Dict[str, Any]:
        """The file's JSON, when the file is a repository.

        Raises ``NoSuchRepository`` otherwise, saying what went wrong.
        """
        # The reason is the whole point of catching: it goes into the message.
        try:
            with open(self.path, "r") as handle:
                tree = json.load(handle)
        except OSError as error:
            raise NoSuchRepository(f"cannot read {self.path}: {error}") from error
        except ValueError as error:
            raise NoSuchRepository(f"{self.path} is not JSON: {error}") from error
        if not RepositoryJson.is_repository(tree):
            raise NoSuchRepository(
                f"{self.path} is JSON but not shaped like a repository: "
                f"{RepositoryJson.why_not_a_repository(tree)}")
        return tree

    def repository(self) -> RepositoryJson:
        """The file's tree, against the file's own folder."""
        return RepositoryJson(self.tree(), self.directory())

    def load_into(self, projectName: str, overwrite: bool = False) -> None:
        """Add every datasource in the file to the project."""
        self.repository().load_into(projectName, overwrite)
