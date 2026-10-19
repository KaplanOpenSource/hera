"""A repository tree, loaded into a project."""
from __future__ import annotations

import os
from typing import Any, Dict

# A repository's top level is toolkit names; under each one, only these sections.
SECTIONS = {"config", "datasource", "measurements", "simulations", "cache", "function"}


class RepositoryJson:
    """One repository tree and the directory its relative paths are against.

    Built for one request and dropped: every call makes a fresh dataToolkit.
    """

    def __init__(self, tree: Dict[str, Any], baseDir: str = ""):
        self.tree = tree
        self.baseDir = os.path.abspath(os.path.expanduser(baseDir))

    @staticmethod
    def is_repository(tree: Any) -> bool:
        """True when the tree is shaped like a repository file."""
        if not isinstance(tree, dict) or not tree:
            return False
        for toolkitValue in tree.values():
            if not isinstance(toolkitValue, dict):
                return False
            for sectionKey in toolkitValue:
                if sectionKey.lower() not in SECTIONS:
                    return False
        return True

    @staticmethod
    def why_not_a_repository(tree: Any) -> str:
        """Why ``is_repository`` said no, in a line, so a refusal can explain itself."""
        if not isinstance(tree, dict):
            return f"the top level is {type(tree).__name__}, not an object of toolkit names"
        if not tree:
            return "it is empty"
        for toolkitName, toolkitValue in tree.items():
            if not isinstance(toolkitValue, dict):
                return f"toolkit {toolkitName} is {type(toolkitValue).__name__}, not an object"
            for sectionKey in toolkitValue:
                if sectionKey.lower() not in SECTIONS:
                    return (f"{toolkitName}.{sectionKey} is not a section: "
                            f"expected one of {', '.join(sorted(SECTIONS))}")
        return ""

    def resolved(self) -> Dict[str, Any]:
        """The same tree with its relative datasource paths made absolute."""
        from hera.utils.data.toolkit import dataToolkit
        return dataToolkit.resolveDataSourcePaths(self.tree, self.baseDir)

    def load_into(self, projectName: str, overwrite: bool = False) -> None:
        """Add every datasource in the tree to the project."""
        from hera.utils.data.toolkit import dataToolkit
        dataToolkit().loadAllDatasourcesInRepositoryJSONToProject(
            projectName=projectName,
            repositoryJSON=self.tree,
            basedir=self.baseDir,
            overwrite=overwrite)
