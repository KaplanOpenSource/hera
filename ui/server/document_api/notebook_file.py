"""The .ipynb file a notebook document points at."""
from __future__ import annotations

import json
import os
from typing import Any, Dict, List

NOTEBOOK_FORMAT = 4
NOTEBOOK_FORMAT_MINOR = 5


class NotebookFile:
    """One notebook path on disk.

    Writing only ever creates a missing file: a notebook the user already has is
    never touched, so pointing a new document at an existing notebook keeps it.
    """

    def __init__(self, path: str):
        self.path = os.path.expanduser(path)

    def exists(self) -> bool:
        return os.path.exists(self.path)

    def write_if_missing(self, projectName: str, toolkitNames: List[str]) -> bool:
        """Write an empty notebook here. Answers False when a file is already there."""
        if self.exists():
            return False
        parent = os.path.dirname(self.path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        with open(self.path, "w") as handle:
            handle.write(json.dumps(self._empty(projectName, toolkitNames), indent=2))
        return True

    @staticmethod
    def _first_cell_source(projectName: str, toolkitNames: List[str]) -> List[str]:
        """The opening cell: the project, and one line per toolkit it uses."""
        lines = [
            "import hera\n",
            f'PROJECT_NAME = "{projectName}"\n',
            "\n",
            "p = hera.datalayer.Project(projectName=PROJECT_NAME)\n",
        ]
        if toolkitNames:
            lines.append("\n")
        for name in toolkitNames:
            # The toolkit name is not always a legal python name.
            varName = "".join(c if c.isalnum() or c == "_" else "_" for c in name)
            lines.append(
                f'{varName}_tk = hera.toolkit.ToolkitHome()'
                f'.getToolkit("{name}", projectName=PROJECT_NAME)\n')
        lines.append("p.getMetadata()\n")
        return lines

    def _empty(self, projectName: str, toolkitNames: List[str]) -> Dict[str, Any]:
        return {
            "nbformat": NOTEBOOK_FORMAT,
            "nbformat_minor": NOTEBOOK_FORMAT_MINOR,
            "metadata": {"kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            }},
            "cells": [{
                "cell_type": "code",
                "metadata": {},
                "source": self._first_cell_source(projectName, toolkitNames),
                "outputs": [],
                "execution_count": None,
            }],
        }
