"""Which of a repository's toolkits can actually be built."""
from __future__ import annotations

import traceback
from typing import Any, Dict, List


class ToolkitCheck:
    """The toolkits a repository tree names, tried one at a time.

    hera's loader skips a toolkit it cannot build with an info line and carries
    on, so a tree whose toolkits all fail loads cleanly and adds nothing. Trying
    them here is what turns that silence into a reason.
    """

    def __init__(self, tree: Dict[str, Any], projectName: str):
        self.tree = tree
        self.projectName = projectName

    def problems(self) -> List[str]:
        """One line per toolkit that cannot be built, in the tree's own order.

        Reporting the failures is the whole purpose, so every one is caught and
        none is thrown away: the message is answered and the traceback printed
        to the server console.
        """
        from hera.toolkit import ToolkitHome
        home = ToolkitHome(projectName=self.projectName)
        problems = []
        for toolkitName in self.tree:
            try:
                home.getToolkit(toolkitName=toolkitName)
            except Exception as error:
                problems.append(f"{toolkitName}: {type(error).__name__}: {error}")
                print(f"repository load: toolkit {toolkitName} cannot be built:\n"
                      f"{traceback.format_exc()}")
        return problems
