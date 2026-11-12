"""The registered repositories: the rows the default project keeps."""
from __future__ import annotations

from typing import Any, Dict, List


class RepositoryTable:
    """Every repository registered in the default project.

    Built for one request and dropped: every call makes a fresh dataToolkit.
    """

    def rows(self) -> List[Dict[str, Any]]:
        """One row per registered repository, with its name and its file."""
        from hera.utils.data.toolkit import dataToolkit
        table = dataToolkit().getRepositoryTable()
        if table.empty:
            return []
        return table.to_dict(orient="records")

    def load_all_into(self, projectName: str, overwrite: bool = False) -> None:
        """Add the datasources of every registered repository to the project."""
        from hera.utils.data.toolkit import dataToolkit
        dataToolkit().loadAllDatasourcesInAllRepositoriesToProject(
            projectName=projectName, overwrite=overwrite)
