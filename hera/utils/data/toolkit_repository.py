import os
import inspect
from typing import Optional, Tuple
import pandas as pd

from hera.datalayer import Project


class ToolkitRepository:
    """
    Repository helper for registering and discovering toolkits as data-sources.

    Conventions for a dynamic toolkit record:
      - Saved as a *Measurements document* with:
          type="ToolkitDataSource", dataFormat="Class"
      - 'resource'  -> a directory to be put on sys.path (usually the *package* dir)
      - 'desc'      -> {
            "toolkit": <datasource_name>,
            "datasourceName": <datasource_name>,
            "version": [major, minor, patch],
            "classpath": "pkg.mod.ClassName",
            "parameters": {...}              # kwargs for the class __init__
        }
    """

    def __init__(self, project_name: str):
        """
        Initialize a ToolkitRepository for the given project.

        Parameters
        ----------
        project_name : str
            The project name to manage toolkit registrations for.
        """
        self._project = Project(projectName=project_name)

    # ---------------------------
    # Registration
    # ---------------------------
    def registerToolkit(
        self,
        *,
        toolkitclass=None,
        datasource_name: str,
        params: Optional[dict] = None,
        version: Tuple[int, int, int] = (0, 0, 1),
        overwrite: bool = False,
        resource: Optional[str] = None,
        classpath: Optional[str] = None,
    ):
        """
        Register a toolkit class as a ToolkitDataSource measurements-document.

        You may pass either:
        - toolkitclass (actual class object) and let this method derive 'resource' and 'classpath'
        - or pass 'resource' + 'classpath' explicitly.

        Returns: the created document (mongoengine object)
        """
        if toolkitclass is None and (resource is None or classpath is None):
            raise TypeError(
                "registerToolkit requires either 'toolkitclass' OR both 'resource' and 'classpath'."
            )

        # Derive resource & classpath from the class object when provided
        if toolkitclass is not None:
            # Compute classpath
            if classpath is None:
                classpath = f"{toolkitclass.__module__}.{toolkitclass.__name__}"

            # Compute resource (directory that contains the package/module)
            mod = inspect.getmodule(toolkitclass)
            if mod is None or not hasattr(mod, "__file__"):
                raise ValueError("Cannot locate module for the provided class.")
            file_path = inspect.getfile(toolkitclass)  # e.g. /tmp/.../demo_pkg/mytoolkit.py
            pkg_dir = os.path.dirname(file_path)       # .../demo_pkg
            resource = resource or pkg_dir

        params = params or {}
        version_list = list(version) if isinstance(version, tuple) else version

        desc = {
            "toolkit": datasource_name,
            "datasourceName": datasource_name,
            "version": version_list,
            "classpath": classpath,
            "parameters": params,
        }

        # If overwrite=True, remove previous records of this datasource
        if overwrite:
            try:
                # delete by datasourceName (most common)
                self._project.deleteMeasurementsDocuments(
                    type="ToolkitDataSource", datasourceName=datasource_name
                )
            except Exception:
                # some implementations accept a free-text filter dict; try a broader cleanup
                try:
                    docs = self._project.getMeasurementsDocuments(type="ToolkitDataSource")
                    for d in docs:
                        dsn = (d.desc or {}).get("datasourceName")
                        tk = (d.desc or {}).get("toolkit")
                        if dsn == datasource_name or tk == datasource_name:
                            d.delete()
                except Exception:
                    # Best effort; if deletion API differs – we just continue
                    pass

        # Save the document
        doc = self._project.addMeasurementsDocument(
            resource=resource,
            dataFormat="Class",
            type="ToolkitDataSource",
            desc=desc,
        )
        return doc

    # ---------------------------
    # Lookups used by ToolkitHome
    # ---------------------------
    def getToolkitDocument(self, toolkit_name: str):
        """
        Find a dynamic toolkit document by name (either desc.datasourceName or desc.toolkit).
        Returns the mongoengine document or None.
        """
        # First: direct filter on datasourceName (works on most implementations)
        try:
            q = self._project.getMeasurementsDocuments(
                type="ToolkitDataSource", datasourceName=toolkit_name
            )
            if q and len(q) > 0:
                return q[0]
        except Exception:
            # fall through to broader search below
            pass

        # Second: scan all ToolkitDataSource docs and match by desc fields
        try:
            q = self._project.getMeasurementsDocuments(type="ToolkitDataSource")
            for d in q:
                desc = d.desc or {}
                if desc.get("datasourceName") == toolkit_name or desc.get("toolkit") == toolkit_name:
                    return d
        except Exception:
            pass

        # Optional: also look in DataSource collection if your project uses it
        try:
            q = self._project.getDataSourceDocuments(datasourceName=toolkit_name)
            if q and len(q) > 0:
                return q[0]
        except Exception:
            pass

        return None

    def getToolkitTable(self) -> pd.DataFrame:
        """
        Return a DataFrame with dynamic toolkits for the project.
        Columns: toolkit, cls, source, type, description
        """
        rows = []
        try:
            q = self._project.getMeasurementsDocuments(type="ToolkitDataSource")
        except Exception:
            q = []

        for d in q:
            desc = d.desc or {}
            toolkit = desc.get("datasourceName") or desc.get("toolkit") or ""
            cls = desc.get("classpath", "")
            rows.append(
                {
                    "toolkit": toolkit,
                    "cls": cls,
                    "source": desc.get("source", "external"),
                    "type": desc.get("toolkitType", desc.get("type", "measurements")),
                    "description": desc.get("description", ""),
                }
            )

        return pd.DataFrame(rows)
