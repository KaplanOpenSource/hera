"""Everything the UI needs at startup, gathered in one call (see ui/client/INIT_ENTRYPOINT_PLAN.md).

Each piece is collected on its own, so one broken piece (a toolkit document
without desc.classpath, a project whose documents fail to load) leaves the rest
usable instead of blanking the UI.
"""
from __future__ import annotations

import logging
import traceback
from typing import Any, Callable, Dict, List, Optional

# Reuse uvicorn's logger so these lines match the rest of the server output.
logger = logging.getLogger("uvicorn.error")

# Shown first in the project list; the client sorts the rest by name.
DEFAULT_PROJECT = "defaultProject"


def _collect(pieces: Dict[str, Any], problems: Dict[str, str], name: str, read: Callable[[], Any]) -> Any:
    """Run one reader, keeping its failure out of the other pieces."""
    try:
        value = read()
    except Exception as exc:
        logger.error("init: %s failed: %s", name, traceback.format_exc())
        problems[name] = f"{type(exc).__name__}: {exc}"
        return None
    pieces[name] = value
    return value


def _project_names() -> List[Dict[str, str]]:
    from hera.datalayer.project import getProjectList
    names = getProjectList()
    first = [{"name": n} for n in names if n == DEFAULT_PROJECT]
    rest = [{"name": n} for n in names if n != DEFAULT_PROJECT]
    rest.sort(key=lambda p: p["name"].lower())
    return first + rest


def _datatypes() -> Dict[str, str]:
    from hera import datalayer
    return {key: value for key, value in vars(datalayer.datatypes).items()
            if not callable(value) and not key.startswith("__") and key[0] == key[0].upper()}


def _toolkit_documents() -> List[Dict[str, Any]]:
    from hera import toolkitHome
    return toolkitHome.getToolkitDocuments()


def _project_documents(projectName: str) -> Dict[str, Any]:
    from hera.datalayer import All
    docs = All.getDocumentsAsDict(projectName, with_id=True)
    return {"name": projectName, "documents": docs["documents"]}


def _server_info() -> Dict[str, Any]:
    import getpass
    import os
    import hera
    return {
        "heraVersion": hera.__version__,
        "username": getpass.getuser(),
        "inDocker": os.path.exists("/.dockerenv"),
    }


def choose_project(projects: Optional[List[Dict[str, str]]], wanted: str) -> Optional[str]:
    """The project to open: the wanted one if it exists, else the first, else none."""
    names = [p["name"] for p in projects or []]
    if wanted and wanted in names:
        return wanted
    return names[0] if names else None


def collect_init(wanted_project: str = "") -> Dict[str, Any]:
    """The startup payload: project list, datatypes, toolkits, one project, server info.

    `problems` names the pieces that failed, if any; those pieces are absent.
    """
    pieces: Dict[str, Any] = {"ready": True}
    problems: Dict[str, str] = {}

    projects = _collect(pieces, problems, "projects", _project_names)
    _collect(pieces, problems, "datatypes", _datatypes)
    _collect(pieces, problems, "toolkitDocs", _toolkit_documents)
    _collect(pieces, problems, "server", _server_info)

    chosen = choose_project(projects, wanted_project)
    if chosen is not None:
        _collect(pieces, problems, "project", lambda: _project_documents(chosen))

    if problems:
        pieces["problems"] = problems
    return pieces
