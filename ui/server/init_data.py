"""Everything the UI needs at startup, gathered in one call (see ui/client/INIT_ENTRYPOINT_PLAN.md).

Each piece is collected on its own, so one broken piece (a toolkit document
without desc.classpath, a project whose documents fail to load) leaves the rest
usable instead of blanking the UI.
"""
from __future__ import annotations

import logging
import traceback
from typing import Any, Callable, Dict, List, Optional

# The same code the named project routes use, so the two cannot drift apart.
from project_api.project_api import DEFAULT_PROJECT, ProjectApi
from project_api.toolkit_documents import toolkit_documents

# Reuse uvicorn's logger so these lines match the rest of the server output.
logger = logging.getLogger("uvicorn.error")


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


def _datatypes() -> Dict[str, str]:
    from hera import datalayer
    return {key: value for key, value in vars(datalayer.datatypes).items()
            if not callable(value) and not key.startswith("__") and key[0] == key[0].upper()}


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

    projects = _collect(pieces, problems, "projects", ProjectApi.names)
    _collect(pieces, problems, "datatypes", _datatypes)
    _collect(pieces, problems, "toolkitDocs", toolkit_documents)
    _collect(pieces, problems, "server", _server_info)

    chosen = choose_project(projects, wanted_project)
    if chosen is not None:
        _collect(pieces, problems, "project", lambda: ProjectApi(chosen).documents())

    if problems:
        pieces["problems"] = problems
    return pieces
