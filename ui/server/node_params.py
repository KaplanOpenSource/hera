"""Validation of a Hermes node's parameter values.

Calls the node executer's own ``testParamValues`` (see
``Hermes/hermes/Resources/executers/abstractExecuter.py``), so the message the
user sees is Hermes's, not a copy of its rules.

Unlike the node catalog, this really imports the executer class and the node's
checks reach hera and MongoDB, so it only works in the warmed-up server process.
"""
from __future__ import annotations

import importlib
import re
import sys
from pathlib import Path

# A node type is the Resources directory path with dots, e.g.
# "RiskAssessment.calculateThresholds". The type arrives from the client, so only
# this shape is ever imported.
_NODE_TYPE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*(\.[A-Za-z_][A-Za-z0-9_]*)*$")

# node_params.py lives in <repo>/ui/server, so the Hermes submodule is two up.
_HERMES_ROOT = Path(__file__).resolve().parents[2] / "Hermes"


def _executer_class(node_type: str):
    """The executer class for a node type, or None when the type has no executer."""
    if not _NODE_TYPE.match(node_type):
        return None
    node_dir = _HERMES_ROOT / "hermes" / "Resources" / Path(*node_type.split("."))
    if not (node_dir / "executer.py").is_file():
        return None
    if str(_HERMES_ROOT) not in sys.path:
        sys.path.insert(0, str(_HERMES_ROOT))
    module = importlib.import_module(f"hermes.Resources.{node_type}.executer")
    # The class is named after the node's own directory.
    return getattr(module, node_type.rsplit(".", 1)[-1], None)


def validate_node_params(node_type: str, params: dict) -> dict:
    """Test a node's parameter values against its own executer.

    Returns ``{"ok": bool, "message": str}``. ok with an empty message also covers
    "nothing to say": an unknown type, or a node that defines no checks of its
    own. Only a class that overrides ``testParamValues`` is called, because the
    base stub is declared static yet takes ``self``.
    """
    executer = _executer_class(node_type)
    if executer is None or "testParamValues" not in executer.__dict__:
        return {"ok": True, "message": ""}
    ok, message = executer.testParamValues(params)
    return {"ok": bool(ok), "message": message or ""}
