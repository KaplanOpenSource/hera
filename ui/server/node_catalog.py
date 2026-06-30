"""Hermes node catalog — the node types available to the workflow editor.

Loads the vendored Hermes node lookup (``Hermes/hermes/utils/node_lookup.py``)
straight from its file — it only needs jinja2 + the stdlib — so it works even
when ``hermes`` is not installed and depends on nothing in the hera package.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Optional

_CATALOG_CACHE: Optional[list] = None


def get_node_catalog() -> list:
    """Available Hermes node types and the parameters each one accepts.

    Each entry is ``{"type": <dotted name>, "parameters": [{"name", "is_required",
    "source"}, ...]}``, sorted by type. Parameter values (defaults) are not
    included — node_lookup carries only names. Cached in-process; the Resources
    tree is static at runtime.
    """
    global _CATALOG_CACHE
    if _CATALOG_CACHE is None:
        # node_catalog.py lives in <repo>/ui/server, so the Hermes submodule is
        # two levels up.
        hermes_package = Path(__file__).resolve().parents[2] / "Hermes" / "hermes"
        node_lookup_path = hermes_package / "utils" / "node_lookup.py"
        if not node_lookup_path.is_file():
            raise FileNotFoundError(f"Hermes node lookup not found at {node_lookup_path}; is the Hermes submodule checked out?")

        spec = importlib.util.spec_from_file_location("hermes_node_lookup", node_lookup_path)
        if spec is None or spec.loader is None:
            raise ImportError(f"Could not load Hermes node lookup from {node_lookup_path}")
        node_lookup = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(node_lookup)

        # node_lookup's own default Resources path is computed relative to its
        # module and lands one directory too high, so pass the real location.
        resources_root = str(hermes_package / "Resources")
        catalog = [info.to_dict() for info in node_lookup.get_all_node_types(resources_root).values()]
        catalog.sort(key=lambda entry: entry["type"])
        _CATALOG_CACHE = catalog
    return _CATALOG_CACHE
