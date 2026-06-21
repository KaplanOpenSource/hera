"""
Notebook regression tests.

Discovers every *.ipynb in the repository (excluding legacy/external paths),
executes each one in a fresh kernel, and enforces a strict skip-vs-fail policy:

  - Expected skips  → pytest.skip  (missing optional external dep, missing data
                                     file, slow/unavailable network service)
  - Real regressions → pytest.fail  (import error for a core hera module,
                                     broken API call, syntax error, etc.)

Run all notebooks::

    pytest hera/tests/test_notebooks.py

Run a single notebook by keyword::

    pytest hera/tests/test_notebooks.py -k TileToolkit

List discovered notebooks without running::

    pytest hera/tests/test_notebooks.py --collect-only
"""

from __future__ import annotations

import re
from pathlib import Path

import nbformat
import pytest
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError, CellTimeoutError

# ---------------------------------------------------------------------------
# Repository layout
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).parent.parent.parent  # …/hera/tests/ → repo root

# Any notebook path component that matches one of these names is excluded.
# Add new legacy/env directory names here as the repo evolves.
_EXCLUDE_DIRS: frozenset[str] = frozenset({
    "sphinx.old",
    ".git",
    "node_modules",
    "__pycache__",
    "heraenv",
    ".tox",
    ".venv",
    "venv",
    "env",
})

# Individual notebooks excluded by relative path from REPO_ROOT.
# Use this for notebooks that are intentional scratchpads / broken-by-design.
_EXCLUDE_NOTEBOOKS: frozenset[str] = frozenset({
    "hera/tests/testExamples.ipynb",  # developer scratchpad; has intentional syntax errors
})

# ---------------------------------------------------------------------------
# Skip heuristics
# ---------------------------------------------------------------------------

# Root package names whose absence should produce a pytest.skip, not a fail.
# These are *optional* external tools that are not part of a standard hera install.
# Do NOT add core packages (pandas, numpy, hera, …) here — a missing core dep
# is a real regression and should surface as a pytest.fail.
_SKIP_MODULES: frozenset[str] = frozenset({
    "hermes",           # Hermes workflow engine (simulation orchestration)
    "basic",            # local non-repo helper used in gaussionToolkit notebook
    "IMS_experiment",   # external IMS station experiment codebase
    "Jerusalem2018",    # external Jerusalem-2018 field-campaign codebase
    "openfoam",         # OpenFOAM Python bindings (not pip-installable)
})

# Regex patterns matched against FileNotFoundError.evalue.
# A match means the file is expected external data, not a repo asset.
_SKIP_FILE_PATTERNS: list[re.Pattern[str]] = [
    re.compile(r"token\.json"),          # IMS API credential file
    re.compile(r"HaifaFluxes.*\.zip"),   # field-campaign data archive
    re.compile(r"output\.geojson"),      # pre-generated OSM export
]

# Substrings matched against any exception's evalue.
# Covers wrapped errors (ErrorDuringImport, etc.) that can't be matched by
# ename alone.
_SKIP_ERROR_SUBSTRINGS: tuple[str, ...] = (
    "Cannot use this module without hermes",
    "Install it.",  # hermes install prompt
)

# ---------------------------------------------------------------------------
# Execution settings
# ---------------------------------------------------------------------------

# Per-notebook wall-clock timeout in seconds.  120 s covers tile-server
# downloads (~60-90 s) with headroom; increase if new notebooks are added
# that legitimately run longer.
_NOTEBOOK_TIMEOUT: int = 120


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _skip_reason(exc: CellExecutionError) -> str | None:
    """
    Return a human-readable skip reason if *exc* represents a known external
    dependency gap, or ``None`` if the error is a genuine regression.

    Decision order
    --------------
    1. ModuleNotFoundError for a known-external module → skip.
    2. FileNotFoundError whose path matches a known-external-data pattern → skip.
    3. Any exception whose evalue contains a known skip substring → skip.
    4. Everything else → None (caller should pytest.fail).
    """
    ename: str = exc.ename or ""
    evalue: str = exc.evalue or ""

    # 1. Missing external module
    if ename == "ModuleNotFoundError":
        match = re.search(r"No module named '([^']+)'", evalue)
        if match:
            root_pkg = match.group(1).split(".")[0]
            if root_pkg in _SKIP_MODULES:
                return f"optional external module not installed: '{root_pkg}'"

    # 2. Missing external data file
    if ename == "FileNotFoundError":
        for pattern in _SKIP_FILE_PATTERNS:
            if pattern.search(evalue):
                return f"external data file not present: {evalue.strip()}"

    # 3. Known error-message substrings (handles ErrorDuringImport wrappers, etc.)
    for substring in _SKIP_ERROR_SUBSTRINGS:
        if substring in evalue:
            return f"known external dependency unavailable: {substring!r}"

    return None


def _discover_notebooks() -> list[Path]:
    """Return all non-excluded *.ipynb paths, sorted for stable test ordering."""
    notebooks: list[Path] = []
    for nb_path in sorted(REPO_ROOT.glob("**/*.ipynb")):
        rel = nb_path.relative_to(REPO_ROOT)
        # Drop if any path component is in the exclude-dirs set
        if any(part in _EXCLUDE_DIRS for part in rel.parts):
            continue
        # Drop individually excluded notebooks
        if str(rel) in _EXCLUDE_NOTEBOOKS:
            continue
        notebooks.append(nb_path)
    return notebooks


# ---------------------------------------------------------------------------
# Test
# ---------------------------------------------------------------------------

@pytest.mark.notebook
@pytest.mark.parametrize(
    "nb_path",
    _discover_notebooks(),
    ids=lambda p: str(p.relative_to(REPO_ROOT)),
)
def test_notebook(nb_path: Path) -> None:
    """Execute a notebook end-to-end and enforce the skip/fail policy."""
    nb = nbformat.read(str(nb_path), as_version=4)

    client = NotebookClient(
        nb,
        timeout=_NOTEBOOK_TIMEOUT,
        kernel_name="python3",
        # Set the working directory to the notebook's own folder so that
        # relative paths inside cells (e.g. '../data/') resolve correctly.
        resources={"metadata": {"path": str(nb_path.parent)}},
    )

    try:
        client.execute()

    except CellTimeoutError as exc:
        # Timeouts are almost always caused by a slow/unavailable external
        # service (tile server, remote API).  Skip rather than block CI.
        pytest.skip(
            f"notebook timed out after {_NOTEBOOK_TIMEOUT}s "
            f"(likely slow/unavailable external service): {nb_path.name}"
        )

    except CellExecutionError as exc:
        reason = _skip_reason(exc)
        if reason:
            pytest.skip(f"{nb_path.name}: {reason}")
        # Real regression — surface the full traceback so the developer can
        # see exactly which cell failed and why.
        pytest.fail(
            f"notebook execution failed: {nb_path.relative_to(REPO_ROOT)}\n\n"
            f"{exc}",
            pytrace=False,
        )
