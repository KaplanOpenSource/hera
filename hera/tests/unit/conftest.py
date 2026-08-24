"""Hermetic unit layer: no MongoDB, no S3 data, no network.

The module-level block below runs at collection time, before any test module
is imported.  Its order is load-bearing — see hera/tests/unit/_seam.py.
"""
# --- bootstrap: the order of these five statements is load-bearing --------
#
# 1. HOME moves to a temp directory.  This must precede every import of
#    hera.datalayer, so only the standard library may be touched above this
#    point — hence these lines are inline rather than a call into _seam.
#
#    Known limitation: this does NOT precede `import hera` itself.  Because
#    this directory is a package, pytest imports the chain
#    hera -> hera.tests -> hera.tests.unit before executing this module, and
#    hera/__init__.py calls initialize_logging() eagerly.  Log files
#    therefore still land in the developer's real ~/.pyhera/log — the same
#    place any `import hera` writes them, so nothing test-specific leaks.
#    Everything that matters for isolation (the pyhera connection config and
#    every hera.datalayer import) happens after this point.
import os
import pathlib
import tempfile

_UNIT_HOME = tempfile.mkdtemp(prefix="hera_unit_home_")
os.environ["HOME"] = _UNIT_HOME

# 2. Importing hera.tests.unit is safe now: hera/__init__.py is fully lazy
#    (it exposes datalayer through __getattr__ only), hera/tests/__init__.py
#    is empty, and _seam/_stubs import nothing but the standard library.
from hera.tests.unit import _seam, _stubs

# 3. The pyhera config must exist before hera.datalayer is imported, because
#    hera/datalayer/__init__.py:7-10 builds collection singletons eagerly.
_seam.write_pyhera_config()

# 4. Stub the CI-unavailable packages before any module that imports them.
_stubs.install()

# 5. Now route the datalayer at mongomock.  This is the first statement that
#    causes hera.datalayer to be imported.
_seam.install()
# --------------------------------------------------------------------------

import pytest

UNIT_PROJECT_NAME = "UNIT_TEST_PROJECT"


@pytest.fixture(scope="session", autouse=True)
def _no_trace_guard():
    """Neutralise the parent conftest's MongoDB-purging session guard.

    hera/tests/conftest.py defines a session-scoped autouse fixture of this
    name that enumerates and purges projects in the real MongoDB.  A parent
    conftest's autouse fixtures apply to subdirectories too, so without this
    override every unit run pays pymongo's 30-second server-selection
    timeout twice — measured at 62.95s before this override existed.

    Overriding by name is the whole mechanism: pytest resolves fixtures from
    the closest conftest outward, so this definition wins for everything
    under hera/tests/unit/ and the parent's version never runs here.  The
    unit layer has nothing to purge: its database lives in memory and is
    dropped after every test by _reset_unit_database.
    """
    yield


@pytest.fixture()
def unit_files_directory(tmp_path):
    """A per-test files directory, so nothing is written to the real home."""
    directory = tmp_path / "heraFiles"
    directory.mkdir()
    return str(directory)


@pytest.fixture()
def unit_project(unit_files_directory):
    """A real hera Project backed by the in-memory database."""
    from hera.datalayer import Project

    return Project(projectName=UNIT_PROJECT_NAME, filesDirectory=unit_files_directory)


@pytest.fixture()
def unit_toolkit_factory(unit_files_directory):
    """Build any real toolkit, backed by the in-memory database.

    Usage:
        toolkit = unit_toolkit_factory(toolkitHome.GIS_RASTER_TOPOGRAPHY)
    """
    from hera import toolkitHome

    def _build(toolkitName, projectName=UNIT_PROJECT_NAME):
        return toolkitHome.getToolkit(
            toolkitName,
            projectName=projectName,
            filesDirectory=unit_files_directory,
        )

    return _build


# ---------------------------------------------------------------------------
# Guards: turn the hermetic promise into something enforced
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _block_network(monkeypatch):
    """Fail loudly and instantly on any socket creation.

    Without this, a unit test that reaches for MongoDB waits out pymongo's
    30-second serverSelectionTimeoutMS and reports as merely slow rather
    than broken.  That is what cost this very suite 62.95s per run before
    the parent conftest's session guard was overridden.
    """
    import socket

    def _refuse(*args, **kwargs):
        raise RuntimeError(
            "network access is not permitted in the unit layer; "
            "use the unit_project / unit_toolkit_factory fixtures instead"
        )

    # mongomock never opens a socket, so nothing legitimate needs the real one.
    monkeypatch.setattr(socket, "socket", _refuse)


@pytest.fixture(autouse=True)
def _matplotlib_agg():
    """Force a headless backend and close figures, to stop state leaking."""
    import matplotlib

    matplotlib.use("Agg", force=True)
    yield
    import matplotlib.pyplot as plt

    plt.close("all")


@pytest.fixture(autouse=True)
def _no_home_writes():
    """Assert that no test creates a .hera directory in the isolated home."""
    import pathlib

    hera_dir = pathlib.Path(_UNIT_HOME, ".hera")
    existed = hera_dir.exists()
    yield
    if not existed and hera_dir.exists():
        raise AssertionError(
            f"a test created {hera_dir}; pass filesDirectory (see the "
            "unit_files_directory fixture) instead of relying on the default"
        )


@pytest.fixture(autouse=True)
def _reset_unit_database():
    """Drop the in-memory database after every test.

    Test order must never matter.  A test that only passes because an
    earlier test left a document behind is a bug, not a feature.
    """
    yield
    _seam.reset()


def pytest_collection_modifyitems(session, config, items):
    """Refuse to run the unit layer in the same process as any other test.

    The bootstrap at the top of this file is process-wide: it moves HOME,
    replaces hera.datalayer.document.connectToDatabase and rebinds the
    datalayer's collection singletons.  Any integration test collected
    alongside it therefore talks to mongomock, or authenticates against the
    real MongoDB with the placeholder credentials from the temp pyhera
    config.

    That was not theoretical.  Running `pytest hera/tests -m "not notebook"`
    with this directory collected produced three errors in
    dynamic_loading_tests_pack (subprocesses inheriting the moved HOME ->
    "Authentication failed") that all disappeared once the directory was
    ignored.

    -m "not unit" does NOT prevent this: marker filtering happens after
    collection, and the damage is done at import time.  The directory has to
    be excluded outright, so run the two layers as two processes:

        pytest hera/tests/unit -m unit
        pytest hera/tests --ignore=hera/tests/unit -m "not notebook"

    or use `make test-unit` and `make test`, which already do exactly that.
    """
    unit_root = pathlib.Path(__file__).parent.resolve()
    strangers = sorted(
        {
            str(item.path)
            for item in items
            if unit_root not in pathlib.Path(item.path).resolve().parents
        }
    )
    if strangers:
        raise pytest.UsageError(
            "hera/tests/unit must run in its own pytest process; its conftest "
            "moves HOME and reroutes the datalayer for the whole process.\n"
            f"Also collected {len(strangers)} file(s) outside it, first few: "
            f"{strangers[:3]}\n"
            "Run:  pytest hera/tests/unit -m unit\n"
            "and:  pytest hera/tests --ignore=hera/tests/unit -m 'not notebook'\n"
            "(or use `make test-unit` and `make test`)."
        )
