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
