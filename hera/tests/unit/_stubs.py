"""``sys.modules`` stubs for packages that cannot be installed in CI.

Nineteen hera modules import PyFoam, paraview, FreeCAD, hermes, argos or
evtk at module level.  Registering placeholder modules before those imports
run lets the surrounding pure-Python logic be tested without the external
binaries.

Deliberately NOT stubbed: ``torch``.  modelContainer.py subclasses
``torch.nn.Module``, and subclassing a MagicMock raises TypeError.  Batch 9
addresses that separately.
"""
import sys
import types
from unittest.mock import MagicMock

# Namespace packages: must expose __path__ so that `import a.b` works.
_NAMESPACE_PACKAGES = (
    "PyFoam",
    "PyFoam.RunDictionary",
    "PyFoam.Basics",
    "paraview",
    "paraview.simple",
)

# Leaf modules whose attributes are only ever called, never subclassed.
_MOCK_MODULES = (
    "PyFoam.RunDictionary.ParsedParameterFile",
    "PyFoam.RunDictionary.BoundaryDict",
    "PyFoam.Basics.DataStructures",
    "evtk",
    "evtk.hl",
    "dask.distributed",
)

# Stubbed only when genuinely absent, so a real local install wins.
_STUB_IF_MISSING = ("hermes", "argos", "FreeCAD")


def _namespace_module(name):
    module = types.ModuleType(name)
    module.__path__ = []
    module.__package__ = name
    sys.modules[name] = module
    return module


def install():
    """Register every stub. Idempotent — safe to call more than once."""
    for name in _NAMESPACE_PACKAGES:
        if name not in sys.modules:
            _namespace_module(name)

    for name in _MOCK_MODULES:
        if name not in sys.modules:
            sys.modules[name] = MagicMock()

    for name in _STUB_IF_MISSING:
        if name in sys.modules:
            continue
        try:
            __import__(name)
        except Exception:
            sys.modules[name] = MagicMock()

    # dask.dataframe must be imported before analysis.addDatesColumns runs.
    # Order matters here; see test_toolkit_coverage.py's original note.
    import dask.dataframe  # noqa: F401
