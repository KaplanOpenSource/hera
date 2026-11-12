"""``sys.modules`` stubs for packages that are not installed.

Nineteen hera modules import PyFoam, paraview, FreeCAD, hermes, argos or
evtk at module level.  Registering placeholder modules before those imports
run lets the surrounding pure-Python logic be tested without the external
binaries.

Every stub here is conditional: a placeholder is registered only when the
real module cannot be imported, so an environment that has the package --
CI installs PyFoam and pyevtk from requirements.txt -- tests against the
actual classes and the stub never gets in the way.

Two kinds of placeholder, because the two are not interchangeable:

* a MagicMock is enough for a module whose attributes are only called or
  read;
* a module carrying real (if empty) classes is needed where an attribute is
  used as a BASE class, since ``class Foo(mod.Bar):`` requires ``mod.Bar``
  to be a type.  Subclassing a MagicMock attribute does not raise -- it
  silently produces another MagicMock as the "class", which then fails
  ``isinstance()`` elsewhere with ``TypeError: isinstance() arg 2 must be a
  type``.  ``argos.experimentSetup.dataObjects`` (experiment.py) and
  ``hermes.workflow`` (OFWorkflow.py) both need this.

Deliberately NOT stubbed: ``torch``.  Verified in batch 9: a MagicMock is
not enough, because modelContainer.py reaches into submodules and the
import fails with "No module named 'torch.utils'; 'torch' is not a
package".  torch is declared in requirements.txt, so CI has the real thing.
"""
import sys
import types
from unittest.mock import MagicMock

# name -> the placeholder to register when `import name` fails.
#
#   "namespace"          a package, i.e. it must expose __path__ so that
#                        `import name.submodule` works
#   "mock"               a leaf module; every attribute is a MagicMock
#   (class names tuple)  a leaf module exposing those names as real classes
_STUBS = (
    ("PyFoam", "namespace"),
    ("PyFoam.RunDictionary", "namespace"),
    ("PyFoam.Basics", "namespace"),
    ("PyFoam.RunDictionary.ParsedParameterFile", "mock"),
    ("PyFoam.RunDictionary.BoundaryDict", "mock"),
    ("PyFoam.Basics.DataStructures", "mock"),
    ("paraview", "namespace"),
    ("paraview.simple", "namespace"),
    ("evtk", "namespace"),
    ("evtk.hl", "mock"),
    ("dask.distributed", "mock"),
    ("argos", "namespace"),
    ("argos.experimentSetup", "namespace"),
    ("argos.experimentSetup.dataObjects",
     ("ExperimentZipFile", "TrialSet", "Trial", "EntityType", "Entity")),
    ("hermes", ("workflow",)),
    ("FreeCAD", "mock"),
)


def _isImportable(name):
    if name in sys.modules:
        return True
    try:
        __import__(name)
    except Exception:
        return False
    return True


def _register(name, kind):
    if kind == "mock":
        sys.modules[name] = MagicMock()
        return
    module = types.ModuleType(name)
    if kind == "namespace":
        module.__path__ = []
        module.__package__ = name
    else:
        for className in kind:
            setattr(module, className, type(className, (), {}))
    sys.modules[name] = module


def install():
    """Register every stub. Idempotent — safe to call more than once."""
    for name, kind in _STUBS:
        if not _isImportable(name):
            _register(name, kind)

    # dask.dataframe must be imported before analysis.addDatesColumns runs.
    # Order matters here; see test_toolkit_coverage.py's original note.
    import dask.dataframe  # noqa: F401
