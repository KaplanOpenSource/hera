"""``sys.modules`` stubs for packages that cannot be installed in CI.

Nineteen hera modules import PyFoam, paraview, FreeCAD, hermes, argos or
evtk at module level.  Registering placeholder modules before those imports
run lets the surrounding pure-Python logic be tested without the external
binaries.

``argos.experimentSetup.dataObjects`` needs the same namespace-package
treatment as PyFoam, plus real (if empty) classes rather than a leaf
MagicMock: experiment.py subclasses ``ExperimentZipFile``, ``TrialSet``,
``Trial``, ``EntityType`` and ``Entity`` at module level, and a class
statement's bases must be actual types. Unlike ``argos.experimentSetup.
dataObjects.ExperimentZipFile``, subclassing a bare ``MagicMock()``
attribute (``hermes.workflow`` before this fix) does not raise -- it
silently produces another MagicMock as the "class" instead, which then
fails ``isinstance()`` checks against it elsewhere with a confusing
``TypeError: isinstance() arg 2 must be a type, a tuple of types, or a
union``. ``hermes.workflow`` gets the same real-placeholder-class
treatment for exactly this reason: ``OFWorkflow.py``'s
``abstractWorkflow(hermes.workflow)`` needs a real base class.

Deliberately NOT stubbed: ``torch``.  Verified in batch 9: a leaf MagicMock
is not enough, because modelContainer.py reaches submodules and the import
fails with "No module named 'torch.utils'; 'torch' is not a package".  It
would need the namespace-package treatment PyFoam gets, for torch, torch.nn,
torch.utils and whatever else the module walks into.  That is tractable but
was out of scope; scikit-learn-style optional deps are declared in
requirements.txt, and torch is too, so CI has the real thing.

(An earlier version of this note claimed modelContainer.py subclasses
torch.nn.Module and that inheriting from a MagicMock raises TypeError.  It
does not subclass it -- the reference is in a docstring -- so the reason
recorded above replaces that one.)
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

# Namespace packages stubbed only when genuinely absent, so a real local
# install of the same package wins (unlike _NAMESPACE_PACKAGES above, which
# is unconditional because nothing in this repo ever has a real PyFoam or
# paraview install to prefer).
_NAMESPACE_PACKAGES_IF_MISSING = (
    "argos",
    "argos.experimentSetup",
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

# Leaf modules that must expose real (if empty) classes, because they are
# used as base classes -- `class Foo(mod.Bar):` requires `mod.Bar` to be a
# type, unlike a plain attribute access or call, which a MagicMock covers.
_CLASS_STUB_MODULES = {
    "argos.experimentSetup.dataObjects": (
        "ExperimentZipFile",
        "TrialSet",
        "Trial",
        "EntityType",
        "Entity",
    ),
}

# Top-level modules stubbed only when genuinely absent (real local install
# wins), where a flat MagicMock module is enough for everything except one
# or more attributes that must be real classes -- see _CLASS_STUB_MODULES's
# docstring above for why.
_TOP_LEVEL_CLASS_ATTRS_IF_MISSING = {
    "hermes": ("workflow",),
}

# Stubbed only when genuinely absent, so a real local install wins.
_STUB_IF_MISSING = ("FreeCAD",)


def _namespace_module(name):
    module = types.ModuleType(name)
    module.__path__ = []
    module.__package__ = name
    sys.modules[name] = module
    return module


def _class_stub_module(name, class_names):
    module = types.ModuleType(name)
    for class_name in class_names:
        setattr(module, class_name, type(class_name, (), {}))
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

    if "argos" not in sys.modules:
        try:
            __import__("argos.experimentSetup.dataObjects")
        except Exception:
            for name in _NAMESPACE_PACKAGES_IF_MISSING:
                _namespace_module(name)
            for name, class_names in _CLASS_STUB_MODULES.items():
                _class_stub_module(name, class_names)

    for name, class_names in _TOP_LEVEL_CLASS_ATTRS_IF_MISSING.items():
        if name in sys.modules:
            continue
        try:
            __import__(name)
        except Exception:
            _class_stub_module(name, class_names)

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
