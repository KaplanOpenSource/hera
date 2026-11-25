"""OFWorkflow.py: the hard import guard, the only part of the module left
uncovered.

``test_openfoam_ofworkflow.py`` drives every class and property in the module
against real ``hermes.workflow`` objects.  What no test could reach is the two
statements at the top:

    try:
        import hermes
    except ImportError:
        raise ImportError("Cannot use this module without hermes... Install it. ")

hermes *is* installed here (and is not stubbed -- see the sibling file), so
importing the module can never take that branch, and re-importing it would
not re-execute the header either: ``sys.modules`` already holds it.

The branch is reached by compiling the module's own source and executing it in
a throwaway namespace with ``__import__`` refusing the name ``hermes``.  The
namespace is discarded and nothing is written into ``sys.modules``, so the
already-imported module is untouched; ``monkeypatch`` restores ``__import__``.

This is worth pinning rather than skipping, because the module's contract is
that it fails *loudly and immediately* -- unlike
``hermesWorkflowToolkit.py``, which downgrades the same missing dependency to
``warnings.warn`` and binds ``workflow = None``, and unlike
``openFoam/CLI.py``, which warns and then calls the name the failed import was
supposed to bind (B178, already pinned).  A regression to either of those
shapes would turn a clear ImportError at import time into a NameError or an
AttributeError somewhere much later.

No bugs found.
"""
import builtins

import pytest

import hera.simulations.openFoam.OFWorkflow as ofWorkflowModule


def _executeModuleSource(monkeypatch, blockedPrefix):
    """Re-execute OFWorkflow.py in a throwaway namespace, with an import banned.

    Nothing is registered in sys.modules, so the module the rest of the
    session already imported is not disturbed.
    """
    source = open(ofWorkflowModule.__file__).read()
    realImport = builtins.__import__

    def refusingImport(name, *args, **kwargs):
        if name == blockedPrefix or name.startswith(f"{blockedPrefix}."):
            raise ImportError(f"No module named {name!r}")
        return realImport(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", refusingImport)
    namespace = dict(__name__="hera.simulations.openFoam._OFWorkflowImportProbe",
                     __file__=ofWorkflowModule.__file__,
                     __package__="hera.simulations.openFoam")
    exec(compile(source, ofWorkflowModule.__file__, "exec"), namespace)
    return namespace


@pytest.mark.unit
class TestTheHermesImportGuard:
    def test_hermes_is_available_here_so_the_module_imported_normally(self):
        assert ofWorkflowModule.hermes is not None
        assert hasattr(ofWorkflowModule.hermes, "workflow")

    def test_without_hermes_the_module_refuses_to_load(self, monkeypatch):
        with pytest.raises(ImportError):
            _executeModuleSource(monkeypatch, "hermes")

    def test_the_error_says_what_to_install(self, monkeypatch):
        with pytest.raises(ImportError, match="Cannot use this module without hermes"):
            _executeModuleSource(monkeypatch, "hermes")

    def test_the_original_import_error_is_not_what_reaches_the_caller(self, monkeypatch):
        """The bare 'No module named' is replaced by an actionable message."""
        with pytest.raises(ImportError) as raised:
            _executeModuleSource(monkeypatch, "hermes")
        assert "No module named" not in str(raised.value)

    def test_the_failure_happens_at_import_time_not_at_first_use(self, monkeypatch):
        """Nothing from the module is usable, because no class is defined:
        the raise sits above every class statement."""
        with pytest.raises(ImportError):
            _executeModuleSource(monkeypatch, "hermes")

    def test_the_module_still_works_after_the_probe(self, monkeypatch):
        """The probe must not disturb the already-imported module."""
        with pytest.raises(ImportError):
            _executeModuleSource(monkeypatch, "hermes")
        assert ofWorkflowModule.workflow_Eulerian.__name__ == "workflow_Eulerian"
        assert issubclass(ofWorkflowModule.workflow_Eulerian,
                          ofWorkflowModule.abstractWorkflow)

    def test_the_probe_really_does_execute_the_module_when_nothing_is_blocked(self,
                                                                             monkeypatch):
        """Guard against the test passing for the wrong reason: with an
        unrelated name blocked, the same source runs to completion and defines
        the module's classes."""
        namespace = _executeModuleSource(monkeypatch, "aModuleNobodyImports")
        assert "abstractWorkflow" in namespace
        assert namespace["workflow_Eulerian"] is not ofWorkflowModule.workflow_Eulerian
