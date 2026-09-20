"""absractEulerianSolver_toolkitExtension: flowType and the blockMesh_* setters.

B84/B85/B86, all resolved independently: an earlier pass found flowType's
else-branch referencing the undefined SIMULATIONTYPE_COMPRESSIBLE, and
blockMesh_setBoundFromFile/blockMesh_setDomainHeight both checking an
undefined eulerianWF instead of their real eulerianWorkFlow parameter.
Those were real, accurate findings at the time; the module now correctly
uses FLOWTYPE_COMPRESSIBLE and eulerianWorkFlow throughout. Verified
against a venv pinned to requirements.txt with hera.tests.unit._stubs
installed (which, after a related fix, gives hermes.workflow -- and so
workflow_Eulerian, which subclasses it transitively via OFWorkflow.py's
abstractWorkflow -- a real placeholder class instead of a bare MagicMock,
letting isinstance() checks against it work as intended).
"""
import pytest

from hera.simulations.openFoam.eulerian.abstractEulerianSolver import (
    absractEulerianSolver_toolkitExtension,
)
from hera.simulations.openFoam.OFWorkflow import workflow_Eulerian


class _FakeToolkit:
    FLOWTYPE_INCOMPRESSIBLE = "incompressible"


class _FakeWorkflow(workflow_Eulerian):
    def __init__(self):
        self.calls = []

    def set_blockMesh_blockHeight(self, Z, dz):
        self.calls.append(("height", Z, dz))

    def set_blockMesh_blockBoundaries(self, **kwargs):
        self.calls.append(("bounds", kwargs))


@pytest.fixture()
def ext():
    return absractEulerianSolver_toolkitExtension(
        toolkit=_FakeToolkit(), solverName="mySolver", incompressible=True
    )


@pytest.mark.unit
class TestConstruction:
    def test_it_stores_the_constructor_arguments(self, ext):
        assert ext.solverName == "mySolver"
        assert ext.incompressible is True


@pytest.mark.unit
class TestFlowType:
    def test_incompressible_true_reads_it_from_the_toolkit(self, ext):
        assert ext.flowType == "incompressible"

    def test_incompressible_false_returns_compressible(self):
        ext = absractEulerianSolver_toolkitExtension(
            toolkit=_FakeToolkit(), solverName="s", incompressible=False
        )
        assert ext.flowType == "compressible"


@pytest.mark.unit
class TestBlockMeshSetDomainHeight:
    def test_a_non_eulerian_workflow_is_rejected(self, ext):
        with pytest.raises(ValueError, match="not eulerian"):
            ext.blockMesh_setDomainHeight(eulerianWorkFlow=object(), Z=10, dz=1)

    def test_a_real_eulerian_workflow_gets_its_height_set(self, ext):
        wf = _FakeWorkflow()
        ext.blockMesh_setDomainHeight(eulerianWorkFlow=wf, Z=10, dz=1)
        assert wf.calls == [("height", 10, 1)]


@pytest.mark.unit
class TestBlockMeshSetBoundFromFile:
    def test_a_non_eulerian_workflow_is_rejected_before_reading_the_file(self, ext):
        """The file (a nonexistent path) is never touched -- validation
        happens first."""
        with pytest.raises(ValueError, match="not eulerian"):
            ext.blockMesh_setBoundFromFile(eulerianWorkFlow=object(), fileName="/no/such/file.obj", dx=1, dy=1, dz=1)
