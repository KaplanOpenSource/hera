"""absractEulerianSolver_toolkitExtension: three defects in a row.

* B84: ``flowType``'s else-branch references ``SIMULATIONTYPE_COMPRESSIBLE``,
  a name that does not exist anywhere in the module (the import at the top
  is ``FLOWTYPE_COMPRESSIBLE``). Every call with ``incompressible=False``
  raises ``NameError``.
* B85: ``blockMesh_setBoundFromFile``'s parameter is ``eulerianWorkFlow``,
  but its body checks ``isinstance(eulerianWF, workflow_Eulerian)`` --
  ``eulerianWF`` is not a parameter of this method at all. Every call
  raises ``NameError`` before touching any argument.
* B86: ``blockMesh_setDomainHeight`` is a copy-paste of
  ``blockMesh_setBoundFromFile`` that was never adapted: it also checks
  the undefined ``eulerianWF``, and its body references ``fileName``,
  ``dx`` and ``dy`` -- none of which are parameters of this method
  (``eulerianWorkFlow, Z, dz``). ``Z``, its only real input, is never used.
"""
import pytest

from hera.simulations.openFoam.eulerian.abstractEulerianSolver import (
    absractEulerianSolver_toolkitExtension,
)


class _FakeToolkit:
    FLOWTYPE_INCOMPRESSIBLE = "incompressible"


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
class TestFlowTypeIsBroken:
    def test_incompressible_true_reads_it_from_the_toolkit(self, ext):
        assert ext.flowType == "incompressible"

    @pytest.mark.xfail(
        strict=True,
        reason="B84: the else-branch references SIMULATIONTYPE_COMPRESSIBLE, "
               "a name that exists nowhere in this module. "
               "See the consolidated findings issue.",
    )
    def test_incompressible_false_should_also_work(self):
        ext = absractEulerianSolver_toolkitExtension(
            toolkit=_FakeToolkit(), solverName="s", incompressible=False
        )
        ext.flowType

    def test_incompressible_false_currently_raises(self):
        """Characterisation of B84."""
        ext = absractEulerianSolver_toolkitExtension(
            toolkit=_FakeToolkit(), solverName="s", incompressible=False
        )
        with pytest.raises(NameError, match="SIMULATIONTYPE_COMPRESSIBLE"):
            ext.flowType


@pytest.mark.unit
class TestBlockMeshSetBoundFromFileIsBroken:
    @pytest.mark.xfail(
        strict=True,
        reason="B85: the body checks isinstance(eulerianWF, ...), but the "
               "method's parameter is named eulerianWorkFlow -- eulerianWF "
               "is not defined anywhere. See the consolidated findings issue.",
    )
    def test_it_should_validate_its_actual_parameter(self, ext):
        ext.blockMesh_setBoundFromFile(eulerianWorkFlow=object(), fileName="x", dx=1, dy=1, dz=1)

    def test_it_currently_raises_nameerror_before_touching_any_argument(self, ext):
        """Characterisation of B85."""
        with pytest.raises(NameError, match="eulerianWF"):
            ext.blockMesh_setBoundFromFile(eulerianWorkFlow=object(), fileName="x", dx=1, dy=1, dz=1)


@pytest.mark.unit
class TestBlockMeshSetDomainHeightIsBroken:
    @pytest.mark.xfail(
        strict=True,
        reason="B86: copy-pasted from blockMesh_setBoundFromFile without "
               "adaptation -- checks the same undefined eulerianWF, and "
               "references fileName/dx/dy, none of which are parameters "
               "of this method (eulerianWorkFlow, Z, dz). Z is never used. "
               "See the consolidated findings issue.",
    )
    def test_it_should_set_the_domain_height(self, ext):
        ext.blockMesh_setDomainHeight(eulerianWorkFlow=object(), Z=10, dz=1)

    def test_it_currently_raises_nameerror_before_touching_z(self, ext):
        """Characterisation of B86."""
        with pytest.raises(NameError, match="eulerianWF"):
            ext.blockMesh_setDomainHeight(eulerianWorkFlow=object(), Z=10, dz=1)
