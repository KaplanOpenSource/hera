"""absractEulerianSolver_toolkitExtension: the three methods left uncovered,
plus simpleFoam's two.

``test_openfoam_eulerian_solver.py`` covers ``flowType``,
``blockMesh_setDomainHeight`` and the rejection path of
``blockMesh_setBoundFromFile``.  What was left is
``blockMesh_setBoundFromBounds``, ``blockMesh_setBoundFromFile``'s success
path, ``writeFieldInCase`` and ``IC_getHydrostaticPressure`` -- and, in the
sibling module, ``simpleFoam_toolkitExtension``'s constructor and its
``buildSpecializedField`` hook, which no test imported at all.

The workflow stand-in subclasses the real ``workflow_Eulerian`` (hermes is
installed, so that class is real and ``isinstance`` against it works) and
records the keyword arguments it is handed, which is the only observable
effect these setters have.  ``getObjFileBoundaries`` reads an OBJ through
FreeCAD, which is not installable here, so it is monkeypatched on the module
and asserted as a forwarding seam: which file name went in, which corners
came back out, and that they arrived at the workflow unchanged alongside
dx/dy/dz.

Three of the four uncovered methods turn out to be dead on arrival.

Bugs pinned here (each with a strict xfail for the intended behaviour and a
passing characterisation of what happens today):

* B287 ``blockMesh_setBoundFromBounds`` passes the *bare name*
  ``blockMesh_setBoundFromBounds`` to ``get_classMethod_logger`` instead of
  the string ``"blockMesh_setBoundFromBounds"``.  A method name is not in
  scope inside its own body, so the very first statement raises NameError --
  before the workflow-type check, so even the error path is unreachable.  The
  three sibling ``blockMesh_*``/``IC_*`` methods all pass a string.
* B288 ``writeFieldInCase`` calls
  ``self.toolkit.OFObjectHome.getField(...)``.  ``OFObjectHome`` has no
  ``getField``; the closest are ``getEmptyField`` and
  ``getEmptyFieldFromCase``.  It also ignores four of its own six parameters
  (``caseDirectory``, ``components``, ``xarrayData`` are never read) and
  carries two stacked docstrings.
* B289 ``IC_getHydrostaticPressure`` calls
  ``self.toolkit.OFObjectHome.getFieldFromCase(...)`` on its first line;
  ``OFObjectHome`` has no ``getFieldFromCase`` either (the reader is
  ``readFieldFromCase``), so the hydrostatic-pressure initial condition
  cannot be produced at all -- and ``hera-openFoam ... IC hydrostatic``
  (CLI.IC_hydrostaticPressure) is wired straight to it.

Rough edge characterised but not pinned: ``simpleFoam_toolkitExtension``'s
``buildSpecializedField`` is a declared-but-unimplemented hook -- its
docstring describes updating the solver name, the field names and the
turbulence node, and its body is ``pass``.  There is no one obviously
intended behaviour to xfail against, so it is only characterised as a no-op.
"""
import pytest

from hera.simulations.openFoam.OFWorkflow import workflow_Eulerian
from hera.simulations.openFoam.eulerian import abstractEulerianSolver as solverModule
from hera.simulations.openFoam.eulerian.abstractEulerianSolver import (
    absractEulerianSolver_toolkitExtension,
)
from hera.simulations.openFoam.eulerian.simpleFoam import simpleFoam_toolkitExtension
from hera.simulations.openFoam.preprocessOFObjects import OFObjectHome

_CORNERS = dict(minx=-1.0, maxx=11.0, miny=-2.0, maxy=12.0, minz=0.0, maxz=20.0)


class _RecordingWorkflow(workflow_Eulerian):
    """A real workflow_Eulerian that records what the setters hand it."""

    def __init__(self):
        self.boundaryCalls = []
        self.heightCalls = []

    def set_blockMesh_blockBoundaries(self, **kwargs):
        self.boundaryCalls.append(kwargs)

    def set_blockMesh_blockHeight(self, Z, dz):
        self.heightCalls.append(dict(Z=Z, dz=dz))


class _ToolkitWithARealObjectHome:
    """The shape openFoam/toolkit.py has: an OFObjectHome instance attribute."""

    FLOWTYPE_INCOMPRESSIBLE = "incompressible"

    def __init__(self):
        self.OFObjectHome = OFObjectHome()
        self.meshCalls = []

    def getMesh(self, caseDirectory):
        self.meshCalls.append(caseDirectory)
        raise AssertionError("getMesh must not be reached in these tests")


@pytest.fixture()
def toolkit():
    return _ToolkitWithARealObjectHome()


@pytest.fixture()
def extension(toolkit):
    return absractEulerianSolver_toolkitExtension(toolkit=toolkit, solverName="mySolver",
                                                  incompressible=True)


# ---------------------------------------------------------------------------
# blockMesh_setBoundFromBounds -- B287
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestBlockMeshSetBoundFromBounds:
    @pytest.mark.xfail(
        strict=True,
        reason="B287: the first statement is "
               "get_classMethod_logger(self, blockMesh_setBoundFromBounds) -- the "
               "method's own name, unquoted, which is not a name in the "
               "function's scope -- so every call raises NameError before the "
               "bounds ever reach the workflow.  The sibling methods pass the "
               "name as a string.  See the consolidated findings issue.",
    )
    def test_the_bounds_should_reach_the_workflow(self, extension):
        workflow = _RecordingWorkflow()
        extension.blockMesh_setBoundFromBounds(workflow, minx=0.0, maxx=10.0, miny=0.0,
                                               maxy=10.0, minz=0.0, maxz=5.0,
                                               dx=1.0, dy=1.0, dz=0.5)
        assert workflow.boundaryCalls == [dict(minx=0.0, maxx=10.0, miny=0.0, maxy=10.0,
                                               minz=0.0, maxz=5.0, dx=1.0, dy=1.0, dz=0.5)]

    def test_it_currently_raises_nameerror_on_its_own_method_name(self, extension):
        """Characterisation of B287."""
        workflow = _RecordingWorkflow()
        with pytest.raises(NameError, match="blockMesh_setBoundFromBounds"):
            extension.blockMesh_setBoundFromBounds(workflow, 0.0, 10.0, 0.0, 10.0, 0.0, 5.0,
                                                   1.0, 1.0, 0.5)
        assert workflow.boundaryCalls == []

    def test_even_a_wrong_workflow_type_dies_on_the_logger_first(self, extension):
        """Characterisation of B287: the ValueError guard is unreachable."""
        with pytest.raises(NameError, match="blockMesh_setBoundFromBounds"):
            extension.blockMesh_setBoundFromBounds(object(), 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                                                   0.1, 0.1, 0.1)


# ---------------------------------------------------------------------------
# blockMesh_setBoundFromFile -- the success path
# ---------------------------------------------------------------------------

@pytest.fixture()
def objBoundaries(monkeypatch):
    """Replace the FreeCAD-backed OBJ reader with a recorder."""
    asked = []

    def reader(fileName):
        asked.append(fileName)
        return dict(_CORNERS)

    monkeypatch.setattr(solverModule, "getObjFileBoundaries", reader)
    return asked


@pytest.mark.unit
class TestBlockMeshSetBoundFromFile:
    def test_the_object_file_is_read_and_its_corners_are_forwarded(self, extension,
                                                                   objBoundaries):
        workflow = _RecordingWorkflow()
        extension.blockMesh_setBoundFromFile(workflow, "/cases/terrain.obj",
                                             dx=2.0, dy=3.0, dz=4.0)
        assert objBoundaries == ["/cases/terrain.obj"]
        assert workflow.boundaryCalls == [dict(dx=2.0, dy=3.0, dz=4.0, **_CORNERS)]

    def test_the_spacings_are_not_derived_from_the_object(self, extension, objBoundaries):
        workflow = _RecordingWorkflow()
        extension.blockMesh_setBoundFromFile(workflow, "/cases/terrain.obj",
                                             dx=0.5, dy=0.5, dz=0.5)
        forwarded = workflow.boundaryCalls[0]
        assert (forwarded["dx"], forwarded["dy"], forwarded["dz"]) == (0.5, 0.5, 0.5)

    def test_a_non_eulerian_workflow_is_rejected_before_the_file_is_read(self, extension,
                                                                         objBoundaries):
        with pytest.raises(ValueError, match="not eulerian"):
            extension.blockMesh_setBoundFromFile(object(), "/cases/terrain.obj",
                                                 dx=1.0, dy=1.0, dz=1.0)
        assert objBoundaries == []


# ---------------------------------------------------------------------------
# writeFieldInCase -- B288
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestWriteFieldInCase:
    @pytest.mark.xfail(
        strict=True,
        reason="B288: the body is a single call to "
               "self.toolkit.OFObjectHome.getField(...), and OFObjectHome "
               "defines no getField -- its factories are getEmptyField, "
               "getEmptyFieldFromCase and readFieldFromCase -- so the method "
               "raises AttributeError for every input.  Its caseDirectory, "
               "components and xarrayData parameters are never read either.  "
               "See the consolidated findings issue.",
    )
    def test_it_should_return_a_field_object(self, extension):
        field = extension.writeFieldInCase(fieldName="p", caseDirectory="/cases/one",
                                           components=None, xarrayData=None)
        assert field.name == "p"

    def test_it_currently_asks_the_object_home_for_a_method_it_does_not_have(self, extension):
        """Characterisation of B288."""
        with pytest.raises(AttributeError, match="getField"):
            extension.writeFieldInCase(fieldName="p", caseDirectory="/cases/one",
                                       components=None, xarrayData=None)

    def test_the_object_home_offers_the_neighbouring_factories_instead(self, toolkit):
        """Characterisation of B288: what the name should probably have been."""
        assert not hasattr(toolkit.OFObjectHome, "getField")
        assert hasattr(toolkit.OFObjectHome, "getEmptyField")
        assert hasattr(toolkit.OFObjectHome, "getEmptyFieldFromCase")


# ---------------------------------------------------------------------------
# IC_getHydrostaticPressure -- B289
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestICGetHydrostaticPressure:
    @pytest.mark.xfail(
        strict=True,
        reason="B289: the first statement calls "
               "self.toolkit.OFObjectHome.getFieldFromCase(...), and "
               "OFObjectHome defines no getFieldFromCase -- the reader is "
               "readFieldFromCase -- so the hydrostatic initial condition "
               "cannot be built at all, and CLI.IC_hydrostaticPressure, which "
               "calls straight into it, cannot work.  See the consolidated "
               "findings issue.",
    )
    def test_it_should_return_a_pressure_field(self, extension):
        field = extension.IC_getHydrostaticPressure(caseDirectory="/cases/one")
        assert field.name == "p"

    def test_it_currently_asks_the_object_home_for_a_method_it_does_not_have(self, extension):
        """Characterisation of B289."""
        with pytest.raises(AttributeError, match="getFieldFromCase"):
            extension.IC_getHydrostaticPressure(caseDirectory="/cases/one")

    def test_the_mesh_is_never_reached(self, extension, toolkit):
        """Characterisation of B289: it dies on line one, so nothing else
        in the method -- not the mesh read, not the fixedValue patch loop --
        ever runs."""
        with pytest.raises(AttributeError):
            extension.IC_getHydrostaticPressure(caseDirectory="/cases/one")
        assert toolkit.meshCalls == []

    def test_the_ground_pressure_argument_cannot_change_the_outcome(self, extension):
        """Characterisation of B289."""
        with pytest.raises(AttributeError, match="getFieldFromCase"):
            extension.IC_getHydrostaticPressure(caseDirectory="/cases/one",
                                                groundPressure=90000)


# ---------------------------------------------------------------------------
# simpleFoam_toolkitExtension
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestSimpleFoamExtension:
    @pytest.fixture()
    def simpleFoam(self, toolkit):
        return simpleFoam_toolkitExtension(toolkit=toolkit)

    def test_it_declares_the_simplefoam_solver(self, simpleFoam):
        assert simpleFoam.solverName == "simpleFOAM"

    def test_it_is_incompressible_by_construction(self, simpleFoam):
        assert simpleFoam.incompressible is True

    def test_the_flow_type_therefore_comes_from_the_toolkit(self, simpleFoam, toolkit):
        assert simpleFoam.flowType == toolkit.FLOWTYPE_INCOMPRESSIBLE

    def test_the_toolkit_is_kept_as_given(self, simpleFoam, toolkit):
        assert simpleFoam.toolkit is toolkit

    def test_it_inherits_the_abstract_extension(self, simpleFoam):
        assert isinstance(simpleFoam, absractEulerianSolver_toolkitExtension)

    def test_the_cache_doctype_is_the_shared_eulerian_one(self, simpleFoam):
        assert simpleFoam.DOCTYPE_EULERIAN_CACHE == \
            absractEulerianSolver_toolkitExtension.DOCTYPE_EULERIAN_CACHE

    @pytest.mark.parametrize("turbulenceType", ["laminar", "kEpsilon", "LES"])
    def test_build_specialized_field_is_still_a_no_op(self, simpleFoam, turbulenceType):
        """Characterisation: the docstring describes updating the solver name,
        the field names and the turbulence node; the body is `pass`."""
        workflow = _RecordingWorkflow()
        assert simpleFoam.buildSpecializedField(workflow=workflow,
                                                turbulenceType=turbulenceType) is None
        assert workflow.boundaryCalls == []
        assert workflow.heightCalls == []
