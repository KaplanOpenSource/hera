"""hera/simulations/openFoam/CLI.py: the `hera-openFoam` command handlers.

Every handler in that module takes the argparse Namespace produced by
`hera/bin/hera-openFoam` and drives the OpenFOAM toolkit.  These tests
build the Namespace by hand -- mirroring the dests that script actually
declares -- and hand the handler a recording stand-in for the toolkit, so
what is under test is the CLI plumbing only: which toolkit call each
subcommand makes, which arguments it forwards, which branch it takes on
the flags, what it writes, what it prints, and what it refuses to do on
bad input.  Nothing here touches MongoDB, OpenFOAM, Luigi, Slurm, hermes,
FreeCAD or paraview.

Deliberately not tested:
  * the toolkit methods themselves (createEmptyCase, compareWorkflows,
    getCaseResults, createDispersionFlowField, ...) -- they are the
    simulation engine, and they are precisely the seam these tests stub
    out.  The stand-ins return canned objects, so the assertions are all
    about *how the handler used the toolkit*, never about numbers coming
    back from it;
  * `hermes.utils.workflowAssembly.handler_buildExecute`, which the two
    buildExecute handlers delegate to unchanged: hermes is not installable
    in CI (see hera/tests/unit/_stubs.py), so only the delegation is
    asserted, against a fake hermes module;
  * the argparse wiring in hera/bin/hera-openFoam, which is a script and
    not importable as a module.  Where a bug below is about a mismatch
    between handler and parser, the parser side is quoted from that script
    rather than executed;
  * the first, shadowed definitions of `foam_mesh_blockMesh`,
    `foam_mesh_setDomainHeight` and `IC_hydrostaticPressure` (CLI.py:259,
    263, 267).  They are unreachable dead code -- see B166 -- and no
    test can call them.

Bugs pinned here (each with a strict xfail for the intended behaviour plus
a passing characterisation of what actually happens):

  * B166: `foam_mesh_blockMesh` and `foam_mesh_setDomainHeight` are each
    defined twice.  The second definition -- a bare `pass` -- wins, so the
    `raise NotImplementedError` the author wrote is unreachable and
    `hera-openFoam <solver> blockMesh setBoundsFromFile ...` silently
    succeeds while doing nothing.
  * B167: `if 'projectName' not in arguments` never fires.  argparse
    declares `--projectName` with `default=None`, so the attribute always
    exists and the caseConfiguration.json fallback behind that test is
    dead code; `projectName=None` is forwarded to getToolkit instead.
    (`stochasticLagrangian_dispersion_create` in the same file gets this
    right: `('projectName' in arguments) and (arguments.projectName is not
    None)`.)
  * B168: `Foam_createEmpty` (and the live `IC_hydrostaticPressure`)
    catch *any* failure to load caseConfiguration.json, report
    "Configuration file ... not found", and overwrite the file with
    `{"projectName": null}` -- so a caseConfiguration.json that merely has
    a syntax error is destroyed rather than reported.
  * B169: `stochasticLagrangian_dispersionFlow_create` catches
    FileExistsError, logs it with `logger.error` and returns normally, so
    the command exits 0 having created nothing.
  * B170: `stochasticLagrangian_postProcess_toParquet` computes `cache`
    on both branches and then never forwards it to `getCaseResults`, whose
    default is `cache=True` -- the "directory, not in the DB" branch
    explicitly sets `cache = False` and is ignored.
  * B171: `..._postProcess_toParquet` / `..._postProcess_toVTK` leave
    `outputName` unbound when the dispersion name is neither in the DB nor
    a directory on disk, and die with UnboundLocalError instead of the
    error message the surrounding code clearly intends.
  * B172: `foam_solver_simulations_list --format pandas --file out`
    tries to write a DataFrame *object* to a text file (TypeError).
  * B173: `stochasticLagrangian_dispersionFlow_list` never assigns
    `output` in its `--format pandas` branch -- the default -- and then
    prints it: UnboundLocalError after the groups have been printed.
  * B174: the same handler calls `tk.workflowCompare(...)`, a method
    that exists nowhere in hera; the real name is `compareWorkflows`.
  * B175: `stochasticLagrangian_dispersion_create --overwrite` bypasses
    the DB lookup of the *flow field* as well as the removal of the case
    directory, so re-creating a dispersion whose flow field lives in the
    DB fails with "<name> not found!".
  * B176: `stochasticLagrangian_postProcess_toVTK` reads
    `arguments.outputdir`, but the parser declares `--outputDirectory`
    (dest `outputDirectory`).  The flag is silently ignored and the VTK
    files always land under the process cwd.
  * B177: `stochasticLagrangian_source_makeEscapedMassFile` asks for
    `toolkitHome.OF_LSM`, which is not a ToolkitHome constant
    (AttributeError), and passes the project name positionally into
    getToolkit's `filesDirectory` slot.  Its subparser
    (`createReleaseRateFile`) declares none of the seven arguments the
    handler reads, so the command cannot work at all.
  * B178: `foam_solver_template_buildExecute` and `hermes_buildExecute`
    turn a missing hermes into `warnings.warn(...)` and then call the name
    the failed import was supposed to bind: NameError.

Rough edges characterised but not pinned as bugs, because there is no one
obviously-intended behaviour to xfail against: the live
`IC_hydrostaticPressure` computes `simulationType` and reads the whole
mesh into `cellCenters` and uses neither; `foam_templates_node_list` is
wired to no parser at all and still requires `arguments.solver`;
`foam_solver_template_buildExecute` forces `arguments.force = True`,
overriding `--allowDuplicate`; `--parameters` is optional in the
dispersionFlow-create parser but mandatory in the handler; and
`--outputDirectory` is inert for `postProcess toParquet`.
"""
import json
import os
import sys
import types
from argparse import Namespace

import pandas
import pytest

from hera.simulations.openFoam.CLI import (
    Foam_createEmpty,
    Foam_parser_FieldDescription,
    IC_hydrostaticPressure,
    foam_BC,
    foam_IC,
    foam_mesh_blockMesh,
    foam_mesh_setDomainHeight,
    foam_snappyhexmesh_addobject,
    foam_snappyhexmesh_setLocationInDomain,
    foam_solver_simulations_list,
    foam_solver_template_buildExecute,
    foam_solver_template_create,
    foam_solver_templates_list,
    foam_templates_node_list,
    hermes_buildExecute,
    objects_createVerticesAndBoundary,
    stochasticLagrangian_dispersionFlow_create,
    stochasticLagrangian_dispersionFlow_list,
    stochasticLagrangian_dispersionFlow_writeEmptyTemplate,
    stochasticLagrangian_dispersion_create,
    stochasticLagrangian_postProcess_toParquet,
    stochasticLagrangian_postProcess_toVTK,
    stochasticLagrangian_source_cylinder,
    stochasticLagrangian_source_makeEscapedMassFile,
)

PROJECT = "UNIT_OPENFOAM_PROJECT"
SOLVER = "buoyantReactingFoam"


# ---------------------------------------------------------------------------
# stand-ins
# ---------------------------------------------------------------------------

class _Recorder:
    """Shared call log.  Sub-objects record into their parent's list."""

    def __init__(self, calls=None, prefix=""):
        self.calls = [] if calls is None else calls
        self._prefix = prefix

    def _record(self, _callName, *args, **kwargs):
        # The leading underscore keeps this from colliding with a recorded
        # keyword literally called "name".
        self.calls.append((self._prefix + _callName, args, kwargs))

    def names(self):
        return [call[0] for call in self.calls]

    def call(self, name):
        matching = [call for call in self.calls if call[0] == name]
        assert len(matching) == 1, f"{name} called {len(matching)} times: {self.calls}"
        return matching[0]


class _FakeStochasticLagrangian(_Recorder):
    """The `toolkit.stochasticLagrangian` extension the CLI reaches for."""

    def __init__(self, calls, caseResults=None, createFlowRaises=None):
        super().__init__(calls, "stochasticLagrangian.")
        self.caseResults = caseResults
        self.createFlowRaises = createFlowRaises

    def createDispersionFlowField(self, **kwargs):
        self._record("createDispersionFlowField", **kwargs)
        if self.createFlowRaises is not None:
            raise self.createFlowRaises

    def createAndLinkDispersionCaseDirectory(self, dispersionDirectory,
                                             dispersionFlowDirectory=None):
        self._record("createAndLinkDispersionCaseDirectory", dispersionDirectory,
                     dispersionFlowDirectory=dispersionFlowDirectory)

    def writeParticlePositionFile(self, **kwargs):
        self._record("writeParticlePositionFile", **kwargs)

    def getCaseResults(self, caseDescriptor, **kwargs):
        self._record("getCaseResults", caseDescriptor, **kwargs)
        return self.caseResults


class _FakePresentation(_Recorder):
    def toUnstructuredVTK(self, **kwargs):
        self._record("toUnstructuredVTK", **kwargs)


class _FakeSolverExtension(_Recorder):
    """`getattr(toolkit, arguments.solver)` -- e.g. tk.buoyantReactingFoam."""

    def __init__(self, calls, pressureField):
        super().__init__(calls, "solver.")
        self.pressureField = pressureField

    def IC_getHydrostaticPressure(self, caseDirectory, **kwargs):
        self._record("IC_getHydrostaticPressure", caseDirectory, **kwargs)
        return self.pressureField


class _FakePressureField(_Recorder):
    def writeToCase(self, caseDirectory, **kwargs):
        self._record("writeToCase", caseDirectory, **kwargs)


class _FakeDocument:
    """Stand-in for a stored-workflow document."""

    def __init__(self, workflowName="dispersion_0", groupName="dispersion",
                 resource="/data/dispersion_0"):
        self.desc = {"workflowName": workflowName, "groupName": groupName}
        self.resource = resource


class _FakeOFToolkit(_Recorder):
    """Records every call the CLI makes, and returns canned answers.

    Only the surface hera/simulations/openFoam/CLI.py actually touches is
    implemented, with the real OFToolkit's constant values, so that a
    handler reaching for something the real toolkit does not have (see
    B174) fails here the same way it fails in production.
    """

    FLOWTYPE_COMPRESSIBLE = "compressible"
    FLOWTYPE_INCOMPRESSIBLE = "incompressible"
    DOCTYPE_OF_FLOWDISPERSION = "flowDispersion"

    def __init__(self, projectName=PROJECT, templates=None, nodeTemplates=None,
                 comparison=None, solverDocuments=(), documents=(),
                 dataSourceData=None, caseResults=None,
                 createFlowRaises=None, solverName=SOLVER):
        super().__init__()
        self.projectName = projectName
        self.templates = pandas.DataFrame() if templates is None else templates
        self.nodeTemplates = pandas.DataFrame() if nodeTemplates is None else nodeTemplates
        self.comparison = comparison
        self.solverDocuments = list(solverDocuments)
        self.documents = list(documents)
        self.dataSourceData = {} if dataSourceData is None else dataSourceData
        self.stochasticLagrangian = _FakeStochasticLagrangian(
            self.calls, caseResults=caseResults, createFlowRaises=createFlowRaises
        )
        self.presentation = _FakePresentation(self.calls, "presentation.")
        self.pressureField = _FakePressureField(self.calls, "pressureField.")
        setattr(self, solverName, _FakeSolverExtension(self.calls, self.pressureField))

    # -- the toolkit surface the CLI uses ----------------------------------
    def createEmptyCase(self, caseDirectory, fieldList=None, flowType=None,
                        additionalFieldsDescription=None):
        self._record("createEmptyCase", caseDirectory=caseDirectory,
                     fieldList=fieldList, flowType=flowType,
                     additionalFieldsDescription=additionalFieldsDescription)

    def listHermesSolverTemplates(self, solverName):
        self._record("listHermesSolverTemplates", solverName)
        return self.templates

    def listHermesNodesTemplates(self):
        self._record("listHermesNodesTemplates")
        return self.nodeTemplates

    def getDataSourceData(self, datasourceName, **query):
        self._record("getDataSourceData", datasourceName, **query)
        return self.dataSourceData

    def getWorkflowListOfSolvers(self, solverName, **query):
        self._record("getWorkflowListOfSolvers", solverName, **query)
        return self.solverDocuments

    def compareWorkflows(self, Workflow, longFormat=False, transpose=False):
        self._record("compareWorkflows", Workflow, longFormat=longFormat,
                     transpose=transpose)
        return self.comparison

    def getWorkflowDocumentFromDB(self, nameOrWorkflow, doctype=None, **query):
        self._record("getWorkflowDocumentFromDB", nameOrWorkflow, doctype=doctype,
                     **query)
        return self.documents

    def getMesh(self, caseDirectory, **kwargs):
        self._record("getMesh", caseDirectory, **kwargs)
        return "mesh-cell-centers"


class _FakeOFToolkitWithWorkflowCompare(_FakeOFToolkit):
    """Adds the `workflowCompare` method the real toolkit lacks (B174).

    Needed to reach anything past CLI.py:330 in
    stochasticLagrangian_dispersionFlow_list.
    """

    def workflowCompare(self, workflowsType=None):
        self._record("workflowCompare", workflowsType=workflowsType)
        return self.comparison


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def install_toolkit(monkeypatch):
    """Make `toolkitHome.getToolkit` hand the CLI our stand-in.

    Returns a callable that installs the toolkit and yields the list of
    getToolkit calls, so a test can assert which toolkit of which project
    was requested -- and how many times.

    The patch is applied to the *class*: patching the singleton instance
    would leave a permanent instance attribute behind on teardown.
    """
    from hera import toolkitHome

    def _install(toolkit):
        asked = []

        def _getToolkit(self, toolkitName=None, filesDirectory=None, **kwargs):
            asked.append(dict(toolkitName=toolkitName,
                              filesDirectory=filesDirectory,
                              kwargs=kwargs))
            return toolkit

        monkeypatch.setattr(type(toolkitHome), "getToolkit", _getToolkit)
        return asked

    return _install


@pytest.fixture()
def in_tmp_cwd(tmp_path, monkeypatch):
    """Most handlers here read or write files relative to the process cwd."""
    monkeypatch.chdir(tmp_path)
    return tmp_path


@pytest.fixture()
def fake_hermes(monkeypatch):
    """Replace the whole `hermes` package with a recording fake.

    hermes is stubbed out in CI (_stubs.py) and a real local install must
    not be executed either, so the module chain is replaced wholesale.
    """
    recorded = []

    def _handler_buildExecute(arguments):
        recorded.append(arguments)

    hermes = types.ModuleType("hermes")
    hermes.__path__ = []
    utils = types.ModuleType("hermes.utils")
    utils.__path__ = []
    assembly = types.ModuleType("hermes.utils.workflowAssembly")
    assembly.handler_buildExecute = _handler_buildExecute

    monkeypatch.setitem(sys.modules, "hermes", hermes)
    monkeypatch.setitem(sys.modules, "hermes.utils", utils)
    monkeypatch.setitem(sys.modules, "hermes.utils.workflowAssembly", assembly)
    return recorded


@pytest.fixture()
def no_hermes(monkeypatch):
    """Make `from hermes.utils.workflowAssembly import ...` raise ImportError."""
    monkeypatch.setitem(sys.modules, "hermes.utils.workflowAssembly", None)


@pytest.fixture()
def recorded_workflow_add(monkeypatch):
    """Intercept the `workflow_add` that buildExecute imports at call time."""
    import hera.simulations.CLI as simulationsCLI

    recorded = []
    monkeypatch.setattr(simulationsCLI, "workflow_add", recorded.append)
    return recorded


def _templatesFrame(names=("windLog", "hydro"), datasources=None):
    """The frame shape listHermesSolverTemplates really returns: the
    template name is the *index* (named "name"), and the datasource is a
    column, which is why the handler can `query("name==@tname")`."""
    datasources = [f"{name}_source" for name in names] if datasources is None else datasources
    return pandas.DataFrame(
        {"name": list(names), "datasourceName": list(datasources)}
    ).set_index("name")


# ---------------------------------------------------------------------------
# the two handlers that only write a file: no toolkit, no mocking
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestFoamParserFieldDescription:
    def test_every_requested_field_gets_a_dimensionless_template_entry(self, tmp_path):
        target = tmp_path / "fields.json"
        Foam_parser_FieldDescription(
            Namespace(fileName=str(target), fields=["U", "T"])
        )
        written = json.loads(target.read_text())
        assert sorted(written) == ["T", "U"]
        assert written["U"] == {"dimensions": "[0 0 0 0 0 0 0]", "componentNames": None}

    def test_no_fields_at_all_yields_a_single_example_entry(self, tmp_path):
        target = tmp_path / "fields.json"
        Foam_parser_FieldDescription(Namespace(fileName=str(target), fields=None))
        assert list(json.loads(target.read_text())) == ["exampleField"]

    def test_an_empty_field_list_also_yields_the_example_entry(self, tmp_path):
        """`--fields` with no values parses as [], not None."""
        target = tmp_path / "fields.json"
        Foam_parser_FieldDescription(Namespace(fileName=str(target), fields=[]))
        assert list(json.loads(target.read_text())) == ["exampleField"]

    def test_the_file_name_is_taken_verbatim_and_overwrites_what_is_there(self, tmp_path):
        target = tmp_path / "fields.json"
        target.write_text("PREVIOUS CONTENT")
        Foam_parser_FieldDescription(Namespace(fileName=str(target), fields=["U"]))
        assert list(json.loads(target.read_text())) == ["U"]

    def test_a_bare_file_name_is_written_relative_to_the_current_directory(self, in_tmp_cwd):
        Foam_parser_FieldDescription(Namespace(fileName="fields.json", fields=["U"]))
        assert (in_tmp_cwd / "fields.json").exists()


@pytest.mark.unit
class TestDispersionFlowWriteEmptyTemplate:
    def test_the_template_carries_the_documented_placeholder_skeleton(self, tmp_path):
        target = tmp_path / "flow.json"
        stochasticLagrangian_dispersionFlow_writeEmptyTemplate(
            Namespace(templateFile=str(target))
        )
        written = json.loads(target.read_text())
        assert sorted(written) == ["dispersionDuration", "dispersionFields", "originalFlow"]
        assert written["originalFlow"]["linkMeshSymbolically"] is True
        assert written["originalFlow"]["time"] == {
            "temporalType": "steadyState|dynamic",
            "timestep": "< time >",
        }
        assert written["dispersionFields"] == {"<field name>": "<constant value>"}

    def test_the_template_is_valid_json_written_where_asked(self, in_tmp_cwd):
        stochasticLagrangian_dispersionFlow_writeEmptyTemplate(
            Namespace(templateFile="myDispersion.json")
        )
        assert json.loads((in_tmp_cwd / "myDispersion.json").read_text())

    def test_an_existing_template_file_is_overwritten(self, tmp_path):
        target = tmp_path / "flow.json"
        target.write_text("PREVIOUS CONTENT")
        stochasticLagrangian_dispersionFlow_writeEmptyTemplate(
            Namespace(templateFile=str(target))
        )
        assert "originalFlow" in json.loads(target.read_text())


# ---------------------------------------------------------------------------
# the placeholder handlers
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestPlaceholderHandlers:
    """Five handlers whose whole body is `pass`, plus the two that shadow a
    NotImplementedError (B166)."""

    @pytest.mark.parametrize(
        "handler",
        [
            foam_snappyhexmesh_addobject,
            foam_snappyhexmesh_setLocationInDomain,
            foam_IC,
            foam_BC,
        ],
        ids=["snappyhexmesh_addobject", "snappyhexmesh_setLocationInDomain",
             "foam_IC", "foam_BC"],
    )
    def test_an_unimplemented_handler_accepts_anything_and_returns_none(self, handler):
        assert handler(Namespace(projectName=PROJECT, solver=SOLVER)) is None

    @pytest.mark.parametrize(
        "handler",
        [foam_mesh_blockMesh, foam_mesh_setDomainHeight],
        ids=["blockMesh", "setDomainHeight"],
    )
    @pytest.mark.xfail(
        strict=True,
        reason="B166: CLI.py defines foam_mesh_blockMesh (line 259) and "
               "foam_mesh_setDomainHeight (line 263) with a "
               "`raise NotImplementedError`, then defines both again at "
               "lines 640 and 654 with a bare `pass`.  The second "
               "definition wins, so the guard is unreachable and "
               "`hera-openFoam <solver> blockMesh setBoundsFromFile` -- "
               "which requires --dx --dy --dz --Z -- exits 0 having done "
               "nothing.  See the consolidated findings issue.",
    )
    def test_an_unimplemented_mesh_command_should_say_so(self, handler):
        with pytest.raises(NotImplementedError):
            handler(Namespace(templateName="t", fileName="f", dx=1, dy=1, dz=1,
                              Z=100, projectName=PROJECT, solver=SOLVER))

    @pytest.mark.parametrize(
        "handler",
        [foam_mesh_blockMesh, foam_mesh_setDomainHeight],
        ids=["blockMesh", "setDomainHeight"],
    )
    def test_a_mesh_command_currently_succeeds_silently(self, handler):
        """Characterisation of B166."""
        assert handler(Namespace(templateName="t", fileName="f", dx=1, dy=1,
                                 dz=1, Z=100, projectName=PROJECT,
                                 solver=SOLVER)) is None


# ---------------------------------------------------------------------------
# case createEmpty
# ---------------------------------------------------------------------------

def _createEmptyArgs(**overrides):
    arguments = dict(caseDirectory="myCase", fields=["U", "p"], projectName=PROJECT,
                     fieldsDescription=None, incompressible=True, solver=SOLVER)
    arguments.update(overrides)
    return Namespace(**arguments)


@pytest.mark.unit
class TestFoamCreateEmpty:
    def test_the_case_directory_fields_and_description_are_forwarded(self, install_toolkit):
        toolkit = _FakeOFToolkit()
        asked = install_toolkit(toolkit)
        Foam_createEmpty(_createEmptyArgs(fieldsDescription="extra.json"))
        assert asked[0]["toolkitName"] == "OpenFOAM"
        assert asked[0]["kwargs"]["projectName"] == PROJECT
        name, args, kwargs = toolkit.call("createEmptyCase")
        assert kwargs["caseDirectory"] == "myCase"
        assert kwargs["fieldList"] == ["U", "p"]
        assert kwargs["additionalFieldsDescription"] == "extra.json"

    def test_the_incompressible_flag_selects_the_incompressible_flow_type(self, install_toolkit):
        toolkit = _FakeOFToolkit()
        install_toolkit(toolkit)
        Foam_createEmpty(_createEmptyArgs(incompressible=True))
        assert toolkit.call("createEmptyCase")[2]["flowType"] == toolkit.FLOWTYPE_INCOMPRESSIBLE

    def test_without_the_incompressible_flag_the_flow_type_is_compressible(self, install_toolkit):
        toolkit = _FakeOFToolkit()
        install_toolkit(toolkit)
        Foam_createEmpty(_createEmptyArgs(incompressible=False))
        assert toolkit.call("createEmptyCase")[2]["flowType"] == toolkit.FLOWTYPE_COMPRESSIBLE

    def test_an_absent_project_name_attribute_is_read_from_the_case_configuration(
        self, install_toolkit, in_tmp_cwd
    ):
        (in_tmp_cwd / "caseConfiguration.json").write_text(
            json.dumps({"projectName": "FROM_CASE_CONFIG"})
        )
        arguments = _createEmptyArgs()
        del arguments.projectName
        asked = install_toolkit(_FakeOFToolkit())
        Foam_createEmpty(arguments)
        assert asked[0]["kwargs"]["projectName"] == "FROM_CASE_CONFIG"

    def test_a_named_configuration_file_is_honoured_when_present(
        self, install_toolkit, in_tmp_cwd
    ):
        (in_tmp_cwd / "other.json").write_text(json.dumps({"projectName": "FROM_OTHER"}))
        arguments = _createEmptyArgs(configurationFile="other.json")
        del arguments.projectName
        asked = install_toolkit(_FakeOFToolkit())
        Foam_createEmpty(arguments)
        assert asked[0]["kwargs"]["projectName"] == "FROM_OTHER"

    def test_a_missing_configuration_file_raises_and_leaves_a_skeleton_behind(
        self, install_toolkit, in_tmp_cwd
    ):
        arguments = _createEmptyArgs()
        del arguments.projectName
        install_toolkit(_FakeOFToolkit())
        with pytest.raises(ValueError, match="not found"):
            Foam_createEmpty(arguments)
        assert json.loads((in_tmp_cwd / "caseConfiguration.json").read_text()) == {
            "projectName": None
        }


@pytest.mark.unit
class TestFoamCreateEmptyProjectNameFallbackIsDead:
    @pytest.mark.xfail(
        strict=True,
        reason="B167: Foam_createEmpty guards its caseConfiguration.json "
               "fallback with `if 'projectName' not in arguments`, but "
               "hera/bin/hera-openFoam declares --projectName with "
               "default=None, so the attribute is always present and the "
               "fallback is unreachable: projectName=None is handed to "
               "getToolkit instead of being read from the case "
               "configuration.  The correct idiom is used a few functions "
               "down in the same file.  See the consolidated findings issue.",
    )
    def test_a_none_project_name_should_fall_back_to_the_case_configuration(
        self, install_toolkit, in_tmp_cwd
    ):
        (in_tmp_cwd / "caseConfiguration.json").write_text(
            json.dumps({"projectName": "FROM_CASE_CONFIG"})
        )
        asked = install_toolkit(_FakeOFToolkit())
        Foam_createEmpty(_createEmptyArgs(projectName=None))
        assert asked[0]["kwargs"]["projectName"] == "FROM_CASE_CONFIG"

    def test_a_none_project_name_is_currently_forwarded_verbatim(
        self, install_toolkit, in_tmp_cwd
    ):
        """Characterisation of B167."""
        (in_tmp_cwd / "caseConfiguration.json").write_text(
            json.dumps({"projectName": "FROM_CASE_CONFIG"})
        )
        asked = install_toolkit(_FakeOFToolkit())
        Foam_createEmpty(_createEmptyArgs(projectName=None))
        assert asked[0]["kwargs"]["projectName"] is None


@pytest.mark.unit
class TestFoamCreateEmptyClobbersABrokenConfiguration:
    @pytest.mark.xfail(
        strict=True,
        reason="B168: Foam_createEmpty wraps loadJSON in a bare `except:` "
               "and treats every failure as 'file not found', overwriting "
               "the existing caseConfiguration.json with "
               "{'projectName': null}.  A configuration file with a JSON "
               "syntax error -- or any other read error -- is therefore "
               "destroyed, and the message misreports why.  See the "
               "consolidated findings issue.",
    )
    def test_an_unparsable_configuration_file_should_not_be_overwritten(
        self, install_toolkit, in_tmp_cwd
    ):
        configuration = in_tmp_cwd / "caseConfiguration.json"
        configuration.write_text('{"projectName": "REAL_PROJECT",,}')
        arguments = _createEmptyArgs()
        del arguments.projectName
        install_toolkit(_FakeOFToolkit())
        with pytest.raises(ValueError):
            Foam_createEmpty(arguments)
        assert "REAL_PROJECT" in configuration.read_text()

    def test_an_unparsable_configuration_file_is_currently_destroyed(
        self, install_toolkit, in_tmp_cwd
    ):
        """Characterisation of B168."""
        configuration = in_tmp_cwd / "caseConfiguration.json"
        configuration.write_text('{"projectName": "REAL_PROJECT",,}')
        arguments = _createEmptyArgs()
        del arguments.projectName
        install_toolkit(_FakeOFToolkit())
        with pytest.raises(ValueError, match="not found"):
            Foam_createEmpty(arguments)
        assert json.loads(configuration.read_text()) == {"projectName": None}


# ---------------------------------------------------------------------------
# template listing / creation
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestFoamSolverTemplatesList:
    def test_the_solver_is_forwarded_and_the_templates_are_printed(
        self, install_toolkit, capsys
    ):
        toolkit = _FakeOFToolkit(templates=_templatesFrame())
        asked = install_toolkit(toolkit)
        foam_solver_templates_list(Namespace(projectName=PROJECT, solver=SOLVER))
        assert asked[0]["kwargs"]["projectName"] == PROJECT
        assert toolkit.call("listHermesSolverTemplates")[1] == (SOLVER,)
        printed = capsys.readouterr().out
        assert f"The templates for project {PROJECT} with solver {SOLVER}" in printed
        assert "windLog" in printed

    def test_the_title_is_underlined_to_its_own_width(self, install_toolkit, capsys):
        toolkit = _FakeOFToolkit(templates=_templatesFrame())
        install_toolkit(toolkit)
        foam_solver_templates_list(Namespace(projectName=PROJECT, solver=SOLVER))
        lines = [line for line in capsys.readouterr().out.splitlines() if line]
        title = f"The templates for project {PROJECT} with solver {SOLVER}"
        assert lines[0] == "-" * len(title)
        assert lines[1] == title
        assert lines[2] == "-" * len(title)

    def test_an_absent_project_name_attribute_becomes_a_none_project(self, install_toolkit):
        """Intentional here: since hera 2.13.2 the toolkit resolves a None
        project name from the case file itself."""
        toolkit = _FakeOFToolkit()
        asked = install_toolkit(toolkit)
        foam_solver_templates_list(Namespace(solver=SOLVER))
        assert asked[0]["kwargs"]["projectName"] is None

    def test_the_printed_project_name_comes_from_the_toolkit_not_the_arguments(
        self, install_toolkit, capsys
    ):
        toolkit = _FakeOFToolkit(projectName="RESOLVED_BY_TOOLKIT")
        install_toolkit(toolkit)
        foam_solver_templates_list(Namespace(projectName=None, solver=SOLVER))
        assert "RESOLVED_BY_TOOLKIT" in capsys.readouterr().out


@pytest.mark.unit
class TestFoamSolverTemplateCreate:
    def _arguments(self, **overrides):
        arguments = dict(projectName=PROJECT, solver=SOLVER, templateName="windLog",
                         groupName=None, projectPath=None)
        arguments.update(overrides)
        return Namespace(**arguments)

    def test_the_template_datasource_is_written_as_group_underscore_one(
        self, install_toolkit, tmp_path
    ):
        toolkit = _FakeOFToolkit(templates=_templatesFrame(),
                                 dataSourceData={"nodes": {"a": 1}})
        install_toolkit(toolkit)
        foam_solver_template_create(
            self._arguments(groupName="myGroup", projectPath=str(tmp_path))
        )
        target = tmp_path / "myGroup_1.json"
        assert json.loads(target.read_text()) == {"nodes": {"a": 1}}

    def test_the_datasource_name_is_looked_up_by_the_template_name(
        self, install_toolkit, tmp_path
    ):
        toolkit = _FakeOFToolkit(templates=_templatesFrame(),
                                 dataSourceData={"nodes": {}})
        install_toolkit(toolkit)
        foam_solver_template_create(self._arguments(projectPath=str(tmp_path)))
        assert toolkit.call("getDataSourceData")[1] == ("windLog_source",)

    def test_the_group_name_defaults_to_the_template_name(self, install_toolkit, tmp_path):
        toolkit = _FakeOFToolkit(templates=_templatesFrame(), dataSourceData={})
        install_toolkit(toolkit)
        foam_solver_template_create(self._arguments(projectPath=str(tmp_path)))
        assert (tmp_path / "windLog_1.json").exists()

    def test_a_missing_project_path_writes_into_the_current_directory(
        self, install_toolkit, in_tmp_cwd
    ):
        toolkit = _FakeOFToolkit(templates=_templatesFrame(), dataSourceData={})
        install_toolkit(toolkit)
        foam_solver_template_create(self._arguments(projectPath=None))
        assert (in_tmp_cwd / "windLog_1.json").exists()

    def test_an_unknown_template_is_refused_with_the_known_ones_listed(
        self, install_toolkit, in_tmp_cwd
    ):
        toolkit = _FakeOFToolkit(templates=_templatesFrame(), dataSourceData={})
        install_toolkit(toolkit)
        with pytest.raises(ValueError, match="notATemplate is not known"):
            foam_solver_template_create(self._arguments(templateName="notATemplate"))
        assert not list(in_tmp_cwd.iterdir())

    def test_an_unknown_template_is_refused_before_any_datasource_is_read(
        self, install_toolkit, in_tmp_cwd
    ):
        toolkit = _FakeOFToolkit(templates=_templatesFrame(), dataSourceData={})
        install_toolkit(toolkit)
        with pytest.raises(ValueError):
            foam_solver_template_create(self._arguments(templateName="notATemplate"))
        assert "getDataSourceData" not in toolkit.names()

    def test_an_absent_project_name_attribute_becomes_a_none_project(
        self, install_toolkit, in_tmp_cwd
    ):
        toolkit = _FakeOFToolkit(templates=_templatesFrame(), dataSourceData={})
        asked = install_toolkit(toolkit)
        arguments = self._arguments()
        del arguments.projectName
        foam_solver_template_create(arguments)
        assert asked[0]["kwargs"]["projectName"] is None

    def test_a_project_without_any_templates_refuses_every_name(
        self, install_toolkit, in_tmp_cwd
    ):
        """listHermesSolverTemplates returns an empty frame when nothing is
        registered, and an empty index rejects any template name."""
        toolkit = _FakeOFToolkit(templates=pandas.DataFrame())
        install_toolkit(toolkit)
        with pytest.raises(ValueError, match="is not known"):
            foam_solver_template_create(self._arguments())


@pytest.mark.unit
class TestFoamTemplatesNodeList:
    def test_the_node_templates_are_printed_under_a_title(self, install_toolkit, capsys):
        nodes = pandas.DataFrame({"nodeName": ["copyDirectory"], "component": ["Node"]}).set_index("nodeName")
        toolkit = _FakeOFToolkit(nodeTemplates=nodes)
        install_toolkit(toolkit)
        foam_templates_node_list(Namespace(projectName=PROJECT, solver=SOLVER))
        printed = capsys.readouterr().out
        assert f"The templates of nodes for project {PROJECT}" in printed
        assert "copyDirectory" in printed
        assert toolkit.names().count("listHermesNodesTemplates") == 1

    def test_an_absent_project_name_attribute_becomes_a_none_project(self, install_toolkit):
        toolkit = _FakeOFToolkit()
        asked = install_toolkit(toolkit)
        foam_templates_node_list(Namespace(solver=SOLVER))
        assert asked[0]["kwargs"]["projectName"] is None

    def test_it_still_demands_a_solver_it_does_not_use(self, install_toolkit):
        """Characterisation: listHermesNodesTemplates takes no solver, but
        the title interpolates arguments.solver, so a Namespace without it
        crashes.  No parser in hera/bin/hera-openFoam reaches this handler
        (the line that would is commented out), so the attribute is not
        guaranteed by anything."""
        install_toolkit(_FakeOFToolkit())
        with pytest.raises(AttributeError, match="solver"):
            foam_templates_node_list(Namespace(projectName=PROJECT))


# ---------------------------------------------------------------------------
# simulation listing
# ---------------------------------------------------------------------------

def _listArgs(**overrides):
    arguments = dict(projectName=PROJECT, solver=SOLVER, longFormat=False,
                     transpose=False, format="pandas", file=None)
    arguments.update(overrides)
    return Namespace(**arguments)


def _comparison():
    return pandas.DataFrame({"parameter": ["dt"], "wf_0": [0.1], "wf_1": [0.2]})


@pytest.mark.unit
class TestFoamSolverSimulationsList:
    def test_no_stored_cases_reports_the_solver_and_the_project(
        self, install_toolkit, capsys
    ):
        install_toolkit(_FakeOFToolkit(solverDocuments=[]))
        foam_solver_simulations_list(_listArgs())
        assert f"There are no cases for {SOLVER} in {PROJECT}" in capsys.readouterr().out

    def test_the_group_of_every_document_is_compared_with_the_table_flags(
        self, install_toolkit, capsys
    ):
        toolkit = _FakeOFToolkit(solverDocuments=[_FakeDocument(groupName="grp")],
                                 comparison=_comparison())
        install_toolkit(toolkit)
        foam_solver_simulations_list(_listArgs(longFormat=True, transpose=True))
        name, args, kwargs = toolkit.call("compareWorkflows")
        assert args == (["grp"],)
        assert kwargs == dict(longFormat=True, transpose=True)
        assert "Group name: grp" in capsys.readouterr().out

    def test_every_distinct_group_is_compared_once(self, install_toolkit):
        documents = [_FakeDocument(groupName="a"), _FakeDocument(groupName="a"),
                     _FakeDocument(groupName="b")]
        toolkit = _FakeOFToolkit(solverDocuments=documents, comparison=_comparison())
        install_toolkit(toolkit)
        foam_solver_simulations_list(_listArgs())
        compared = [call[1][0] for call in toolkit.calls if call[0] == "compareWorkflows"]
        assert sorted(group for (group,) in compared) == ["a", "b"]

    def test_the_solver_is_forwarded_to_the_document_query(self, install_toolkit):
        toolkit = _FakeOFToolkit(solverDocuments=[])
        install_toolkit(toolkit)
        foam_solver_simulations_list(_listArgs())
        assert toolkit.call("getWorkflowListOfSolvers")[1] == (SOLVER,)

    def test_the_json_format_prints_reindented_json(self, install_toolkit, capsys):
        toolkit = _FakeOFToolkit(solverDocuments=[_FakeDocument(groupName="grp")],
                                 comparison=_comparison())
        install_toolkit(toolkit)
        foam_solver_simulations_list(_listArgs(format="json"))
        printed = capsys.readouterr().out
        payload = printed[printed.index("{"):printed.rindex("}") + 1]
        assert json.loads(payload)["parameter"] == {"0": "dt"}

    @pytest.mark.parametrize(
        "outputFormat,extension",
        [("json", "json"), ("latex", "tex"), ("csv", "csv")],
    )
    def test_a_file_without_a_suffix_gets_the_one_that_matches_the_format(
        self, install_toolkit, in_tmp_cwd, outputFormat, extension
    ):
        toolkit = _FakeOFToolkit(solverDocuments=[_FakeDocument(groupName="grp")],
                                 comparison=_comparison())
        install_toolkit(toolkit)
        foam_solver_simulations_list(_listArgs(format=outputFormat, file="out"))
        assert (in_tmp_cwd / f"out.{extension}").exists()

    def test_a_file_name_that_already_has_a_suffix_is_used_as_is(
        self, install_toolkit, in_tmp_cwd
    ):
        toolkit = _FakeOFToolkit(solverDocuments=[_FakeDocument(groupName="grp")],
                                 comparison=_comparison())
        install_toolkit(toolkit)
        foam_solver_simulations_list(_listArgs(format="csv", file="out.data"))
        assert (in_tmp_cwd / "out.data").exists()
        assert "dt" in (in_tmp_cwd / "out.data").read_text()

    def test_nothing_is_written_when_no_file_was_requested(
        self, install_toolkit, in_tmp_cwd
    ):
        toolkit = _FakeOFToolkit(solverDocuments=[_FakeDocument(groupName="grp")],
                                 comparison=_comparison())
        install_toolkit(toolkit)
        foam_solver_simulations_list(_listArgs(format="csv", file=None))
        assert not list(in_tmp_cwd.iterdir())

    def test_an_empty_comparison_says_so_and_writes_no_file(
        self, install_toolkit, in_tmp_cwd, capsys
    ):
        toolkit = _FakeOFToolkit(solverDocuments=[_FakeDocument(groupName="grp")],
                                 comparison=pandas.DataFrame())
        install_toolkit(toolkit)
        foam_solver_simulations_list(_listArgs(format="csv", file="out"))
        assert "Could not found any workflows to compare" in capsys.readouterr().out
        assert not list(in_tmp_cwd.iterdir())

    def test_the_message_for_an_empty_project_reports_the_raw_argument(
        self, install_toolkit, capsys
    ):
        """Characterisation of B167's knock-on effect: the handler prints
        `projectName`, which is whatever came off the command line -- so a
        default --projectName prints the word None rather than the project
        the toolkit actually resolved."""
        install_toolkit(_FakeOFToolkit(projectName="RESOLVED", solverDocuments=[]))
        foam_solver_simulations_list(_listArgs(projectName=None))
        assert f"There are no cases for {SOLVER} in None" in capsys.readouterr().out


@pytest.mark.unit
class TestFoamSolverSimulationsListPandasToFile:
    @pytest.mark.xfail(
        strict=True,
        reason="B172: with the default --format pandas, "
               "foam_solver_simulations_list sets `output = res` -- the "
               "DataFrame itself -- and then, if --file was given, calls "
               "outputFile.write(output) on a text-mode file, which raises "
               "TypeError.  Every other format converts to a string first. "
               "See the consolidated findings issue.",
    )
    def test_the_default_format_should_be_writable_to_a_file(
        self, install_toolkit, in_tmp_cwd
    ):
        toolkit = _FakeOFToolkit(solverDocuments=[_FakeDocument(groupName="grp")],
                                 comparison=_comparison())
        install_toolkit(toolkit)
        foam_solver_simulations_list(_listArgs(format="pandas", file="out"))
        assert (in_tmp_cwd / "out.txt").exists()

    def test_the_default_format_with_a_file_currently_raises_a_type_error(
        self, install_toolkit, in_tmp_cwd
    ):
        """Characterisation of B172."""
        toolkit = _FakeOFToolkit(solverDocuments=[_FakeDocument(groupName="grp")],
                                 comparison=_comparison())
        install_toolkit(toolkit)
        with pytest.raises(TypeError):
            foam_solver_simulations_list(_listArgs(format="pandas", file="out"))


# ---------------------------------------------------------------------------
# dispersionFlow create
# ---------------------------------------------------------------------------

def _dispersionFlowArgs(**overrides):
    arguments = dict(projectName=PROJECT, OriginalFlowField="baseFlow",
                     dispersionFlowName="myDispersionFlow",
                     dispersionFlowParams='{"dispersionFields": {"C": 0}}',
                     dispersionDuration="3600", overwrite=False)
    arguments.update(overrides)
    return Namespace(**arguments)


@pytest.mark.unit
class TestDispersionFlowCreate:
    def test_the_flow_name_source_flow_duration_and_overwrite_are_forwarded(
        self, install_toolkit
    ):
        toolkit = _FakeOFToolkit()
        asked = install_toolkit(toolkit)
        stochasticLagrangian_dispersionFlow_create(_dispersionFlowArgs(overwrite=True))
        assert asked[0]["kwargs"]["projectName"] == PROJECT
        name, args, kwargs = toolkit.call("stochasticLagrangian.createDispersionFlowField")
        assert kwargs["flowName"] == "myDispersionFlow"
        assert kwargs["OriginalFlowField"] == "baseFlow"
        assert kwargs["dispersionDuration"] == "3600"
        assert kwargs["overwrite"] is True

    def test_the_parameters_string_is_parsed_into_the_flow_data_dictionary(
        self, install_toolkit
    ):
        toolkit = _FakeOFToolkit()
        install_toolkit(toolkit)
        stochasticLagrangian_dispersionFlow_create(_dispersionFlowArgs())
        kwargs = toolkit.call("stochasticLagrangian.createDispersionFlowField")[2]
        assert kwargs["flowData"] == {"dispersionFields": {"C": 0}}

    def test_python_style_parameters_from_the_shell_are_translated_to_json(
        self, install_toolkit
    ):
        """Shell quoting makes single quotes and Python literals the norm,
        so the handler rewrites ' -> ", True -> true, None -> null."""
        toolkit = _FakeOFToolkit()
        install_toolkit(toolkit)
        stochasticLagrangian_dispersionFlow_create(
            _dispersionFlowArgs(
                dispersionFlowParams="{'linkMeshSymbolically': True, "
                                     "'timestep': None, 'steady': False}"
            )
        )
        kwargs = toolkit.call("stochasticLagrangian.createDispersionFlowField")[2]
        assert kwargs["flowData"] == {"linkMeshSymbolically": True,
                                      "timestep": None,
                                      "steady": False}

    def test_a_missing_parameters_flag_currently_raises_an_attribute_error(
        self, install_toolkit
    ):
        """Characterisation: `--parameters` is declared without required=True
        and without a default, so it arrives as None -- and the handler
        calls .replace() on it straight away."""
        install_toolkit(_FakeOFToolkit())
        with pytest.raises(AttributeError, match="replace"):
            stochasticLagrangian_dispersionFlow_create(
                _dispersionFlowArgs(dispersionFlowParams=None)
            )


@pytest.mark.unit
class TestDispersionFlowCreateSwallowsFileExists:
    @pytest.mark.xfail(
        strict=True,
        reason="B169: stochasticLagrangian_dispersionFlow_create catches "
               "the FileExistsError raised for an existing flow field, "
               "logs it with logger.error and then returns normally, so "
               "`hera-openFoam stochasticLagrangian dispersionFlow create` "
               "exits 0 having created nothing -- indistinguishable from "
               "success for any script that checks the exit status.  See "
               "the consolidated findings issue.",
    )
    def test_an_existing_flow_field_should_fail_the_command(self, install_toolkit):
        toolkit = _FakeOFToolkit(createFlowRaises=FileExistsError("already there"))
        install_toolkit(toolkit)
        with pytest.raises(Exception):
            stochasticLagrangian_dispersionFlow_create(_dispersionFlowArgs())

    def test_an_existing_flow_field_currently_returns_quietly(self, install_toolkit):
        """Characterisation of B169."""
        toolkit = _FakeOFToolkit(createFlowRaises=FileExistsError("already there"))
        install_toolkit(toolkit)
        assert stochasticLagrangian_dispersionFlow_create(_dispersionFlowArgs()) is None
        assert "stochasticLagrangian.createDispersionFlowField" in toolkit.names()

    def test_other_failures_are_not_swallowed(self, install_toolkit):
        """Only FileExistsError is caught -- a genuine failure propagates."""
        toolkit = _FakeOFToolkit(createFlowRaises=RuntimeError("boom"))
        install_toolkit(toolkit)
        with pytest.raises(RuntimeError, match="boom"):
            stochasticLagrangian_dispersionFlow_create(_dispersionFlowArgs())


# ---------------------------------------------------------------------------
# dispersionFlow list
# ---------------------------------------------------------------------------

def _dispersionFlowListArgs(**overrides):
    arguments = dict(projectName=PROJECT, format="json", file=None)
    arguments.update(overrides)
    return Namespace(**arguments)


def _groupedComparison(groups=("grp",)):
    return pandas.DataFrame({"groupName": list(groups),
                             "parameter": ["dt"] * len(groups),
                             "value": [0.1] * len(groups)})


@pytest.mark.unit
class TestDispersionFlowList:
    def test_the_workflow_type_is_fixed_to_the_stochastic_lagrangian_solver(
        self, install_toolkit
    ):
        toolkit = _FakeOFToolkitWithWorkflowCompare(comparison=_groupedComparison())
        asked = install_toolkit(toolkit)
        stochasticLagrangian_dispersionFlow_list(_dispersionFlowListArgs())
        assert asked[0]["toolkitName"] == "OpenFOAM"
        assert toolkit.call("workflowCompare")[2] == {
            "workflowsType": "stochasticLagrangianSolver"
        }

    def test_the_json_format_prints_reindented_json(self, install_toolkit, capsys):
        toolkit = _FakeOFToolkitWithWorkflowCompare(comparison=_groupedComparison())
        install_toolkit(toolkit)
        stochasticLagrangian_dispersionFlow_list(_dispersionFlowListArgs(format="json"))
        printed = capsys.readouterr().out
        payload = printed[printed.index("{"):printed.rindex("}") + 1]
        assert json.loads(payload)["groupName"] == {"0": "grp"}

    @pytest.mark.parametrize(
        "outputFormat,extension",
        [("json", "json"), ("latex", "tex"), ("csv", "csv")],
    )
    def test_a_file_without_a_suffix_gets_the_one_that_matches_the_format(
        self, install_toolkit, in_tmp_cwd, outputFormat, extension
    ):
        toolkit = _FakeOFToolkitWithWorkflowCompare(comparison=_groupedComparison())
        install_toolkit(toolkit)
        stochasticLagrangian_dispersionFlow_list(
            _dispersionFlowListArgs(format=outputFormat, file="flows")
        )
        assert (in_tmp_cwd / f"flows.{extension}").exists()

    def test_an_empty_comparison_says_so_and_writes_nothing(
        self, install_toolkit, in_tmp_cwd, capsys
    ):
        toolkit = _FakeOFToolkitWithWorkflowCompare(
            comparison=pandas.DataFrame({"groupName": []})
        )
        install_toolkit(toolkit)
        stochasticLagrangian_dispersionFlow_list(
            _dispersionFlowListArgs(format="pandas", file="flows")
        )
        assert (f"Could not found any workflows to compare in project {PROJECT}"
                in capsys.readouterr().out)
        assert not list(in_tmp_cwd.iterdir())

    def test_an_absent_project_name_attribute_is_read_from_the_case_configuration(
        self, install_toolkit, in_tmp_cwd
    ):
        (in_tmp_cwd / "caseConfiguration.json").write_text(
            json.dumps({"projectName": "FROM_CASE_CONFIG"})
        )
        toolkit = _FakeOFToolkitWithWorkflowCompare(comparison=_groupedComparison())
        asked = install_toolkit(toolkit)
        arguments = _dispersionFlowListArgs()
        del arguments.projectName
        stochasticLagrangian_dispersionFlow_list(arguments)
        assert asked[0]["kwargs"]["projectName"] == "FROM_CASE_CONFIG"

    def test_a_none_project_name_is_forwarded_verbatim(self, install_toolkit, in_tmp_cwd):
        """Characterisation of B167 in this handler: the caseConfiguration
        fallback is guarded by `'projectName' not in arguments`, which the
        argparse default never satisfies."""
        (in_tmp_cwd / "caseConfiguration.json").write_text(
            json.dumps({"projectName": "FROM_CASE_CONFIG"})
        )
        toolkit = _FakeOFToolkitWithWorkflowCompare(comparison=_groupedComparison())
        asked = install_toolkit(toolkit)
        stochasticLagrangian_dispersionFlow_list(_dispersionFlowListArgs(projectName=None))
        assert asked[0]["kwargs"]["projectName"] is None


@pytest.mark.unit
class TestDispersionFlowListPandasFormat:
    @pytest.mark.xfail(
        strict=True,
        reason="B173: the `--format pandas` branch of "
               "stochasticLagrangian_dispersionFlow_list prints each group "
               "but never assigns `output`, unlike the other three "
               "branches.  Execution then falls through to `print(output)` "
               "and dies with UnboundLocalError -- and pandas is the "
               "default format, so the subcommand fails whenever it has "
               "anything to show.  See the consolidated findings issue.",
    )
    def test_the_default_format_should_list_the_flows_without_crashing(
        self, install_toolkit
    ):
        toolkit = _FakeOFToolkitWithWorkflowCompare(comparison=_groupedComparison())
        install_toolkit(toolkit)
        stochasticLagrangian_dispersionFlow_list(_dispersionFlowListArgs(format="pandas"))

    def test_the_default_format_currently_dies_after_printing_the_groups(
        self, install_toolkit, capsys
    ):
        """Characterisation of B173."""
        toolkit = _FakeOFToolkitWithWorkflowCompare(comparison=_groupedComparison())
        install_toolkit(toolkit)
        with pytest.raises(UnboundLocalError, match="output"):
            stochasticLagrangian_dispersionFlow_list(
                _dispersionFlowListArgs(format="pandas")
            )
        assert "Group name grp" in capsys.readouterr().out


@pytest.mark.unit
class TestDispersionFlowListCallsAMethodThatDoesNotExist:
    @pytest.mark.xfail(
        strict=True,
        reason="B174: stochasticLagrangian_dispersionFlow_list calls "
               "tk.workflowCompare(workflowsType=...).  No class in hera "
               "defines workflowCompare -- the OpenFOAM toolkit's method is "
               "compareWorkflows(Workflow, longFormat, transpose) -- so the "
               "subcommand cannot get past that line.  See the consolidated "
               "findings issue.",
    )
    def test_the_openfoam_toolkit_should_have_the_method_the_handler_calls(self):
        from hera.simulations.openFoam.toolkit import OFToolkit

        assert hasattr(OFToolkit, "workflowCompare")

    def test_the_openfoam_toolkit_only_has_compareworkflows(self):
        """Characterisation of B174."""
        from hera.simulations.openFoam.toolkit import OFToolkit

        assert not hasattr(OFToolkit, "workflowCompare")
        assert hasattr(OFToolkit, "compareWorkflows")

    def test_a_toolkit_with_the_real_surface_raises_an_attribute_error(
        self, install_toolkit
    ):
        """Characterisation of B174."""
        install_toolkit(_FakeOFToolkit(comparison=_groupedComparison()))
        with pytest.raises(AttributeError, match="workflowCompare"):
            stochasticLagrangian_dispersionFlow_list(_dispersionFlowListArgs())


# ---------------------------------------------------------------------------
# dispersion create
# ---------------------------------------------------------------------------

def _dispersionCreateArgs(**overrides):
    arguments = dict(projectName=PROJECT, dispersionName="myDispersion",
                     dispersionFlowField="myFlow", overwrite=False)
    arguments.update(overrides)
    return Namespace(**arguments)


@pytest.mark.unit
class TestDispersionCreate:
    def test_a_flow_field_found_in_the_db_is_linked_by_its_resource(
        self, install_toolkit, in_tmp_cwd
    ):
        toolkit = _FakeOFToolkit(documents=[_FakeDocument(resource="/data/flow_0")])
        install_toolkit(toolkit)
        stochasticLagrangian_dispersion_create(_dispersionCreateArgs())
        name, args, kwargs = toolkit.call(
            "stochasticLagrangian.createAndLinkDispersionCaseDirectory"
        )
        assert args == ("myDispersion",)
        assert kwargs == {"dispersionFlowDirectory": "/data/flow_0"}

    def test_the_flow_field_is_looked_up_as_a_dispersion_flow_document(
        self, install_toolkit, in_tmp_cwd
    ):
        toolkit = _FakeOFToolkit(documents=[_FakeDocument()])
        install_toolkit(toolkit)
        stochasticLagrangian_dispersion_create(_dispersionCreateArgs())
        name, args, kwargs = toolkit.call("getWorkflowDocumentFromDB")
        assert args == ("myFlow",)
        assert kwargs["doctype"] == "flowDispersion"

    def test_a_flow_field_absent_from_the_db_is_used_as_a_directory(
        self, install_toolkit, in_tmp_cwd
    ):
        flowDirectory = in_tmp_cwd / "flowOnDisk"
        flowDirectory.mkdir()
        toolkit = _FakeOFToolkit(documents=[])
        install_toolkit(toolkit)
        stochasticLagrangian_dispersion_create(
            _dispersionCreateArgs(dispersionFlowField="flowOnDisk")
        )
        kwargs = toolkit.call(
            "stochasticLagrangian.createAndLinkDispersionCaseDirectory"
        )[2]
        assert kwargs["dispersionFlowDirectory"] == str(flowDirectory)
        assert os.path.isabs(kwargs["dispersionFlowDirectory"])

    def test_a_flow_field_that_is_neither_in_the_db_nor_on_disk_is_refused(
        self, install_toolkit, in_tmp_cwd
    ):
        toolkit = _FakeOFToolkit(documents=[])
        install_toolkit(toolkit)
        with pytest.raises(ValueError, match="nowhere not found"):
            stochasticLagrangian_dispersion_create(
                _dispersionCreateArgs(dispersionFlowField="nowhere")
            )
        assert "stochasticLagrangian.createAndLinkDispersionCaseDirectory" \
            not in toolkit.names()

    def test_an_existing_case_directory_without_overwrite_is_refused(
        self, install_toolkit, in_tmp_cwd
    ):
        (in_tmp_cwd / "myDispersion").mkdir()
        toolkit = _FakeOFToolkit(documents=[_FakeDocument()])
        install_toolkit(toolkit)
        with pytest.raises(ValueError, match="Use --overwrite"):
            stochasticLagrangian_dispersion_create(_dispersionCreateArgs())
        assert (in_tmp_cwd / "myDispersion").exists()

    def test_overwrite_removes_the_existing_case_directory_before_recreating(
        self, install_toolkit, in_tmp_cwd
    ):
        existing = in_tmp_cwd / "myDispersion"
        existing.mkdir()
        (existing / "stale.txt").write_text("stale")
        (in_tmp_cwd / "flowOnDisk").mkdir()
        toolkit = _FakeOFToolkit(documents=[])
        install_toolkit(toolkit)
        stochasticLagrangian_dispersion_create(
            _dispersionCreateArgs(dispersionFlowField="flowOnDisk", overwrite=True)
        )
        assert not existing.exists()
        assert "stochasticLagrangian.createAndLinkDispersionCaseDirectory" in toolkit.names()

    def test_a_none_project_name_falls_back_to_the_case_configuration(
        self, install_toolkit, in_tmp_cwd
    ):
        """This handler gets the check right -- `'projectName' in arguments
        and arguments.projectName is not None` -- which is what B167 is
        missing elsewhere."""
        (in_tmp_cwd / "caseConfiguration.json").write_text(
            json.dumps({"projectName": "FROM_CASE_CONFIG"})
        )
        asked = install_toolkit(_FakeOFToolkit(documents=[_FakeDocument()]))
        stochasticLagrangian_dispersion_create(_dispersionCreateArgs(projectName=None))
        assert asked[0]["kwargs"]["projectName"] == "FROM_CASE_CONFIG"


@pytest.mark.unit
class TestDispersionCreateOverwriteIgnoresTheDatabase:
    @pytest.mark.xfail(
        strict=True,
        reason="B175: stochasticLagrangian_dispersion_create tests "
               "`if len(doc) == 0 or arguments.overwrite:` before deciding "
               "how to resolve the dispersion *flow field*, so --overwrite "
               "-- which is documented as overwriting the dispersion case "
               "directory -- also discards the DB lookup and insists the "
               "flow field be a directory on disk.  Recreating a dispersion "
               "whose flow field lives in the DB therefore fails with "
               "'<name> not found!'.  See the consolidated findings issue.",
    )
    def test_overwrite_should_still_use_a_flow_field_stored_in_the_db(
        self, install_toolkit, in_tmp_cwd
    ):
        toolkit = _FakeOFToolkit(documents=[_FakeDocument(resource="/data/flow_0")])
        install_toolkit(toolkit)
        stochasticLagrangian_dispersion_create(_dispersionCreateArgs(overwrite=True))
        kwargs = toolkit.call(
            "stochasticLagrangian.createAndLinkDispersionCaseDirectory"
        )[2]
        assert kwargs["dispersionFlowDirectory"] == "/data/flow_0"

    def test_overwrite_currently_rejects_a_flow_field_that_is_only_in_the_db(
        self, install_toolkit, in_tmp_cwd
    ):
        """Characterisation of B175."""
        toolkit = _FakeOFToolkit(documents=[_FakeDocument(resource="/data/flow_0")])
        install_toolkit(toolkit)
        with pytest.raises(ValueError, match="myFlow not found"):
            stochasticLagrangian_dispersion_create(_dispersionCreateArgs(overwrite=True))


# ---------------------------------------------------------------------------
# injector position file
# ---------------------------------------------------------------------------

def _cylinderArgs(**overrides):
    arguments = dict(projectName=PROJECT, dispersionName="myDispersion",
                     center=["10", "20", "30"], radius="5", height="2",
                     particles="1000")
    arguments.update(overrides)
    return Namespace(**arguments)


@pytest.mark.unit
class TestSourceCylinder:
    def test_the_geometry_is_converted_to_numbers_and_forwarded(self, install_toolkit):
        toolkit = _FakeOFToolkit()
        asked = install_toolkit(toolkit)
        stochasticLagrangian_source_cylinder(_cylinderArgs())
        assert asked[0]["kwargs"]["projectName"] == PROJECT
        name, args, kwargs = toolkit.call("stochasticLagrangian.writeParticlePositionFile")
        assert kwargs == dict(type="Cylinder", dispersionName="myDispersion",
                              x=10.0, y=20.0, z=30.0, radius=5.0, height=2.0,
                              nParticles=1000)

    def test_the_particle_count_is_an_integer_and_the_rest_are_floats(self, install_toolkit):
        toolkit = _FakeOFToolkit()
        install_toolkit(toolkit)
        stochasticLagrangian_source_cylinder(_cylinderArgs(particles="7", radius="0.5"))
        kwargs = toolkit.call("stochasticLagrangian.writeParticlePositionFile")[2]
        assert isinstance(kwargs["nParticles"], int)
        assert isinstance(kwargs["radius"], float)
        assert kwargs["radius"] == 0.5

    def test_a_none_project_name_falls_back_to_the_case_configuration(
        self, install_toolkit, in_tmp_cwd
    ):
        (in_tmp_cwd / "caseConfiguration.json").write_text(
            json.dumps({"projectName": "FROM_CASE_CONFIG"})
        )
        asked = install_toolkit(_FakeOFToolkit())
        stochasticLagrangian_source_cylinder(_cylinderArgs(projectName=None))
        assert asked[0]["kwargs"]["projectName"] == "FROM_CASE_CONFIG"

    def test_a_named_configuration_file_is_honoured(self, install_toolkit, in_tmp_cwd):
        (in_tmp_cwd / "other.json").write_text(json.dumps({"projectName": "FROM_OTHER"}))
        asked = install_toolkit(_FakeOFToolkit())
        stochasticLagrangian_source_cylinder(
            _cylinderArgs(projectName=None, configurationFile="other.json")
        )
        assert asked[0]["kwargs"]["projectName"] == "FROM_OTHER"

    def test_a_non_numeric_centre_is_rejected_before_the_toolkit_is_touched(
        self, install_toolkit
    ):
        toolkit = _FakeOFToolkit()
        install_toolkit(toolkit)
        with pytest.raises(ValueError):
            stochasticLagrangian_source_cylinder(_cylinderArgs(center=["a", "b", "c"]))
        assert toolkit.calls == []

    def test_a_centre_with_too_few_components_is_an_index_error(self, install_toolkit):
        """Characterisation: --center is `nargs="*"`, so argparse accepts
        any number of values and the handler indexes [0], [1], [2] blind."""
        install_toolkit(_FakeOFToolkit())
        with pytest.raises(IndexError):
            stochasticLagrangian_source_cylinder(_cylinderArgs(center=["10"]))


# ---------------------------------------------------------------------------
# escaped-mass file
# ---------------------------------------------------------------------------

class _FakeLSMAnalysis(_Recorder):
    def __init__(self, calls, data):
        super().__init__(calls, "analysis.")
        self.data = data

    def getMassFromLog(self, logFile=None, solver=None):
        self._record("getMassFromLog", logFile=logFile, solver=solver)
        return self.data


class _FakeLSMToolkit(_Recorder):
    def __init__(self, data):
        super().__init__()
        self.analysis = _FakeLSMAnalysis(self.calls, data)


def _massLog():
    return pandas.DataFrame({
        "time": [1.0, 2.0, 3.0],
        "mass": [0.0, 10.0, 25.0],
        "filterType": ["outlet"] * 3,
        "action": ["escaped"] * 3,
    })


def _massArgs(casePath, **overrides):
    arguments = dict(casePath=casePath, patch="outlet", massFileName=None,
                     dt=None, logFile="log.solver", solver="myFoam",
                     action="escaped")
    arguments.update(overrides)
    return Namespace(**arguments)


@pytest.mark.unit
class TestSourceMakeEscapedMassFile:
    @pytest.mark.xfail(
        strict=True,
        reason="B177: stochasticLagrangian_source_makeEscapedMassFile "
               "asks for toolkitHome.OF_LSM, which is not a ToolkitHome "
               "constant (only the internal _toolkits registry has an "
               "'OF_LSM' key), so the command dies with AttributeError "
               "before doing anything.  On top of that the project name is "
               "passed positionally -- getToolkit(toolkitName, "
               "filesDirectory, **kwargs) -- so 'tmpProject' would land in "
               "filesDirectory, and the createReleaseRateFile subparser "
               "declares none of the seven attributes the handler reads. "
               "See the consolidated findings issue.",
    )
    def test_the_toolkit_name_constant_the_handler_uses_should_exist(self):
        from hera import toolkitHome

        assert hasattr(toolkitHome, "OF_LSM")

    def test_the_handler_currently_dies_on_the_missing_constant(self, tmp_path):
        """Characterisation of B177."""
        with pytest.raises(AttributeError, match="OF_LSM"):
            stochasticLagrangian_source_makeEscapedMassFile(_massArgs(str(tmp_path)))

    def test_the_project_name_is_passed_where_the_files_directory_belongs(
        self, install_toolkit, monkeypatch, tmp_path
    ):
        """Characterisation of B177: with the missing constant supplied,
        "tmpProject" is consumed by getToolkit's `filesDirectory`
        parameter."""
        from hera import toolkitHome

        monkeypatch.setattr(type(toolkitHome), "OF_LSM", "OF_LSM", raising=False)
        (tmp_path / "constant").mkdir()
        asked = install_toolkit(_FakeLSMToolkit(_massLog()))
        stochasticLagrangian_source_makeEscapedMassFile(_massArgs(str(tmp_path)))
        assert asked[0]["toolkitName"] == "OF_LSM"
        assert asked[0]["filesDirectory"] == "tmpProject"
        assert "projectName" not in asked[0]["kwargs"]

    def test_the_mass_file_is_written_into_the_constant_directory_of_the_case(
        self, install_toolkit, monkeypatch, tmp_path
    ):
        from hera import toolkitHome

        monkeypatch.setattr(type(toolkitHome), "OF_LSM", "OF_LSM", raising=False)
        (tmp_path / "constant").mkdir()
        toolkit = _FakeLSMToolkit(_massLog())
        install_toolkit(toolkit)
        stochasticLagrangian_source_makeEscapedMassFile(_massArgs(str(tmp_path)))
        written = (tmp_path / "constant" / "outletMass").read_text()
        assert "class       scalarField;" in written
        assert written.rstrip().endswith(")")
        assert toolkit.call("analysis.getMassFromLog")[2] == dict(
            logFile="log.solver", solver="myFoam"
        )

    def test_an_explicit_mass_file_name_overrides_the_patch_derived_one(
        self, install_toolkit, monkeypatch, tmp_path, capsys
    ):
        from hera import toolkitHome

        monkeypatch.setattr(type(toolkitHome), "OF_LSM", "OF_LSM", raising=False)
        (tmp_path / "constant").mkdir()
        install_toolkit(_FakeLSMToolkit(_massLog()))
        stochasticLagrangian_source_makeEscapedMassFile(
            _massArgs(str(tmp_path), massFileName="customMass")
        )
        assert (tmp_path / "constant" / "customMass").exists()
        assert str(tmp_path / "constant" / "customMass") in capsys.readouterr().out

    def test_the_record_count_and_the_mass_differences_come_from_the_log(
        self, install_toolkit, monkeypatch, tmp_path
    ):
        from hera import toolkitHome

        monkeypatch.setattr(type(toolkitHome), "OF_LSM", "OF_LSM", raising=False)
        (tmp_path / "constant").mkdir()
        install_toolkit(_FakeLSMToolkit(_massLog()))
        stochasticLagrangian_source_makeEscapedMassFile(_massArgs(str(tmp_path)))
        body = (tmp_path / "constant" / "outletMass").read_text().split("//")[-1]
        values = [line for line in body.splitlines() if line.strip()
                  and line.strip() not in ("(", ")")]
        # 3 log rows -> a count of 3 and three diffMass values, the first 0.
        assert values[0] == "3"
        assert [float(value) for value in values[1:]] == [0.0, 10.0, 15.0]

    def test_only_the_requested_patch_and_action_are_kept(
        self, install_toolkit, monkeypatch, tmp_path
    ):
        from hera import toolkitHome

        monkeypatch.setattr(type(toolkitHome), "OF_LSM", "OF_LSM", raising=False)
        (tmp_path / "constant").mkdir()
        data = pandas.concat(
            [
                _massLog(),
                pandas.DataFrame({"time": [1.0], "mass": [99.0],
                                  "filterType": ["wall"], "action": ["escaped"]}),
            ],
            ignore_index=True,
        )
        install_toolkit(_FakeLSMToolkit(data))
        stochasticLagrangian_source_makeEscapedMassFile(_massArgs(str(tmp_path)))
        assert "99" not in (tmp_path / "constant" / "outletMass").read_text()


# ---------------------------------------------------------------------------
# postProcess toParquet
# ---------------------------------------------------------------------------

def _toParquetArgs(**overrides):
    arguments = dict(projectName=PROJECT, dispersionName="myDispersion",
                     overwrite=False, cloudName="kinematicCloud",
                     outputDirectory="VTK")
    arguments.update(overrides)
    return Namespace(**arguments)


@pytest.mark.unit
class TestPostProcessToParquet:
    def test_a_dispersion_found_in_the_db_is_read_by_its_workflow_name(
        self, install_toolkit
    ):
        toolkit = _FakeOFToolkit(documents=[_FakeDocument(workflowName="myDispersion_3")],
                                 caseResults="the-data")
        install_toolkit(toolkit)
        stochasticLagrangian_postProcess_toParquet(_toParquetArgs(overwrite=True))
        name, args, kwargs = toolkit.call("stochasticLagrangian.getCaseResults")
        assert args == ("myDispersion_3",)
        assert kwargs["withVelocity"] is True
        assert kwargs["withMass"] is True
        assert kwargs["overwrite"] is True
        assert kwargs["cloudName"] == "kinematicCloud"

    def test_the_cloud_name_is_forwarded(self, install_toolkit):
        toolkit = _FakeOFToolkit(documents=[_FakeDocument()], caseResults="data")
        install_toolkit(toolkit)
        stochasticLagrangian_postProcess_toParquet(_toParquetArgs(cloudName="myCloud"))
        assert toolkit.call("stochasticLagrangian.getCaseResults")[2]["cloudName"] == "myCloud"

    def test_a_dispersion_that_is_only_a_directory_is_read_by_its_path(
        self, install_toolkit, in_tmp_cwd
    ):
        (in_tmp_cwd / "caseOnDisk").mkdir()
        toolkit = _FakeOFToolkit(documents=[], caseResults="data")
        install_toolkit(toolkit)
        stochasticLagrangian_postProcess_toParquet(
            _toParquetArgs(dispersionName="caseOnDisk")
        )
        assert toolkit.call("stochasticLagrangian.getCaseResults")[1] == ("caseOnDisk",)

    def test_the_toolkit_is_fetched_twice_for_the_same_project(self, install_toolkit):
        """Characterisation: the handler calls getToolkit at CLI.py:527 and
        again at CLI.py:543, rebinding the same name to an identical
        toolkit."""
        toolkit = _FakeOFToolkit(documents=[_FakeDocument()], caseResults="data")
        asked = install_toolkit(toolkit)
        stochasticLagrangian_postProcess_toParquet(_toParquetArgs())
        assert len(asked) == 2
        assert asked[0] == asked[1]

    def test_a_none_project_name_falls_back_to_the_case_configuration(
        self, install_toolkit, in_tmp_cwd
    ):
        (in_tmp_cwd / "caseConfiguration.json").write_text(
            json.dumps({"projectName": "FROM_CASE_CONFIG"})
        )
        asked = install_toolkit(_FakeOFToolkit(documents=[_FakeDocument()],
                                               caseResults="data"))
        stochasticLagrangian_postProcess_toParquet(_toParquetArgs(projectName=None))
        assert asked[0]["kwargs"]["projectName"] == "FROM_CASE_CONFIG"


@pytest.mark.unit
class TestPostProcessToParquetDropsTheCacheFlag:
    @pytest.mark.xfail(
        strict=True,
        reason="B170: stochasticLagrangian_postProcess_toParquet sets "
               "`cache = False` for a case found only as a directory and "
               "`cache = True` for one found in the DB, then calls "
               "getCaseResults without passing cache at all -- and "
               "getCaseResults defaults to cache=True.  The directory "
               "branch's decision not to cache is silently discarded.  See "
               "the consolidated findings issue.",
    )
    def test_the_cache_decision_should_reach_the_toolkit(
        self, install_toolkit, in_tmp_cwd
    ):
        (in_tmp_cwd / "caseOnDisk").mkdir()
        toolkit = _FakeOFToolkit(documents=[], caseResults="data")
        install_toolkit(toolkit)
        stochasticLagrangian_postProcess_toParquet(
            _toParquetArgs(dispersionName="caseOnDisk")
        )
        assert toolkit.call("stochasticLagrangian.getCaseResults")[2]["cache"] is False

    def test_the_cache_decision_is_currently_computed_and_dropped(
        self, install_toolkit, in_tmp_cwd
    ):
        """Characterisation of B170."""
        (in_tmp_cwd / "caseOnDisk").mkdir()
        toolkit = _FakeOFToolkit(documents=[], caseResults="data")
        install_toolkit(toolkit)
        stochasticLagrangian_postProcess_toParquet(
            _toParquetArgs(dispersionName="caseOnDisk")
        )
        assert "cache" not in toolkit.call("stochasticLagrangian.getCaseResults")[2]


@pytest.mark.unit
class TestPostProcessUnknownDispersionName:
    @pytest.mark.parametrize(
        "handler",
        [stochasticLagrangian_postProcess_toParquet,
         stochasticLagrangian_postProcess_toVTK],
        ids=["toParquet", "toVTK"],
    )
    @pytest.mark.xfail(
        strict=True,
        reason="B171: both postProcess handlers assign outputName only "
               "inside `if os.path.isdir(...)` (and, in toParquet, `cache` "
               "with it).  A dispersion name that is neither a stored "
               "workflow nor a directory leaves the name unbound, so the "
               "handler raises UnboundLocalError instead of the 'not found' "
               "error the surrounding log messages clearly intend.  See the "
               "consolidated findings issue.",
    )
    def test_an_unknown_dispersion_name_should_be_reported_clearly(
        self, install_toolkit, in_tmp_cwd, handler
    ):
        install_toolkit(_FakeOFToolkit(documents=[], caseResults="data"))
        with pytest.raises(ValueError):
            handler(_toParquetArgs(dispersionName="nowhere"))

    @pytest.mark.parametrize(
        "handler",
        [stochasticLagrangian_postProcess_toParquet,
         stochasticLagrangian_postProcess_toVTK],
        ids=["toParquet", "toVTK"],
    )
    def test_an_unknown_dispersion_name_currently_raises_unboundlocalerror(
        self, install_toolkit, in_tmp_cwd, handler
    ):
        """Characterisation of B171."""
        install_toolkit(_FakeOFToolkit(documents=[], caseResults="data"))
        with pytest.raises(UnboundLocalError, match="outputName"):
            handler(_toParquetArgs(dispersionName="nowhere"))


# ---------------------------------------------------------------------------
# postProcess toVTK
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestPostProcessToVTK:
    def test_the_case_results_are_handed_to_the_vtk_writer(
        self, install_toolkit, in_tmp_cwd
    ):
        toolkit = _FakeOFToolkit(documents=[_FakeDocument(workflowName="myDispersion_3")],
                                 caseResults="the-data")
        install_toolkit(toolkit)
        stochasticLagrangian_postProcess_toVTK(_toParquetArgs(overwrite=True))
        name, args, kwargs = toolkit.call("presentation.toUnstructuredVTK")
        assert kwargs["data"] == "the-data"
        assert kwargs["filename"] == "kinematicCloud"
        assert kwargs["overwrite"] is True
        assert (kwargs["xcoord"], kwargs["ycoord"], kwargs["zcoord"]) == ("x", "y", "z")

    def test_the_results_are_read_without_caching(self, install_toolkit, in_tmp_cwd):
        toolkit = _FakeOFToolkit(documents=[_FakeDocument()], caseResults="data")
        install_toolkit(toolkit)
        stochasticLagrangian_postProcess_toVTK(_toParquetArgs())
        kwargs = toolkit.call("stochasticLagrangian.getCaseResults")[2]
        assert kwargs["cache"] is False
        assert kwargs["withVelocity"] is True
        assert kwargs["withMass"] is True

    def test_the_output_directory_is_created_and_named_after_the_workflow(
        self, install_toolkit, in_tmp_cwd
    ):
        toolkit = _FakeOFToolkit(documents=[_FakeDocument(workflowName="myDispersion_3")],
                                 caseResults="data")
        install_toolkit(toolkit)
        stochasticLagrangian_postProcess_toVTK(_toParquetArgs())
        expected = in_tmp_cwd / "VTK" / "myDispersion_3"
        assert expected.is_dir()
        assert toolkit.call("presentation.toUnstructuredVTK")[2]["outputdirectory"] \
            == str(expected)

    def test_an_existing_output_directory_is_reused(self, install_toolkit, in_tmp_cwd):
        (in_tmp_cwd / "VTK" / "dispersion_0").mkdir(parents=True)
        toolkit = _FakeOFToolkit(documents=[_FakeDocument()], caseResults="data")
        install_toolkit(toolkit)
        stochasticLagrangian_postProcess_toVTK(_toParquetArgs())
        assert "presentation.toUnstructuredVTK" in toolkit.names()

    def test_a_dispersion_that_is_only_a_directory_names_the_output_after_it(
        self, install_toolkit, in_tmp_cwd
    ):
        (in_tmp_cwd / "caseOnDisk").mkdir()
        toolkit = _FakeOFToolkit(documents=[], caseResults="data")
        install_toolkit(toolkit)
        stochasticLagrangian_postProcess_toVTK(
            _toParquetArgs(dispersionName="caseOnDisk")
        )
        assert (in_tmp_cwd / "VTK" / "caseOnDisk").is_dir()

    def test_the_attribute_the_handler_actually_reads_does_redirect_the_output(
        self, install_toolkit, in_tmp_cwd, tmp_path
    ):
        """The `outputdir` spelling works -- which is what makes B176 a
        typo rather than a missing feature."""
        elsewhere = tmp_path / "elsewhere"
        toolkit = _FakeOFToolkit(documents=[_FakeDocument()], caseResults="data")
        install_toolkit(toolkit)
        arguments = _toParquetArgs()
        arguments.outputdir = str(elsewhere)
        stochasticLagrangian_postProcess_toVTK(arguments)
        assert (elsewhere / "VTK" / "dispersion_0").is_dir()

    def test_a_none_output_directory_attribute_falls_back_to_the_cwd(
        self, install_toolkit, in_tmp_cwd
    ):
        toolkit = _FakeOFToolkit(documents=[_FakeDocument()], caseResults="data")
        install_toolkit(toolkit)
        arguments = _toParquetArgs()
        arguments.outputdir = None
        stochasticLagrangian_postProcess_toVTK(arguments)
        assert (in_tmp_cwd / "VTK" / "dispersion_0").is_dir()


@pytest.mark.unit
class TestPostProcessToVTKIgnoresTheOutputDirectoryFlag:
    @pytest.mark.xfail(
        strict=True,
        reason="B176: stochasticLagrangian_postProcess_toVTK reads "
               "`arguments.outputdir`, but hera/bin/hera-openFoam declares "
               "the flag as --outputDirectory (dest 'outputDirectory', "
               "default 'VTK').  The spellings never match, so the branch "
               "always falls through to os.getcwd() and --outputDirectory "
               "is silently ignored.  See the consolidated findings issue.",
    )
    def test_the_output_directory_flag_should_be_honoured(
        self, install_toolkit, in_tmp_cwd, tmp_path
    ):
        elsewhere = tmp_path / "elsewhere"
        toolkit = _FakeOFToolkit(documents=[_FakeDocument()], caseResults="data")
        install_toolkit(toolkit)
        stochasticLagrangian_postProcess_toVTK(
            _toParquetArgs(outputDirectory=str(elsewhere))
        )
        assert (elsewhere / "VTK" / "dispersion_0").is_dir()

    def test_the_output_currently_always_lands_under_the_working_directory(
        self, install_toolkit, in_tmp_cwd, tmp_path
    ):
        """Characterisation of B176."""
        elsewhere = tmp_path / "elsewhere"
        toolkit = _FakeOFToolkit(documents=[_FakeDocument()], caseResults="data")
        install_toolkit(toolkit)
        stochasticLagrangian_postProcess_toVTK(
            _toParquetArgs(outputDirectory=str(elsewhere))
        )
        assert not elsewhere.exists()
        assert (in_tmp_cwd / "VTK" / "dispersion_0").is_dir()


# ---------------------------------------------------------------------------
# obj/stl vertices and boundaries
# ---------------------------------------------------------------------------

class _FakeRegion:
    def __init__(self, name):
        self.Name = name


@pytest.fixture()
def fake_freecad(monkeypatch):
    """Install the FreeCAD/Mesh pair the handler imports, plus a canned
    bounding box.  FreeCAD is not installable in CI (_stubs.py replaces it
    with a MagicMock and does not provide `Mesh` at all)."""
    recorded = {}

    freecad = types.ModuleType("FreeCAD")
    mesh = types.ModuleType("Mesh")

    class _Document:
        def findObjects(self):
            return [_FakeRegion("inlet"), _FakeRegion("walls")]

    def _getDocument(name):
        recorded["document"] = name
        return _Document()

    def _open(fileName):
        recorded.setdefault("opened", []).append(fileName)

    freecad.getDocument = _getDocument
    mesh.open = _open

    monkeypatch.setitem(sys.modules, "FreeCAD", freecad)
    monkeypatch.setitem(sys.modules, "Mesh", mesh)

    import hera.utils.freeCAD as freeCADutils

    def _boundaries(fileName):
        recorded["boundaries"] = fileName
        return dict(XMin=0.0, XMax=10.0, YMin=0.0, YMax=20.0, ZMin=0.0, ZMax=30.0)

    monkeypatch.setattr(freeCADutils, "getObjFileBoundaries", _boundaries)
    return recorded


def _splitObjectOutput(printed):
    vertices, boundaries = printed.split("----- Boundary conditions -------------")
    vertexJSON = vertices[vertices.index("{"):vertices.rindex("}") + 1]
    boundaryJSON = boundaries[boundaries.index("{"):boundaries.rindex("}") + 1]
    return json.loads(vertexJSON), json.loads(boundaryJSON)


@pytest.mark.unit
class TestObjectsCreateVerticesAndBoundary:
    def test_without_freecad_the_command_refuses_to_run(self, tmp_path):
        """FreeCAD is absent from this environment (the unit stub provides a
        MagicMock for `FreeCAD` but nothing for `Mesh`), which is the branch
        every CI run takes."""
        with pytest.raises(ImportError, match="freecad is not installed"):
            objects_createVerticesAndBoundary(
                Namespace(objectFile=str(tmp_path / "shape.obj"), fields=["U"])
            )

    def test_the_bounding_box_is_printed_as_eight_outward_shifted_vertices(
        self, fake_freecad, capsys
    ):
        objects_createVerticesAndBoundary(
            Namespace(objectFile="shape.obj", fields=["U"])
        )
        vertices, _ = _splitObjectOutput(capsys.readouterr().out)
        assert len(vertices["vertices"]) == 8
        assert vertices["vertices"][0] == [-0.1, -0.1, -0.1]
        assert vertices["vertices"][2] == [10.1, 20.1, -0.1]
        assert vertices["vertices"][6] == [10.1, 20.1, 30.1]

    def test_every_field_gets_a_zerogradient_boundary_for_every_region(
        self, fake_freecad, capsys
    ):
        objects_createVerticesAndBoundary(
            Namespace(objectFile="shape.obj", fields=["U", "T"])
        )
        _, boundaries = _splitObjectOutput(capsys.readouterr().out)
        assert sorted(boundaries) == ["T", "U"]
        assert boundaries["U"]["boundaryField"] == {
            "inlet": {"type": "zeroGradient"},
            "walls": {"type": "zeroGradient"},
        }

    def test_the_object_file_is_opened_and_measured(self, fake_freecad, capsys):
        objects_createVerticesAndBoundary(
            Namespace(objectFile="shape.obj", fields=["U"])
        )
        capsys.readouterr()
        assert fake_freecad["opened"] == ["shape.obj"]
        assert fake_freecad["document"] == "Unnamed"
        assert fake_freecad["boundaries"] == "shape.obj"

    def test_no_fields_at_all_is_a_type_error(self, fake_freecad):
        """Characterisation: --fields has no default, so omitting it makes
        the handler iterate over None."""
        with pytest.raises(TypeError):
            objects_createVerticesAndBoundary(
                Namespace(objectFile="shape.obj", fields=None)
            )


# ---------------------------------------------------------------------------
# initial conditions
# ---------------------------------------------------------------------------

def _hydrostaticArgs(**overrides):
    arguments = dict(caseDirectory="myCase", startTime=0, solver=SOLVER,
                     incompressible=False, projectName=PROJECT)
    arguments.update(overrides)
    return Namespace(**arguments)


@pytest.mark.unit
class TestICHydrostaticPressure:
    def test_the_pressure_field_is_written_to_the_case_at_the_start_time(
        self, install_toolkit
    ):
        toolkit = _FakeOFToolkit()
        install_toolkit(toolkit)
        IC_hydrostaticPressure(_hydrostaticArgs(startTime="100"))
        assert toolkit.call("solver.IC_getHydrostaticPressure")[1] == ("myCase",)
        name, args, kwargs = toolkit.call("pressureField.writeToCase")
        assert args == ("myCase",)
        assert kwargs == {"timeOrLocation": "100"}

    def test_the_solver_extension_is_chosen_by_name(self, install_toolkit):
        toolkit = _FakeOFToolkit(solverName="pimpleFoam")
        install_toolkit(toolkit)
        IC_hydrostaticPressure(_hydrostaticArgs(solver="pimpleFoam"))
        assert "solver.IC_getHydrostaticPressure" in toolkit.names()

    def test_a_solver_the_toolkit_does_not_extend_is_an_attribute_error(
        self, install_toolkit
    ):
        toolkit = _FakeOFToolkit(solverName="pimpleFoam")
        install_toolkit(toolkit)
        with pytest.raises(AttributeError, match="notASolver"):
            IC_hydrostaticPressure(_hydrostaticArgs(solver="notASolver"))

    def test_the_mesh_is_read_and_the_flow_type_computed_but_neither_is_used(
        self, install_toolkit
    ):
        """Characterisation: the handler builds `simulationType` from
        --incompressible and reads the whole mesh into `cellCenters`, then
        passes neither on -- IC_getHydrostaticPressure takes only
        (caseDirectory, fieldName, groundPressure).  The getMesh call is
        pure wasted I/O."""
        toolkit = _FakeOFToolkit()
        install_toolkit(toolkit)
        IC_hydrostaticPressure(_hydrostaticArgs())
        assert toolkit.call("getMesh")[1] == ("myCase",)
        assert toolkit.call("solver.IC_getHydrostaticPressure")[2] == {}
        assert "flowType" not in toolkit.call("pressureField.writeToCase")[2]

    def test_an_absent_project_name_attribute_is_read_from_the_case_configuration(
        self, install_toolkit, in_tmp_cwd
    ):
        (in_tmp_cwd / "caseConfiguration.json").write_text(
            json.dumps({"projectName": "FROM_CASE_CONFIG"})
        )
        asked = install_toolkit(_FakeOFToolkit())
        arguments = _hydrostaticArgs()
        del arguments.projectName
        IC_hydrostaticPressure(arguments)
        assert asked[0]["kwargs"]["projectName"] == "FROM_CASE_CONFIG"

    def test_a_missing_configuration_file_raises_and_leaves_a_skeleton_behind(
        self, install_toolkit, in_tmp_cwd
    ):
        install_toolkit(_FakeOFToolkit())
        arguments = _hydrostaticArgs()
        del arguments.projectName
        with pytest.raises(ValueError, match="not found"):
            IC_hydrostaticPressure(arguments)
        assert json.loads((in_tmp_cwd / "caseConfiguration.json").read_text()) == {
            "projectName": None
        }

    def test_the_parser_never_declares_a_project_name_for_this_subcommand(self):
        """Characterisation of the other half of B167: the
        `IC hydrostaticPressure` subparser declares only caseDirectory and
        --startTime, so here the caseConfiguration branch is the *only*
        one that can run -- the same test spelled the same way is dead code
        in Foam_createEmpty, where --projectName does exist."""
        arguments = _hydrostaticArgs()
        del arguments.projectName
        assert "projectName" not in arguments
        assert "configurationFile" not in arguments


# ---------------------------------------------------------------------------
# buildExecute
# ---------------------------------------------------------------------------

def _buildExecuteArgs(**overrides):
    arguments = dict(workflow="myFlow.json", projectName=PROJECT, workflowGroup=None,
                     overwrite=False, assignName=False, force=False, noDB=False,
                     execute=True, scheduler="local", scheduler_host=None,
                     scheduler_port=None, dispatch_id=None)
    arguments.update(overrides)
    return Namespace(**arguments)


@pytest.mark.unit
class TestFoamSolverTemplateBuildExecute:
    def test_the_workflow_is_added_to_the_database_and_then_executed(
        self, fake_hermes, recorded_workflow_add
    ):
        arguments = _buildExecuteArgs(noDB=False)
        foam_solver_template_buildExecute(arguments)
        assert recorded_workflow_add == [arguments]
        assert fake_hermes == [arguments]

    def test_the_nodb_flag_skips_the_database_and_only_executes(
        self, fake_hermes, recorded_workflow_add
    ):
        arguments = _buildExecuteArgs(noDB=True)
        foam_solver_template_buildExecute(arguments)
        assert recorded_workflow_add == []
        assert fake_hermes == [arguments]

    def test_the_scheduler_options_are_passed_through_untouched(
        self, fake_hermes, recorded_workflow_add
    ):
        arguments = _buildExecuteArgs(scheduler="central", scheduler_host="head",
                                      scheduler_port=8082, dispatch_id="abc")
        foam_solver_template_buildExecute(arguments)
        forwarded = fake_hermes[0]
        assert forwarded.scheduler == "central"
        assert forwarded.scheduler_host == "head"
        assert forwarded.scheduler_port == 8082
        assert forwarded.dispatch_id == "abc"

    def test_the_allow_duplicate_flag_is_overridden_to_true(
        self, fake_hermes, recorded_workflow_add
    ):
        """Characterisation: --allowDuplicate sets dest `force`, and the
        handler assigns `arguments.force = True` unconditionally, so the
        flag cannot be turned off from the command line."""
        arguments = _buildExecuteArgs(force=False)
        foam_solver_template_buildExecute(arguments)
        assert arguments.force is True
        assert recorded_workflow_add[0].force is True

    def test_the_database_step_happens_before_the_execution(
        self, fake_hermes, monkeypatch
    ):
        import hera.simulations.CLI as simulationsCLI

        order = []
        monkeypatch.setattr(simulationsCLI, "workflow_add",
                            lambda arguments: order.append("add"))
        monkeypatch.setitem(
            sys.modules["hermes.utils.workflowAssembly"].__dict__,
            "handler_buildExecute",
            lambda arguments: order.append("execute"),
        )
        foam_solver_template_buildExecute(_buildExecuteArgs(noDB=False))
        assert order == ["add", "execute"]


@pytest.mark.unit
class TestHermesBuildExecute:
    def test_the_arguments_are_delegated_to_hermes_unchanged(self, fake_hermes):
        arguments = _buildExecuteArgs()
        hermes_buildExecute(arguments)
        assert fake_hermes == [arguments]

    def test_nothing_is_added_to_the_database(self, fake_hermes, recorded_workflow_add):
        """Unlike foam_solver_template_buildExecute, this handler is pure
        delegation -- no workflow_add, no `force` rewrite."""
        arguments = _buildExecuteArgs(force=False)
        hermes_buildExecute(arguments)
        assert recorded_workflow_add == []
        assert arguments.force is False


@pytest.mark.unit
class TestBuildExecuteWithoutHermes:
    @pytest.mark.parametrize(
        "handler",
        [foam_solver_template_buildExecute, hermes_buildExecute],
        ids=["foam_solver_template_buildExecute", "hermes_buildExecute"],
    )
    @pytest.mark.xfail(
        strict=True,
        reason="B178: both buildExecute handlers wrap the "
               "`from hermes.utils.workflowAssembly import "
               "handler_buildExecute` in a try/except ImportError that only "
               "calls warnings.warn, and then call handler_buildExecute "
               "anyway -- so a missing hermes surfaces as "
               "`NameError: name 'handler_buildExecute' is not defined` "
               "after a warning, instead of the ImportError the comment in "
               "the source says was replaced.  See the consolidated "
               "findings issue.",
    )
    def test_a_missing_hermes_should_be_reported_as_an_import_error(
        self, no_hermes, recorded_workflow_add, handler
    ):
        with pytest.raises(ImportError):
            handler(_buildExecuteArgs(noDB=True))

    @pytest.mark.parametrize(
        "handler",
        [foam_solver_template_buildExecute, hermes_buildExecute],
        ids=["foam_solver_template_buildExecute", "hermes_buildExecute"],
    )
    def test_a_missing_hermes_currently_warns_and_then_raises_a_name_error(
        self, no_hermes, recorded_workflow_add, handler
    ):
        """Characterisation of B178."""
        with pytest.warns(UserWarning, match="hermes is not installed"):
            with pytest.raises(NameError, match="handler_buildExecute"):
                handler(_buildExecuteArgs(noDB=True))
