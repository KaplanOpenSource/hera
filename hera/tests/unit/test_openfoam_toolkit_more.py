"""openFoam/toolkit.py: everything in OFToolkit (and its Analysis /
Presentation layers) that is not the VTK-pipeline cache
(``test_openfoam_toolkit_vtk_cache.py``) and not the pure filesystem
helpers ``processorList`` / ``read_points_file`` / ``getMeshExtent`` /
``getTimeList``'s single-processor branch
(``test_openfoam_toolkit_pure.py``).

Covered here: ``runOFSimulation``, ``prepareSlurmWorkflowExecution``,
``getHermesWorkflow_Flow``, ``getMesh``, ``getMeshFromName``,
``getMeshExtentFromName``, ``createEmptyCase``, ``writeEmptyField``,
``template_add``, ``getTimeList``'s decomposed branch, the five
``datasetToOF`` delegations, ``Analysis`` (``datalayer``,
``getVTKPipeline``, ``getFiltersDocuments``, ``executeVTKFiltersAndLoad``)
and ``Presentation`` (``to_paraview_CSV``, ``toUnstructuredVTK``,
``toStructuredVTK``, ``loadLagrangianDataParallel``).

Nothing here shells out: OpenFOAM (``foamJob``, ``Allrun``), Luigi and
Slurm are reached only through monkeypatched seams, and PyFoam / paraview /
evtk are the conftest stubs, so no assertion below depends on a value a
stub invented -- only on which arguments crossed the seam, which branch ran
and which files appeared under ``tmp_path``.

Deliberately not covered:

* ``getMesh``'s actual mesh values and ``caseGeometry`` /
  ``datasetToCaseFields`` / ``datasetToSetFieldsDict``'s own behaviour --
  they need foamlib (absent here, which is why
  ``test_openfoam_dataset2of.py`` skips wholesale) and a real OpenFOAM case
  written by ``postProcess -func writeCellCentres``. Only the toolkit's
  argument forwarding into ``datasetToOF`` is pinned here.
* ``xarrayToSetFieldsDictDomain`` -- already covered at the end of
  ``test_openfoam_dataset2of.py``.
* ``Presentation.toUnstructuredVTK``'s dask branch -- the pandas branch
  covers the same body; the dask one only differs in how the timesteps are
  enumerated.

Bugs pinned below (each with an xfail for the intended behaviour and a
passing characterisation of what happens today):

* B192 ``getMesh`` builds a broken ``foamJob`` argv.
* B193 ``getMeshFromName`` / ``getMeshExtentFromName`` drop
  ``readParallel`` and ``time``.
* B194 ``prepareSlurmWorkflowExecution`` cannot take a workflow dict.
* B195 ``prepareSlurmWorkflowExecution`` demands a variations *dict*,
  which ``JSONVariations`` cannot read.
* B196 ``Analysis.getFiltersDocuments``'s not-found guard tests the wrong
  thing.
* B197 ``Presentation.toUnstructuredVTK``'s overwrite guard checks a path
  evtk never writes.
* B198 ``Presentation.toStructuredVTK`` calls a method its analysis layer
  does not have.
* B199 ``Presentation.loadLagrangianDataParallel`` drops the first
  requested time step.
* B200 ``Presentation.loadLagrangianDataParallel`` cannot take numeric
  time steps.
"""
import json
import os
import subprocess
import types

import pandas
import pytest

from hera.simulations.openFoam import toolkit as ofToolkitModule
from hera.simulations.openFoam.preprocessOFObjects import datasetToOF
from hera.simulations.openFoam.preprocessOFObjects.OFObjectHome import OFObjectHome
from hera.simulations.openFoam.toolkit import OFToolkit

_REQUIRED_NODES = ["controlDict", "fvSolution", "fvSchemes", "fileWriter", "Parameters"]
_ALPHA_PATH = "$.workflow.nodes.Parameters.Execution.input_parameters.alpha"


def _workflowJSON(alpha=1.0, solver="simpleFoam"):
    """A minimal real hermes workflow -- the same shape as the one in
    test_hermes_workflow_toolkit_db.py, which is what addWorkflowToGroup
    and workflow_Eulerian accept."""
    nodes = {
        nodeName: {"type": "general.Parameters", "Execution": {"input_parameters": {}}}
        for nodeName in _REQUIRED_NODES
    }
    nodes["Parameters"]["Execution"]["input_parameters"] = {"alpha": alpha}
    return {"workflow": {"nodeList": list(_REQUIRED_NODES), "nodes": nodes, "solver": solver}}


@pytest.fixture()
def of(unit_toolkit_factory):
    from hera import toolkitHome

    return unit_toolkit_factory(toolkitHome.SIMULATIONS_OPENFOAM)


@pytest.fixture()
def workflowInDB(of):
    """One real workflow document in the in-memory DB, named flow_0000."""
    return of.addWorkflowToGroup(_workflowJSON(), "flow")


class _Recorder:
    """Records every call, returns a fixed value."""

    def __init__(self, returnValue=None):
        self.calls = []
        self.returnValue = returnValue

    def __call__(self, *args, **kwargs):
        self.calls.append((args, kwargs))
        return self.returnValue

    @property
    def lastKwargs(self):
        return self.calls[-1][1]

    @property
    def lastArgs(self):
        return self.calls[-1][0]


# ---------------------------------------------------------------------------
# runOFSimulation
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestRunOFSimulation:
    @pytest.fixture()
    def runSeams(self, of, monkeypatch, tmp_path):
        """Neutralise the two things that leave the process: the Luigi build
        and the ./Allrun of every case."""
        caseDirectory = tmp_path / "case"
        caseDirectory.mkdir()

        class Document:
            desc = dict(workflowName="flow_0000")
            resource = str(caseDirectory)

        execute = _Recorder()
        chdir = _Recorder()
        run = _Recorder()
        monkeypatch.setattr(OFToolkit, "executeWorkflowFromDB", lambda self, *a, **kw: execute(*a, **kw))
        monkeypatch.setattr(OFToolkit, "getWorkflowListDocumentFromDB",
                            lambda self, *a, **kw: [Document(), Document()])
        monkeypatch.setattr(os, "chdir", chdir)
        monkeypatch.setattr(subprocess, "run", run)
        return types.SimpleNamespace(execute=execute, chdir=chdir, run=run,
                                     caseDirectory=str(caseDirectory))

    def test_it_forwards_the_scheduler_selection_to_the_build_step(self, of, runSeams):
        of.runOFSimulation("flow_0000", scheduler="central", schedulerHost="head",
                           schedulerPort=8082, dispatch_id="deadbeef")
        assert runSeams.execute.lastArgs == ("flow_0000",)
        assert runSeams.execute.lastKwargs == dict(scheduler="central", schedulerHost="head",
                                                   schedulerPort=8082, dispatch_id="deadbeef")

    def test_the_local_scheduler_is_the_default(self, of, runSeams):
        of.runOFSimulation("flow_0000")
        assert runSeams.execute.lastKwargs["scheduler"] == "local"

    def test_it_runs_allrun_once_per_document_in_its_own_directory(self, of, runSeams):
        of.runOFSimulation("flow_0000")
        assert [call[0][0] for call in runSeams.chdir.calls] == [runSeams.caseDirectory] * 2
        assert [call[0][0] for call in runSeams.run.calls] == [["./Allrun"]] * 2

    def test_a_failing_allrun_is_not_swallowed(self, of, runSeams):
        of.runOFSimulation("flow_0000")
        assert all(call[1] == dict(check=True) for call in runSeams.run.calls)


# ---------------------------------------------------------------------------
# prepareSlurmWorkflowExecution
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestPrepareSlurmWorkflowExecution:
    @pytest.fixture()
    def slurmSeam(self, monkeypatch):
        recorder = _Recorder()
        monkeypatch.setattr(ofToolkitModule.slurm, "prepareSlurmScriptExecution", recorder)
        return recorder

    def test_every_variation_becomes_a_workflow_in_the_group(self, of, workflowInDB, slurmSeam):
        # The base workflow has alpha=1.0, so 2.0 and 3.0 are two new members.
        of.prepareSlurmWorkflowExecution("flow_0000", [{_ALPHA_PATH: [2, 3]}])
        names = [doc.desc["workflowName"] for doc in of.getWorkflowListDocumentFromDB("flow")]
        assert names == ["flow_0000", "flow_0001", "flow_0002"]

    def test_a_variation_equal_to_the_base_reuses_the_existing_workflow(self, of, workflowInDB, slurmSeam):
        of.prepareSlurmWorkflowExecution("flow_0000", [{_ALPHA_PATH: [1, 2]}])
        assert open(os.path.join(of.filesDirectory, "cases.txt")).read().split() == ["flow_0000",
                                                                                    "flow_0001"]

    def test_the_case_list_file_holds_one_line_per_variation(self, of, workflowInDB, slurmSeam):
        of.prepareSlurmWorkflowExecution("flow_0000", [{_ALPHA_PATH: [2, 3]}])
        caseListPath = os.path.join(of.filesDirectory, "cases.txt")
        assert open(caseListPath).read().split() == ["flow_0001", "flow_0002"]

    def test_the_case_list_file_name_is_honoured(self, of, workflowInDB, slurmSeam):
        of.prepareSlurmWorkflowExecution("flow_0000", [{_ALPHA_PATH: [1]}],
                                         caseListFileName="myCases.txt")
        assert os.path.exists(os.path.join(of.filesDirectory, "myCases.txt"))
        assert slurmSeam.lastKwargs["jobDirListFilePath"] == os.path.join(of.filesDirectory, "myCases.txt")

    def test_it_forwards_the_slurm_options(self, of, workflowInDB, slurmSeam):
        of.prepareSlurmWorkflowExecution("flow_0000", [{_ALPHA_PATH: [1]}],
                                         slurmExecutionFileName="submitThese.sh",
                                         allocateProcessorsPerRun=4, memoryInGB=16,
                                         jobName="myJob", exclusive=True)
        assert slurmSeam.lastKwargs["slurmExecutionFilePath"] == os.path.join(of.filesDirectory,
                                                                             "submitThese.sh")
        assert slurmSeam.lastKwargs["allocateProcessorsPerRun"] == 4
        assert slurmSeam.lastKwargs["memoryInGB"] == 16
        assert slurmSeam.lastKwargs["jobName"] == "myJob"
        assert slurmSeam.lastKwargs["exclusive"] is True

    def test_the_script_syncs_and_builds_every_case(self, of, workflowInDB, slurmSeam):
        of.prepareSlurmWorkflowExecution("flow_0000", [{_ALPHA_PATH: [1]}])
        script = slurmSeam.lastKwargs["script"]
        assert 'hera-workflows sync --force "$dir"' in script
        assert 'hera-workflows buildExecute "$dir"' in script

    def test_add_all_run_controls_the_allrun_block(self, of, workflowInDB, slurmSeam):
        of.prepareSlurmWorkflowExecution("flow_0000", [{_ALPHA_PATH: [1]}], addAllRun=True)
        assert "bash Allrun" in slurmSeam.lastKwargs["script"]
        of.prepareSlurmWorkflowExecution("flow_0000", [{_ALPHA_PATH: [2]}], addAllRun=False)
        assert "bash Allrun" not in slurmSeam.lastKwargs["script"]

    def test_the_variations_may_come_from_a_file(self, of, workflowInDB, slurmSeam, tmp_path):
        variationsFile = tmp_path / "variations.json"
        variationsFile.write_text(json.dumps([{_ALPHA_PATH: [7, 8]}]))
        of.prepareSlurmWorkflowExecution("flow_0000", str(variationsFile))
        assert open(os.path.join(of.filesDirectory, "cases.txt")).read().split() == ["flow_0001", "flow_0002"]

    @pytest.mark.xfail(
        strict=True,
        reason="B194: baseConfiguration is documented as 'basic hermes "
               "workflow to run' (a dict), but the local name `workflow` -- "
               "which the group name and the case-list message are read from "
               "-- is only bound in the `isinstance(baseConfiguration, str)` "
               "branch. A dict therefore always dies with UnboundLocalError. "
               "The same hole swallows an invalid type: the else branch only "
               "calls logger.error and falls through with no return/raise. "
               "See the consolidated findings issue.",
    )
    def test_a_workflow_dict_should_be_accepted_as_the_base_configuration(self, of, slurmSeam):
        of.prepareSlurmWorkflowExecution(_workflowJSON(), [{_ALPHA_PATH: [1]}])

    def test_a_workflow_dict_currently_raises_unboundlocalerror(self, of, slurmSeam):
        """Characterisation of B194."""
        with pytest.raises(UnboundLocalError, match="workflow"):
            of.prepareSlurmWorkflowExecution(_workflowJSON(), [{_ALPHA_PATH: [1]}])

    def test_an_unusable_base_configuration_type_is_only_logged(self, of, slurmSeam):
        """Characterisation of B194: the else branch calls logger.error and
        then falls straight through into JSONVariations, which fails on the
        unusable value instead of the method refusing it up front."""
        with pytest.raises(TypeError, match="not iterable"):
            of.prepareSlurmWorkflowExecution(12345, [{_ALPHA_PATH: [1]}])

    @pytest.mark.xfail(
        strict=True,
        reason="B195: the jsonVariations type check accepts a dict (and "
               "logs an error for anything else), but hera.utils.jsonutils."
               "JSONVariations iterates its second argument as a *list of "
               "variation groups* -- so the accepted dict is iterated as its "
               "keys and dies with AttributeError deep inside "
               "JSONvariationItem, while the list form that actually works "
               "is the one the method complains about. "
               "See the consolidated findings issue.",
    )
    def test_a_variations_dict_should_work(self, of, workflowInDB, slurmSeam):
        of.prepareSlurmWorkflowExecution("flow_0000", {_ALPHA_PATH: [1, 2]})

    def test_a_variations_dict_currently_raises(self, of, workflowInDB, slurmSeam):
        """Characterisation of B195."""
        with pytest.raises(AttributeError, match="items"):
            of.prepareSlurmWorkflowExecution("flow_0000", {_ALPHA_PATH: [1, 2]})


# ---------------------------------------------------------------------------
# getHermesWorkflow_Flow / template_add
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetHermesWorkflowFlow:
    def test_it_builds_an_eulerian_workflow_from_the_file(self, of, monkeypatch):
        built = _Recorder(returnValue="theWorkflow")
        monkeypatch.setattr(ofToolkitModule, "workflow_Eulerian", built)
        assert of.getHermesWorkflow_Flow("/some/workflow.json") == "theWorkflow"
        assert built.lastArgs == ("/some/workflow.json",)

    def test_it_builds_a_real_eulerian_workflow_from_a_workflow_dict(self, of):
        from hera.simulations.openFoam.OFWorkflow import workflow_Eulerian

        assert isinstance(of.getHermesWorkflow_Flow(_workflowJSON()), workflow_Eulerian)


@pytest.mark.unit
class TestTemplateAdd:
    def test_it_is_still_a_stub(self, of):
        """template_add is declared but its body is `pass` -- pinned so that
        giving it a body is a deliberate change, not an accident."""
        assert of.template_add("name", "objFile") is None


# ---------------------------------------------------------------------------
# getMesh
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetMesh:
    @pytest.fixture()
    def runSeam(self, monkeypatch):
        recorder = _Recorder()
        monkeypatch.setattr(subprocess, "run", recorder)
        return recorder

    @pytest.fixture()
    def readFieldSeam(self, monkeypatch):
        recorder = _Recorder(returnValue="theCellCenters")
        monkeypatch.setattr(OFObjectHome, "readFieldFromCase", lambda self, **kw: recorder(**kw))
        return recorder

    def test_existing_cell_centres_are_read_without_running_openfoam(self, of, tmp_path,
                                                                     runSeam, readFieldSeam):
        (tmp_path / "0").mkdir()
        (tmp_path / "0" / "C").write_text("")
        assert of.getMesh(str(tmp_path)) == "theCellCenters"
        assert runSeam.calls == []
        assert readFieldSeam.lastKwargs == dict(fieldName="cellCenters", flowType="incompressible",
                                                caseDirectory=str(tmp_path), timeStep=0, readParallel=True)

    def test_the_requested_time_selects_the_time_directory(self, of, tmp_path, runSeam, readFieldSeam):
        (tmp_path / "100").mkdir()
        (tmp_path / "100" / "C").write_text("")
        of.getMesh(str(tmp_path), time=100)
        assert runSeam.calls == []
        assert readFieldSeam.lastKwargs["timeStep"] == 100

    def test_a_decomposed_case_is_preferred_when_it_exists(self, of, tmp_path, runSeam, readFieldSeam):
        (tmp_path / "processor0" / "0").mkdir(parents=True)
        (tmp_path / "processor0" / "0" / "C").write_text("")
        of.getMesh(str(tmp_path))
        # The decomposed cell centres are there, so nothing has to be computed.
        assert runSeam.calls == []

    def test_the_decomposed_case_is_ignored_when_read_parallel_is_false(self, of, tmp_path,
                                                                        runSeam, readFieldSeam):
        (tmp_path / "processor0" / "0").mkdir(parents=True)
        (tmp_path / "processor0" / "0" / "C").write_text("")
        with pytest.raises(RuntimeError, match="writeCellCentres"):
            of.getMesh(str(tmp_path), readParallel=False)
        # It looked for <case>/0/C, did not find it, and ran foamJob without -parallel.
        assert "-parallel" not in runSeam.lastArgs[0]

    def test_a_decomposed_case_without_cell_centres_runs_foamjob_in_parallel(self, of, tmp_path, runSeam):
        (tmp_path / "processor0").mkdir()
        with pytest.raises(RuntimeError, match="writeCellCentres"):
            of.getMesh(str(tmp_path))
        assert runSeam.lastArgs[0][:2] == ["foamJob", "-parallel"]

    def test_it_raises_when_writecellcentres_produced_nothing(self, of, tmp_path, runSeam):
        with pytest.raises(RuntimeError, match="Check mesh"):
            of.getMesh(str(tmp_path))
        assert runSeam.lastArgs[0][0] == "foamJob"
        assert runSeam.lastKwargs == dict(check=False)

    @pytest.mark.xfail(
        strict=True,
        reason="B192: getMesh passes `f'-case {caseDirectory}'` as ONE "
               "element of the argv list handed to subprocess.run, instead "
               "of the two elements '-case' and the directory. Since the "
               "call uses a list (no shell), foamJob receives a single "
               "argument '-case /path' that it cannot parse -- and when the "
               "case IS the cwd it receives an empty-string argument "
               "instead. See the consolidated findings issue.",
    )
    def test_the_case_option_should_be_two_argv_elements(self, of, tmp_path, runSeam):
        with pytest.raises(RuntimeError):
            of.getMesh(str(tmp_path))
        assert "-case" in runSeam.lastArgs[0]

    def test_the_case_option_is_currently_one_joined_argv_element(self, of, tmp_path, runSeam):
        """Characterisation of B192."""
        with pytest.raises(RuntimeError):
            of.getMesh(str(tmp_path))
        assert f"-case {tmp_path}" in runSeam.lastArgs[0]

    def test_a_case_that_is_the_cwd_passes_an_empty_argument(self, of, tmp_path, runSeam, monkeypatch):
        """Characterisation of B192, second symptom."""
        monkeypatch.chdir(tmp_path)
        with pytest.raises(RuntimeError):
            of.getMesh(str(tmp_path))
        assert "" in runSeam.lastArgs[0]


# ---------------------------------------------------------------------------
# getMeshFromName / getMeshExtentFromName
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetMeshFromName:
    def test_an_unknown_name_gives_none(self, of):
        assert of.getMeshFromName("noSuchWorkflow") is None

    def test_the_documents_data_is_handed_to_getmesh(self, of, workflowInDB, monkeypatch):
        getMesh = _Recorder(returnValue="theMesh")
        monkeypatch.setattr(OFToolkit, "getMesh", lambda self, *a, **kw: getMesh(*a, **kw))
        assert of.getMeshFromName("flow_0000") == "theMesh"
        assert getMesh.lastArgs == (workflowInDB.getData(),)

    @pytest.mark.xfail(
        strict=True,
        reason="B193: getMeshFromName takes readParallel and time and "
               "documents both, but calls self.getMesh(doc.getData()) "
               "positionally with neither, so both are silently ignored and "
               "the mesh is always read from time 0 as a parallel case. "
               "getMeshExtentFromName has the same signature and the same "
               "hole (its callee getMeshExtent takes no such arguments at "
               "all). See the consolidated findings issue.",
    )
    def test_read_parallel_and_time_should_reach_getmesh(self, of, workflowInDB, monkeypatch):
        getMesh = _Recorder()
        monkeypatch.setattr(OFToolkit, "getMesh", lambda self, *a, **kw: getMesh(*a, **kw))
        of.getMeshFromName("flow_0000", readParallel=False, time=300)
        assert getMesh.lastKwargs == dict(readParallel=False, time=300)

    def test_read_parallel_and_time_are_currently_dropped(self, of, workflowInDB, monkeypatch):
        """Characterisation of B193."""
        getMesh = _Recorder()
        monkeypatch.setattr(OFToolkit, "getMesh", lambda self, *a, **kw: getMesh(*a, **kw))
        of.getMeshFromName("flow_0000", readParallel=False, time=300)
        assert getMesh.lastKwargs == {}


@pytest.mark.unit
class TestGetMeshExtentFromName:
    def test_an_unknown_name_gives_none(self, of):
        assert of.getMeshExtentFromName("noSuchWorkflow") is None

    def test_the_documents_data_is_handed_to_getmeshextent(self, of, workflowInDB, monkeypatch):
        getMeshExtent = _Recorder(returnValue="theExtent")
        monkeypatch.setattr(OFToolkit, "getMeshExtent", lambda self, *a, **kw: getMeshExtent(*a, **kw))
        assert of.getMeshExtentFromName("flow_0000") == "theExtent"
        assert getMeshExtent.lastArgs == (workflowInDB.getData(),)

    def test_read_parallel_and_time_are_currently_dropped(self, of, workflowInDB, monkeypatch):
        """Characterisation of B193, second site."""
        getMeshExtent = _Recorder()
        monkeypatch.setattr(OFToolkit, "getMeshExtent", lambda self, *a, **kw: getMeshExtent(*a, **kw))
        of.getMeshExtentFromName("flow_0000", readParallel=False, time=300)
        assert getMeshExtent.lastKwargs == {}


# ---------------------------------------------------------------------------
# createEmptyCase / writeEmptyField
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestCreateEmptyCase:
    def test_it_creates_the_standard_case_skeleton(self, of, tmp_path):
        case = tmp_path / "case"
        of.createEmptyCase(caseDirectory=str(case), fieldList=[],
                           flowType=OFToolkit.FLOWTYPE_INCOMPRESSIBLE)
        for expected in ["constant", "system", "0", "0.orig", "0.parallel",
                         os.path.join("constant", "triSurface")]:
            assert (case / expected).is_dir(), expected

    def test_every_field_is_written_three_times(self, of, tmp_path):
        case = tmp_path / "case"
        of.createEmptyCase(caseDirectory=str(case), fieldList=["U", "p"],
                           flowType=OFToolkit.FLOWTYPE_INCOMPRESSIBLE)
        for timeOrLocation in ["0", "0.orig", "0.parallel"]:
            assert sorted(os.listdir(case / timeOrLocation)) == ["U", "p"]

    def test_it_can_be_called_on_an_existing_case(self, of, tmp_path):
        case = tmp_path / "case"
        case.mkdir()
        (case / "constant").mkdir()
        of.createEmptyCase(caseDirectory=str(case), fieldList=["T"],
                           flowType=OFToolkit.FLOWTYPE_INCOMPRESSIBLE)
        assert (case / "0" / "T").exists()

    def test_a_file_in_the_way_is_reported_rather_than_overwritten(self, of, tmp_path):
        notADirectory = tmp_path / "case"
        notADirectory.write_text("I am a file")
        with pytest.raises(ValueError, match="exists as a file"):
            of.createEmptyCase(caseDirectory=str(notADirectory), fieldList=["U"],
                               flowType=OFToolkit.FLOWTYPE_INCOMPRESSIBLE)
        assert notADirectory.read_text() == "I am a file"

    def test_an_additional_field_definition_makes_the_field_writable(self, of, tmp_path):
        case = tmp_path / "case"
        of.createEmptyCase(caseDirectory=str(case), fieldList=["myTracer"],
                           flowType=OFToolkit.FLOWTYPE_INCOMPRESSIBLE,
                           additionalFieldsDescription=dict(
                               myTracer=dict(dimensions=dict(kg=1, m=-3), fieldType="scalar")))
        assert "myTracer" in of.OFObjectHome.fieldDefinitions
        assert (case / "0" / "myTracer").exists()

    def test_an_additional_field_definition_may_come_from_a_json_file(self, of, tmp_path):
        definitions = tmp_path / "fields.json"
        definitions.write_text(json.dumps(
            dict(myTracer=dict(dimensions=dict(kg=1, m=-3), fieldType="scalar"))))
        case = tmp_path / "case"
        of.createEmptyCase(caseDirectory=str(case), fieldList=["myTracer"],
                           flowType=OFToolkit.FLOWTYPE_INCOMPRESSIBLE,
                           additionalFieldsDescription=str(definitions))
        assert (case / "0.orig" / "myTracer").exists()

    def test_no_additional_field_definition_is_accepted(self, of, tmp_path):
        case = tmp_path / "case"
        of.createEmptyCase(caseDirectory=str(case), fieldList=["U"],
                           flowType=OFToolkit.FLOWTYPE_INCOMPRESSIBLE,
                           additionalFieldsDescription=None)
        assert (case / "0" / "U").exists()

    def test_an_unknown_field_is_reported(self, of, tmp_path):
        with pytest.raises(ValueError, match="noSuchField"):
            of.createEmptyCase(caseDirectory=str(tmp_path / "case"), fieldList=["noSuchField"],
                               flowType=OFToolkit.FLOWTYPE_INCOMPRESSIBLE)


@pytest.mark.unit
class TestWriteEmptyField:
    def test_it_writes_the_field_file_at_the_requested_location(self, of, tmp_path):
        of.writeEmptyField(fieldName="U", flowType=OFToolkit.FLOWTYPE_INCOMPRESSIBLE,
                           caseDirectory=str(tmp_path), timeOrLocation="0.orig")
        assert (tmp_path / "0.orig" / "U").exists()

    def test_the_file_is_named_after_the_field_file_not_the_field(self, of, tmp_path):
        # cellCenters is stored in a file called C.
        of.writeEmptyField(fieldName="cellCenters", flowType=OFToolkit.FLOWTYPE_INCOMPRESSIBLE,
                           caseDirectory=str(tmp_path), timeOrLocation=0)
        assert sorted(os.listdir(tmp_path / "0")) == ["C"]

    def test_an_unknown_field_is_reported(self, of, tmp_path):
        with pytest.raises(ValueError, match="noSuchField"):
            of.writeEmptyField(fieldName="noSuchField", flowType=OFToolkit.FLOWTYPE_INCOMPRESSIBLE,
                               caseDirectory=str(tmp_path))

    @pytest.fixture()
    def fieldSeam(self, monkeypatch):
        """Replace the field object, so the calls made on it are visible."""
        from unittest.mock import MagicMock

        field = MagicMock()
        getEmptyField = _Recorder(returnValue=field)
        monkeypatch.setattr(OFObjectHome, "getEmptyField", lambda self, **kw: getEmptyField(**kw))
        return types.SimpleNamespace(field=field, getEmptyField=getEmptyField)

    def test_the_field_name_and_flow_type_reach_the_object_home(self, of, tmp_path, fieldSeam):
        of.writeEmptyField(fieldName="U", flowType=OFToolkit.FLOWTYPE_COMPRESSIBLE,
                           caseDirectory=str(tmp_path))
        assert fieldSeam.getEmptyField.lastKwargs == dict(fieldName="U", flowType="compressible")

    def test_the_case_directory_and_location_reach_the_field(self, of, tmp_path, fieldSeam):
        of.writeEmptyField(fieldName="U", flowType=OFToolkit.FLOWTYPE_INCOMPRESSIBLE,
                           caseDirectory=str(tmp_path), timeOrLocation="0.parallel")
        fieldSeam.field.writeToCase.assert_called_once_with(caseDirectory=str(tmp_path),
                                                            timeOrLocation="0.parallel")

    def test_the_boundaries_are_read_from_the_case_only_when_asked(self, of, tmp_path, fieldSeam):
        of.writeEmptyField(fieldName="U", flowType=OFToolkit.FLOWTYPE_INCOMPRESSIBLE,
                           caseDirectory=str(tmp_path))
        fieldSeam.field.readBoundariesFromCase.assert_not_called()
        of.writeEmptyField(fieldName="U", flowType=OFToolkit.FLOWTYPE_INCOMPRESSIBLE,
                           caseDirectory=str(tmp_path), readBoundaryFromCase=True)
        fieldSeam.field.readBoundariesFromCase.assert_called_once_with(str(tmp_path), readParallel=True)

    def test_the_processor_boundary_is_added_only_when_asked(self, of, tmp_path, fieldSeam):
        of.writeEmptyField(fieldName="U", flowType=OFToolkit.FLOWTYPE_INCOMPRESSIBLE,
                           caseDirectory=str(tmp_path))
        fieldSeam.field.addProcBoundary.assert_not_called()
        of.writeEmptyField(fieldName="U", flowType=OFToolkit.FLOWTYPE_INCOMPRESSIBLE,
                           caseDirectory=str(tmp_path), writeProcBoundary=True)
        fieldSeam.field.addProcBoundary.assert_called_once_with()


# ---------------------------------------------------------------------------
# The datasetToOF delegations
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestDatasetToOFDelegations:
    """Only the forwarding is pinned: the conversion itself needs foamlib
    and is covered (where installed) by test_openfoam_dataset2of.py."""

    @pytest.fixture()
    def dataset(self):
        import numpy
        import xarray

        return xarray.Dataset(dict(T=(("x", "y"), numpy.zeros((2, 2)))),
                              coords=dict(x=[0.0, 1.0], y=[0.0, 1.0]))

    def test_case_geometry_is_delegated(self, of, monkeypatch):
        recorder = _Recorder(returnValue="theGeometry")
        monkeypatch.setattr(datasetToOF, "caseGeometry", recorder)
        assert of.caseGeometry("/case", time=5, patchList=["ground"], readParallel=False) == "theGeometry"
        assert recorder.lastKwargs == dict(caseDirectory="/case", time=5,
                                           patchList=["ground"], readParallel=False)

    def test_dataset_to_case_fields_is_delegated(self, of, monkeypatch, dataset):
        recorder = _Recorder(returnValue=["/case/0/T"])
        monkeypatch.setattr(datasetToOF, "datasetToCaseFields", recorder)
        assert of.datasetToCaseFields("/case", dataset, dict(T="T"), "x", "y",
                                      writeTime=10) == ["/case/0/T"]
        assert recorder.lastKwargs["caseDirectory"] == "/case"
        assert recorder.lastKwargs["fieldMap"] == dict(T="T")
        assert recorder.lastKwargs["xCoordinate"] == "x"
        assert recorder.lastKwargs["yCoordinate"] == "y"
        assert recorder.lastKwargs["writeTime"] == 10
        assert recorder.lastKwargs["dataset"] is dataset

    def test_netcdf_to_case_fields_opens_the_file_and_delegates(self, of, monkeypatch, dataset):
        import xarray

        openDataset = _Recorder(returnValue=dataset)
        monkeypatch.setattr(xarray, "open_dataset", openDataset)
        delegate = _Recorder(returnValue=["/case/0/T"])
        monkeypatch.setattr(OFToolkit, "datasetToCaseFields", lambda self, **kw: delegate(**kw))
        assert of.netcdfToCaseFields("/case", "/data/meteo.nc", dict(T="T"), "x", "y") == ["/case/0/T"]
        assert openDataset.lastArgs == ("/data/meteo.nc",)
        assert delegate.lastKwargs["dataset"] is dataset

    def test_dataset_to_set_fields_dict_is_delegated(self, of, monkeypatch, dataset):
        recorder = _Recorder(returnValue=["regions"])
        monkeypatch.setattr(datasetToOF, "datasetToSetFieldsDict", recorder)
        assert of.datasetToSetFieldsDict(dataset, dict(T="T"), "x", "y", "z",
                                         includeFaces=False) == ["regions"]
        assert recorder.lastKwargs["verticalCoordinate"] == "z"
        assert recorder.lastKwargs["includeFaces"] is False

    def test_netcdf_to_set_fields_dict_opens_the_file_and_delegates(self, of, monkeypatch, dataset):
        import xarray

        openDataset = _Recorder(returnValue=dataset)
        monkeypatch.setattr(xarray, "open_dataset", openDataset)
        delegate = _Recorder(returnValue=["regions"])
        monkeypatch.setattr(OFToolkit, "datasetToSetFieldsDict", lambda self, **kw: delegate(**kw))
        assert of.netcdfToSetFieldsDict("/data/meteo.nc", dict(T="T"), "x", "y", "z") == ["regions"]
        assert openDataset.lastArgs == ("/data/meteo.nc",)
        assert delegate.lastKwargs["verticalCoordinate"] == "z"


# ---------------------------------------------------------------------------
# getTimeList: the decomposed branch (the single-processor one is B83,
# pinned in test_openfoam_toolkit_pure.py)
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetTimeListDecomposed:
    @pytest.fixture()
    def decomposedCase(self, tmp_path):
        case = tmp_path / "myCase"
        for processorName in ["processor0", "processor1"]:
            for timeName in ["0", "100", "50"]:
                (case / processorName / timeName).mkdir(parents=True)
            (case / processorName / "constant").mkdir()
        return case

    def test_it_lists_the_time_directories_of_the_first_processor(self, of, decomposedCase):
        assert of.getTimeList(str(decomposedCase)) == [0.0, 50.0, 100.0]

    def test_return_first_false_keys_the_result_by_case_name(self, of, decomposedCase):
        assert of.getTimeList(str(decomposedCase), returnFirst=False) == {"myCase": [0.0, 50.0, 100.0]}

    def test_a_workflow_in_the_db_is_looked_up_by_name(self, of, workflowInDB, monkeypatch, tmp_path):
        # The document's data is a file, not a case directory, so the lookup
        # is what is being pinned here: the name resolves to the document.
        case = tmp_path / "fromDB"
        (case / "processor0" / "0").mkdir(parents=True)
        monkeypatch.setattr(type(workflowInDB), "getData", lambda self, **kw: str(case))
        assert of.getTimeList("flow_0000") == [0.0]


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestAnalysis:
    def test_the_datalayer_is_the_toolkit(self, of):
        assert of.analysis.datalayer is of

    def test_the_vtk_pipeline_is_bound_to_the_toolkit(self, of):
        from hera.simulations.openFoam.postProcess.VTKPipeline import VTKPipeLine

        pipeline = of.analysis.getVTKPipeline()
        assert isinstance(pipeline, VTKPipeLine)

    def test_a_new_pipeline_is_returned_every_time(self, of):
        assert of.analysis.getVTKPipeline() is not of.analysis.getVTKPipeline()

    def test_execute_and_load_forwards_both_halves(self, of):
        from unittest.mock import MagicMock

        pipeline = MagicMock()
        of.analysis.executeVTKFiltersAndLoad(pipeline, sourceOrName="reader", timeList=[1, 2],
                                             tsBlockNum=7, overwrite=True, append=True,
                                             overwriteMetadata=True)
        pipeline.execute.assert_called_once_with(sourceOrName="reader", timeList=[1, 2],
                                                 tsBlockNum=7, overwrite=True, append=True)
        pipeline.loadToProject.assert_called_once_with(datalayer=of, overwrite=True)

    def test_execute_and_load_defaults_are_conservative(self, of):
        from unittest.mock import MagicMock

        pipeline = MagicMock()
        of.analysis.executeVTKFiltersAndLoad(pipeline)
        assert pipeline.execute.call_args.kwargs == dict(sourceOrName=None, timeList=None,
                                                         tsBlockNum=50, overwrite=False, append=False)
        assert pipeline.loadToProject.call_args.kwargs["overwrite"] is False

    def test_it_finds_the_filters_of_a_simulation(self, of, workflowInDB, monkeypatch):
        # getFiltersDocuments reads desc['simulationName'], which
        # addWorkflowToGroup does not write, so the document is built here.
        of.addCacheDocument(resource="/x", dataFormat="string", type="vtk_filter",
                            desc=dict(simulationName="flow_0000", filterName="f1"))
        of.addCacheDocument(resource="/y", dataFormat="string", type="vtk_filter",
                            desc=dict(simulationName="flow_0000", filterName="f2"))
        monkeypatch.setattr(OFToolkit, "getWorkflowDocumentFromDB",
                            lambda self, *a, **kw: dict(desc=dict(simulationName="flow_0000")))
        assert len(of.analysis.getFiltersDocuments("flow_0000")) == 2
        assert len(of.analysis.getFiltersDocuments("flow_0000", filterName="f1")) == 1

    @pytest.mark.xfail(
        strict=True,
        reason="B196: getFiltersDocuments guards with `if wrkflow is "
               "None`, but getWorkflowDocumentFromDB returns a LIST (empty "
               "when nothing matches), which is never None. The intended "
               "ValueError naming the missing case is therefore unreachable "
               "and an unknown simulation dies with TypeError: list indices "
               "must be integers. See the consolidated findings issue.",
    )
    def test_an_unknown_simulation_should_be_reported_as_a_valueerror(self, of):
        with pytest.raises(ValueError, match="was not found"):
            of.analysis.getFiltersDocuments("noSuchSimulation")

    def test_an_unknown_simulation_currently_raises_typeerror(self, of):
        """Characterisation of B196."""
        with pytest.raises(TypeError, match="list indices"):
            of.analysis.getFiltersDocuments("noSuchSimulation")


# ---------------------------------------------------------------------------
# Presentation
# ---------------------------------------------------------------------------

def _particleData():
    return pandas.DataFrame(dict(time=[1.0, 1.0, 2.0],
                                 globalX=[0.0, 1.0, 2.0],
                                 globalY=[0.0, 1.0, 2.0],
                                 globalZ=[0.0, 1.0, 2.0],
                                 mass=[1.0, 2.0, 3.0]))


@pytest.mark.unit
class TestPresentationToParaviewCSV:
    @pytest.fixture()
    def outputDirectory(self, tmp_path):
        # to_paraview_CSV writes into an existing directory (it does not
        # create one), and tmp_path itself already holds heraFiles.
        directory = tmp_path / "csv"
        directory.mkdir()
        return directory

    def test_one_file_per_time_step(self, of, outputDirectory):
        of.presentation.to_paraview_CSV(_particleData(), str(outputDirectory), "particles")
        assert sorted(os.listdir(outputDirectory)) == ["particles_1.csv", "particles_2.csv"]

    def test_only_the_global_positions_are_written(self, of, outputDirectory):
        of.presentation.to_paraview_CSV(_particleData(), str(outputDirectory), "particles")
        content = (outputDirectory / "particles_1.csv").read_text().split()
        assert content[0] == "globalX,globalY,globalZ"
        assert content[1:] == ["0.0,0.0,0.0", "1.0,1.0,1.0"]

    def test_the_time_factor_scales_the_file_name(self, of, outputDirectory):
        of.presentation.to_paraview_CSV(_particleData(), str(outputDirectory), "particles",
                                        timeFactor=100)
        assert sorted(os.listdir(outputDirectory)) == ["particles_100.csv", "particles_200.csv"]


@pytest.mark.unit
class TestPresentationToUnstructuredVTK:
    @pytest.fixture()
    def vtkSeam(self, monkeypatch):
        recorder = _Recorder()
        monkeypatch.setattr(ofToolkitModule, "evtk_hl",
                            types.SimpleNamespace(pointsToVTK=recorder, structuredToVTK=_Recorder()))
        return recorder

    def test_one_dataset_is_written_per_time_step(self, of, tmp_path, vtkSeam):
        of.presentation.toUnstructuredVTK(_particleData(), str(tmp_path / "out"), "particles")
        written = [call[0][0] for call in vtkSeam.calls]
        assert written == [str(tmp_path / "out" / "particles_0"), str(tmp_path / "out" / "particles_1")]

    def test_the_output_directory_is_created(self, of, tmp_path, vtkSeam):
        of.presentation.toUnstructuredVTK(_particleData(), str(tmp_path / "out"), "particles")
        assert (tmp_path / "out").is_dir()

    def test_the_coordinates_are_split_out_of_the_field_list(self, of, tmp_path, vtkSeam):
        of.presentation.toUnstructuredVTK(_particleData(), str(tmp_path), "particles")
        path, x, y, z, fields = vtkSeam.lastArgs
        assert sorted(fields) == ["mass", "time"]
        assert list(x) == [2.0] and list(y) == [2.0] and list(z) == [2.0]

    def test_an_explicit_field_list_is_honoured(self, of, tmp_path, vtkSeam):
        of.presentation.toUnstructuredVTK(_particleData(), str(tmp_path), "particles",
                                          fieldList=["mass"])
        assert list(vtkSeam.lastArgs[4]) == ["mass"]

    def test_the_coordinate_columns_are_configurable(self, of, tmp_path, vtkSeam):
        data = _particleData().rename(columns=dict(globalX="x", globalY="y", globalZ="z"))
        of.presentation.toUnstructuredVTK(data, str(tmp_path), "particles",
                                          xcoord="x", ycoord="y", zcoord="z")
        assert sorted(vtkSeam.lastArgs[4]) == ["mass", "time"]

    @pytest.mark.xfail(
        strict=True,
        reason="B197: the overwrite guard tests os.path.exists on the "
               "extension-less base path it hands to evtk, but evtk writes "
               "<base>.vtu. The file that exists is therefore never the file "
               "checked, so the FileExistsError that is supposed to protect "
               "an existing output (and the os.remove that overwrite=True "
               "does) can never fire for a real evtk output. "
               "See the consolidated findings issue.",
    )
    def test_an_existing_output_should_be_refused_without_overwrite(self, of, tmp_path, vtkSeam):
        (tmp_path / "particles_0.vtu").write_text("previous run")
        with pytest.raises(FileExistsError):
            of.presentation.toUnstructuredVTK(_particleData(), str(tmp_path), "particles")

    def test_an_existing_vtu_is_currently_ignored(self, of, tmp_path, vtkSeam):
        """Characterisation of B197."""
        (tmp_path / "particles_0.vtu").write_text("previous run")
        of.presentation.toUnstructuredVTK(_particleData(), str(tmp_path), "particles")
        assert len(vtkSeam.calls) == 2

    def test_the_guard_fires_for_the_extensionless_path_it_checks(self, of, tmp_path, vtkSeam):
        """Characterisation of B197: the guard itself works, on the wrong path."""
        (tmp_path / "particles_0").write_text("previous run")
        with pytest.raises(FileExistsError, match="overwrite"):
            of.presentation.toUnstructuredVTK(_particleData(), str(tmp_path), "particles")
        of.presentation.toUnstructuredVTK(_particleData(), str(tmp_path), "particles", overwrite=True)
        assert not (tmp_path / "particles_0").exists()


@pytest.mark.unit
class TestPresentationToStructuredVTK:
    _extents = dict(xmin=0, xmax=3, ymin=0, ymax=3, zmin=0, zmax=3)

    @pytest.mark.xfail(
        strict=True,
        reason="B198: Presentation.toStructuredVTK calls "
               "self.analysis.calcConcentrationTimeStepFullMesh, but the "
               "Analysis class of hera/simulations/openFoam/toolkit.py has "
               "no such method -- it lives on the *LSM* Analysis "
               "(hera/simulations/openFoam/lagrangian/LSM/toolkit.py). The "
               "method therefore raises AttributeError for every input. "
               "See the consolidated findings issue.",
    )
    def test_it_should_write_a_structured_grid_per_time_step(self, of, tmp_path):
        of.presentation.toStructuredVTK(_particleData(), str(tmp_path), "concentration",
                                        extents=self._extents, dxdydz=1)

    def test_it_currently_raises_attributeerror(self, of, tmp_path):
        """Characterisation of B198."""
        with pytest.raises(AttributeError, match="calcConcentrationTimeStepFullMesh"):
            of.presentation.toStructuredVTK(_particleData(), str(tmp_path), "concentration",
                                            extents=self._extents, dxdydz=1)


@pytest.mark.unit
class TestPresentationLoadLagrangianDataParallel:
    @pytest.fixture()
    def readSeam(self, monkeypatch):
        recorder = _Recorder(returnValue=pandas.DataFrame(dict(time=[1.0], globalX=[0.0])))
        monkeypatch.setattr(ofToolkitModule, "readLagrangianRecord", recorder)
        return recorder

    @pytest.fixture()
    def parallelCase(self, tmp_path):
        for timeName in ["0", "100", "200"]:
            (tmp_path / "processor0" / timeName).mkdir(parents=True)
        return tmp_path

    def test_a_case_without_processor_directories_is_reported(self, of, tmp_path, readSeam):
        with pytest.raises(ValueError, match="processor"):
            of.presentation.loadLagrangianDataParallel(str(tmp_path))

    def test_the_reader_options_are_forwarded_for_every_record(self, of, parallelCase, readSeam):
        of.presentation.loadLagrangianDataParallel(str(parallelCase), withVelocity=True,
                                                   withMass=True, withReleaseTimes=True,
                                                   cloudName="myCloud").compute()
        assert readSeam.lastKwargs == dict(casePath=str(parallelCase), withVelocity=True,
                                           withReleaseTimes=True, withMass=True, cloudName="myCloud")

    def test_every_processor_and_time_directory_is_read(self, of, parallelCase, readSeam):
        of.presentation.loadLagrangianDataParallel(str(parallelCase)).compute()
        assert {call[0][0] for call in readSeam.calls} == {os.path.join("processor0", "100"),
                                                           os.path.join("processor0", "200")}

    @pytest.mark.xfail(
        strict=True,
        reason="B199: the record list is built from timeList[1:], so the "
               "first time step is always dropped -- including when the "
               "caller passed `times` explicitly, which the docstring "
               "describes as 'a list of time steps to extract'. A single "
               "requested time therefore yields an empty delayed list and "
               "dask raises 'Must supply at least one delayed object'. "
               "See the consolidated findings issue.",
    )
    def test_an_explicitly_requested_time_should_be_read(self, of, parallelCase, readSeam):
        of.presentation.loadLagrangianDataParallel(str(parallelCase), times=["100"]).compute()
        assert [call[0][0] for call in readSeam.calls] == [os.path.join("processor0", "100")]

    def test_a_single_requested_time_currently_leaves_nothing_to_read(self, of, parallelCase, readSeam):
        """Characterisation of B199."""
        with pytest.raises(TypeError, match="at least one delayed object"):
            of.presentation.loadLagrangianDataParallel(str(parallelCase), times=["100"])

    def test_two_requested_times_currently_read_only_the_second(self, of, parallelCase, readSeam):
        """Characterisation of B199."""
        of.presentation.loadLagrangianDataParallel(str(parallelCase), times=["100", "200"]).compute()
        assert {call[0][0] for call in readSeam.calls} == {os.path.join("processor0", "200")}

    @pytest.mark.xfail(
        strict=True,
        reason="B200: `times` is documented as 'a list of time steps to "
               "extract' and is pushed through numpy.atleast_1d, which turns "
               "numbers into numpy.int64 -- but each element is then used as "
               "a path component in os.path.join, which only accepts str / "
               "bytes / PathLike. Numeric times (the natural form for a "
               "time step) therefore always raise TypeError; only strings "
               "work. See the consolidated findings issue.",
    )
    def test_numeric_times_should_be_accepted(self, of, parallelCase, readSeam):
        of.presentation.loadLagrangianDataParallel(str(parallelCase), times=[100, 200]).compute()
        assert readSeam.calls != []

    def test_numeric_times_currently_raise(self, of, parallelCase, readSeam):
        """Characterisation of B200."""
        with pytest.raises(TypeError, match="join"):
            of.presentation.loadLagrangianDataParallel(str(parallelCase), times=[100, 200])
