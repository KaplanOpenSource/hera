"""hera/simulations/CLI.py: the `hera-workflows` command handlers.

Every handler here takes the argparse Namespace produced by
`hera/bin/hera-workflows` and drives the hermesWorkflows toolkit.  These
tests build the Namespace by hand (mirroring the dests declared in
hera/bin/hera-workflows) and substitute a recording stand-in for the
toolkit, so what is under test is the CLI plumbing only: which toolkit
call each subcommand makes, which arguments it forwards, which branch it
takes on the flags, what it prints, and what it refuses to do on bad
input.  Nothing touches MongoDB, Luigi, hermes or OpenFOAM.

Deliberately not tested:
  * the toolkit methods themselves (listGroups, compareWorkflows,
    executeWorkflowFromDB, ...) -- they are the workflow engine, covered
    elsewhere, and are the seam these tests stub out;
  * `hermes.utils.workflowAssembly.handler_build` / `handler_execute`,
    which workflow_build/workflow_execute delegate to unchanged: hermes
    is stubbed away in CI (see hera/tests/unit/_stubs.py), so only the
    delegation is asserted, against a fake hermes module;
  * the argparse wiring in hera/bin/hera-workflows, which is a script and
    not importable as a module.

Bugs pinned here (each with an xfail for the intended behaviour plus a
passing characterisation of what actually happens):
  * workflow_delete tells the user to run a script whose name it never
    wrote ("completeRemove.py" vs the "completeDelete.py" it writes);
  * workflow_list reads `arguments.object`, an attribute the `list
    workflows` parser never defines, on the very path that reports "not
    found" -- and then indexes the empty list anyway;
  * create_workflow_variations logs an error for a missing variation
    JSON / an unknown workflow and then carries on into
    FileNotFoundError / IndexError;
  * workflow_compare with the default `--format pandas` plus `--file`
    tries to write a DataFrame object to a text file (TypeError).
"""
import json
import sys
import types
from argparse import Namespace

import pandas
import pytest

from hera.simulations.CLI import (
    WorkflowsGroup_list,
    batch_delete_workflows,
    create_workflow_variations,
    sorround_with_char,
    workflow_add,
    workflow_build,
    workflow_buildExecute,
    workflow_compare,
    workflow_compareToDisk,
    workflow_delete,
    workflow_execute,
    workflow_export,
    workflow_list,
    workflow_sync_to_db,
    workflowNodes_list,
    workflowNodes_listParameters,
    _confirm_project_name,
)

PROJECT = "UNIT_WORKFLOW_PROJECT"


# ---------------------------------------------------------------------------
# stand-ins
# ---------------------------------------------------------------------------

class _FakeNode:
    def __init__(self, parameters):
        self.parameters = parameters


class _FakeHermesWorkflow:
    """Stand-in for a hermes workflow object."""

    def __init__(self, name="wf", resource="", parametersJSON=None, nodes=None):
        self.name = name
        self.resource = resource
        self.Resource_path = resource
        self.parametersJSON = parametersJSON if parametersJSON is not None else {}
        self._nodes = nodes if nodes is not None else {}

    @property
    def nodeList(self):
        return list(self._nodes)

    def items(self):
        return [(name, _FakeNode(params)) for name, params in self._nodes.items()]

    def __getitem__(self, name):
        return _FakeNode(self._nodes[name])


class _FakeDocument:
    """Stand-in for a datalayer document of a stored workflow."""

    def __init__(self, workflowName="wf_0", workflow=None, groupName="wf",
                 resource="/data/wf_0"):
        self.desc = {
            "workflowName": workflowName,
            "groupName": groupName,
            "workflow": workflow if workflow is not None else {"nodes": {}},
        }
        self.resource = resource
        self.deleted = False

    def __getitem__(self, key):
        if key == "resource":
            return self.resource
        if key == "desc":
            return self.desc
        raise KeyError(key)

    def to_mongo(self):
        return {"desc": self.desc, "resource": self.resource}

    def delete(self):
        self.deleted = True


class _FakeWorkflowToolkit:
    """Records every call the CLI makes, and returns canned answers."""

    DESC_GROUPNAME = "groupName"

    def __init__(self, documents=(), hermesWorkflow=None, localWorkflow=None,
                 comparison=None, workflows=()):
        self.documents = list(documents)
        self.hermesWorkflow = hermesWorkflow
        self.localWorkflow = localWorkflow
        self.comparison = comparison
        self.workflows = list(workflows)
        self.calls = []

    # -- helpers used by the assertions ------------------------------------
    def _record(self, _callName, *args, **kwargs):
        # The leading underscore keeps this from colliding with a recorded
        # keyword literally called "name".
        self.calls.append((_callName, args, kwargs))

    def names(self):
        return [call[0] for call in self.calls]

    def call(self, name):
        matching = [call for call in self.calls if call[0] == name]
        assert len(matching) == 1, f"{name} called {len(matching)} times: {self.calls}"
        return matching[0]

    # -- the toolkit surface the CLI uses ----------------------------------
    def listGroups(self, solver=None, workflowName=True):
        self._record("listGroups", solver=solver, workflowName=workflowName)

    def addWorkflowFileInGroup(self, workflowFilePath, *args, **kwargs):
        self._record("addWorkflowFileInGroup", workflowFilePath, *args, **kwargs)
        return workflowFilePath

    def getWorkflowListDocumentFromDB(self, nameOrWorkflow, **query):
        self._record("getWorkflowListDocumentFromDB", nameOrWorkflow, **query)
        return self.documents

    def getHermesWorkflowFromDB(self, nameOrWorkflow, returnFirst=True, **query):
        self._record("getHermesWorkflowFromDB", nameOrWorkflow,
                     returnFirst=returnFirst, **query)
        if returnFirst:
            return self.hermesWorkflow
        return [self.hermesWorkflow] if self.hermesWorkflow is not None else []

    def getHemresWorkflowFromDocument(self, documentList, returnFirst=True):
        self._record("getHemresWorkflowFromDocument", documentList=documentList,
                     returnFirst=returnFirst)
        return self.hermesWorkflow

    def getHermesWorkflowFromJSON(self, workflow, name=None, resource=None):
        self._record("getHermesWorkflowFromJSON", workflow, name=name, resource=resource)
        return self.localWorkflow

    def updateDocumentWorkflow(self, document, workflow):
        self._record("updateDocumentWorkflow", document=document, workflow=workflow)

    def listWorkflows(self, workflowGroup, listNodes=False, listParameters=False):
        self._record("listWorkflows", workflowGroup=workflowGroup,
                     listNodes=listNodes, listParameters=listParameters)
        return self.workflows

    def compareWorkflows(self, Workflow, longFormat=False, transpose=False):
        self._record("compareWorkflows", Workflow, longFormat=longFormat,
                     transpose=transpose)
        return self.comparison

    def deleteWorkflowInGroup(self, workflowGroup, deepDelete=False,
                              resetCounter=True, exclude=()):
        self._record("deleteWorkflowInGroup", workflowGroup, deepDelete,
                     resetCounter=resetCounter, exclude=exclude)

    def executeWorkflowFromDB(self, nameOrWorkflow, **kwargs):
        self._record("executeWorkflowFromDB", nameOrWorkflow, **kwargs)


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def install_toolkit(monkeypatch):
    """Make `toolkitHome.getToolkit` hand the CLI our stand-in.

    Returns a callable that installs the toolkit and yields a dict which is
    filled in with the arguments getToolkit was asked for, so the tests can
    assert that the right toolkit of the right project was requested.
    """
    from hera import toolkitHome
    import hera.simulations.hermesWorkflowToolkit as hermesModule

    def _install(toolkit):
        asked = {}

        def _getToolkit(self, toolkitName=None, filesDirectory=None, **kwargs):
            asked.update(toolkitName=toolkitName, kwargs=kwargs)
            return toolkit

        monkeypatch.setattr(type(toolkitHome), "getToolkit", _getToolkit)
        # Three handlers assert isinstance(wftk, hermesWorkflowToolkit).
        monkeypatch.setattr(hermesModule, "hermesWorkflowToolkit", type(toolkit))
        return asked

    return _install


@pytest.fixture()
def in_tmp_cwd(tmp_path, monkeypatch):
    """Several handlers write files relative to the process cwd."""
    monkeypatch.chdir(tmp_path)
    return tmp_path


@pytest.fixture()
def fake_hermes(monkeypatch):
    """Replace the whole `hermes` package with a recording fake.

    hermes is stubbed out in CI (_stubs.py) and a real local install must
    not be executed either, so the module chain is replaced wholesale.
    """
    recorded = {}

    def _workflow(jsonData, Resource_path=None):
        recorded.update(json=jsonData, Resource_path=Resource_path)
        return recorded["result"]

    def _handler_build(arguments):
        recorded["build"] = arguments

    def _handler_execute(arguments):
        recorded["execute"] = arguments

    hermes = types.ModuleType("hermes")
    hermes.__path__ = []
    hermes.workflow = _workflow
    utils = types.ModuleType("hermes.utils")
    utils.__path__ = []
    assembly = types.ModuleType("hermes.utils.workflowAssembly")
    assembly.handler_build = _handler_build
    assembly.handler_execute = _handler_execute

    monkeypatch.setitem(sys.modules, "hermes", hermes)
    monkeypatch.setitem(sys.modules, "hermes.utils", utils)
    monkeypatch.setitem(sys.modules, "hermes.utils.workflowAssembly", assembly)
    return recorded


# ---------------------------------------------------------------------------
# the pure helpers
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestSorroundWithChar:
    def test_the_text_is_padded_symmetrically_to_the_requested_width(self):
        assert sorround_with_char("ab", 8) == "---ab---"

    def test_the_padding_character_is_configurable(self):
        assert sorround_with_char("ab", 6, char="*") == "**ab**"

    def test_a_text_wider_than_the_requested_width_is_returned_unchanged(self):
        assert sorround_with_char("a very long title", 4) == "a very long title"

    def test_a_text_exactly_as_wide_as_requested_gets_no_padding(self):
        assert sorround_with_char("abcd", 4) == "abcd"

    def test_an_odd_gap_yields_a_line_one_character_short(self):
        """The padding is halved with //, so an odd gap loses a character."""
        assert sorround_with_char("abc", 10) == "---abc---"
        assert len(sorround_with_char("abc", 10)) == 9


@pytest.mark.unit
class TestConfirmProjectName:
    def test_an_explicit_project_name_is_left_untouched(self, in_tmp_cwd):
        arguments = Namespace(projectName=PROJECT)
        _confirm_project_name(arguments, _logger())
        assert arguments.projectName == PROJECT

    def test_a_missing_project_name_is_read_from_the_case_configuration(self, in_tmp_cwd):
        (in_tmp_cwd / "caseConfiguration.json").write_text(
            json.dumps({"projectName": "FROM_CASE_CONFIG"})
        )
        arguments = Namespace(projectName=None)
        _confirm_project_name(arguments, _logger())
        assert arguments.projectName == "FROM_CASE_CONFIG"

    def test_a_missing_case_configuration_file_raises_a_value_error(self, in_tmp_cwd):
        arguments = Namespace(projectName=None)
        with pytest.raises(ValueError):
            _confirm_project_name(arguments, _logger())

    def test_a_case_configuration_without_a_project_name_raises_a_key_error(self, in_tmp_cwd):
        (in_tmp_cwd / "caseConfiguration.json").write_text(json.dumps({"other": 1}))
        arguments = Namespace(projectName=None)
        with pytest.raises(KeyError):
            _confirm_project_name(arguments, _logger())


def _logger():
    import logging

    return logging.getLogger("hera.tests.unit.simulations_cli")


# ---------------------------------------------------------------------------
# list groups / add
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestWorkflowsGroupList:
    def test_the_solver_and_workflowname_flags_are_forwarded_to_the_toolkit(self, install_toolkit):
        toolkit = _FakeWorkflowToolkit()
        asked = install_toolkit(toolkit)
        WorkflowsGroup_list(
            Namespace(projectName=PROJECT, solver="simpleFoam", workflowName=True)
        )
        assert asked["kwargs"]["projectName"] == PROJECT
        assert asked["toolkitName"] == "hermesWorkflows"
        assert toolkit.call("listGroups")[2] == dict(solver="simpleFoam", workflowName=True)

    def test_a_missing_project_name_is_resolved_from_the_case_configuration(
        self, install_toolkit, in_tmp_cwd
    ):
        (in_tmp_cwd / "caseConfiguration.json").write_text(
            json.dumps({"projectName": "FROM_CASE_CONFIG"})
        )
        toolkit = _FakeWorkflowToolkit()
        asked = install_toolkit(toolkit)
        arguments = Namespace(projectName=None, solver=None, workflowName=False)
        WorkflowsGroup_list(arguments)
        assert arguments.projectName == "FROM_CASE_CONFIG"
        assert asked["kwargs"]["projectName"] == "FROM_CASE_CONFIG"


@pytest.mark.unit
class TestWorkflowAdd:
    def test_the_workflow_file_is_handed_to_the_toolkit(self, install_toolkit):
        toolkit = _FakeWorkflowToolkit()
        asked = install_toolkit(toolkit)
        workflow_add(Namespace(projectName=PROJECT, workflow="myFlow.json"))
        assert asked["kwargs"]["projectName"] == PROJECT
        assert toolkit.call("addWorkflowFileInGroup")[1] == ("myFlow.json",)

    def test_the_other_add_flags_are_accepted_and_then_ignored(self, install_toolkit):
        """Characterisation: the `add` parser declares --workflowGroup,
        --overwrite, --assignName, --allowDuplicate and --execute, but the
        handler forwards only the file path -- addWorkflowFileInGroup takes
        no equivalent argument, so every one of those flags is inert."""
        toolkit = _FakeWorkflowToolkit()
        install_toolkit(toolkit)
        workflow_add(
            Namespace(
                projectName=PROJECT,
                workflow="myFlow.json",
                workflowGroup="someGroup",
                overwrite=True,
                assignName=True,
                force=True,
                execute=True,
            )
        )
        name, args, kwargs = toolkit.call("addWorkflowFileInGroup")
        assert args == ("myFlow.json",)
        assert kwargs == {}


# ---------------------------------------------------------------------------
# delete
# ---------------------------------------------------------------------------

def _deleteArgs(**overrides):
    arguments = dict(projectName=PROJECT, workflows=["wf_0"], Export=False,
                     forceOverwrite=False)
    arguments.update(overrides)
    return Namespace(**arguments)


@pytest.mark.unit
class TestWorkflowDelete:
    def test_every_matching_document_is_deleted_from_the_database(
        self, install_toolkit, in_tmp_cwd
    ):
        documents = [_FakeDocument("wf_0"), _FakeDocument("wf_1")]
        toolkit = _FakeWorkflowToolkit(documents=documents)
        install_toolkit(toolkit)
        workflow_delete(_deleteArgs(workflows=["wf_0", "wf_1"]))
        assert [document.deleted for document in documents] == [True, True]
        assert toolkit.call("getWorkflowListDocumentFromDB")[1] == (["wf_0", "wf_1"],)

    def test_without_the_export_flag_no_json_is_written(self, install_toolkit, in_tmp_cwd):
        toolkit = _FakeWorkflowToolkit(documents=[_FakeDocument("wf_0")])
        install_toolkit(toolkit)
        workflow_delete(_deleteArgs())
        assert not (in_tmp_cwd / "wf_0.json").exists()

    def test_with_the_export_flag_the_workflow_is_saved_to_the_current_directory(
        self, install_toolkit, in_tmp_cwd
    ):
        document = _FakeDocument("wf_0", workflow={"nodes": {"A": 1}})
        toolkit = _FakeWorkflowToolkit(documents=[document])
        install_toolkit(toolkit)
        workflow_delete(_deleteArgs(Export=True))
        saved = json.loads((in_tmp_cwd / "wf_0.json").read_text())
        assert saved == {"workflow": {"nodes": {"A": 1}}}
        assert document.deleted is True

    def test_an_existing_export_file_cancels_the_deletion(self, install_toolkit, in_tmp_cwd):
        (in_tmp_cwd / "wf_0.json").write_text("do not touch me")
        document = _FakeDocument("wf_0")
        toolkit = _FakeWorkflowToolkit(documents=[document])
        install_toolkit(toolkit)
        workflow_delete(_deleteArgs(Export=True))
        assert (in_tmp_cwd / "wf_0.json").read_text() == "do not touch me"
        assert document.deleted is False

    def test_forceoverwrite_replaces_the_existing_export_and_deletes_anyway(
        self, install_toolkit, in_tmp_cwd
    ):
        (in_tmp_cwd / "wf_0.json").write_text("stale")
        document = _FakeDocument("wf_0", workflow={"nodes": {}})
        toolkit = _FakeWorkflowToolkit(documents=[document])
        install_toolkit(toolkit)
        workflow_delete(_deleteArgs(Export=True, forceOverwrite=True))
        assert json.loads((in_tmp_cwd / "wf_0.json").read_text()) == {"workflow": {"nodes": {}}}
        assert document.deleted is True

    def test_a_helper_script_with_one_rmtree_per_removed_resource_is_written(
        self, install_toolkit, in_tmp_cwd
    ):
        toolkit = _FakeWorkflowToolkit(
            documents=[_FakeDocument("wf_0", resource="/data/wf_0"),
                       _FakeDocument("wf_1", resource="/data/wf_1")]
        )
        install_toolkit(toolkit)
        workflow_delete(_deleteArgs(workflows=["wf_0", "wf_1"]))
        script = (in_tmp_cwd / "completeDelete.py").read_text()
        assert script.splitlines() == [
            "shutil.rmtree('/data/wf_0')",
            "shutil.rmtree('/data/wf_1')",
        ]

    def test_a_workflow_whose_export_was_skipped_is_absent_from_the_helper_script(
        self, install_toolkit, in_tmp_cwd
    ):
        (in_tmp_cwd / "wf_0.json").write_text("blocking")
        toolkit = _FakeWorkflowToolkit(documents=[_FakeDocument("wf_0")])
        install_toolkit(toolkit)
        workflow_delete(_deleteArgs(Export=True))
        assert (in_tmp_cwd / "completeDelete.py").read_text() == ""

    @pytest.mark.xfail(
        strict=True,
        reason="B131: workflow_delete writes its cleanup script to "
               "'completeDelete.py' but prints 'python completeRemove.py'; "
               "following the instruction fails with 'No such file'. "
               "See the consolidated findings issue.",
    )
    def test_the_script_name_it_prints_is_the_one_it_wrote(
        self, install_toolkit, in_tmp_cwd, capsys
    ):
        install_toolkit(_FakeWorkflowToolkit(documents=[_FakeDocument("wf_0")]))
        workflow_delete(_deleteArgs())
        printed = capsys.readouterr().out
        assert "completeDelete.py" in printed

    def test_it_currently_points_the_user_at_a_file_that_does_not_exist(
        self, install_toolkit, in_tmp_cwd, capsys
    ):
        """Characterisation of B131."""
        install_toolkit(_FakeWorkflowToolkit(documents=[_FakeDocument("wf_0")]))
        workflow_delete(_deleteArgs())
        printed = capsys.readouterr().out
        assert "python completeRemove.py" in printed
        assert (in_tmp_cwd / "completeDelete.py").exists()
        assert not (in_tmp_cwd / "completeRemove.py").exists()


# ---------------------------------------------------------------------------
# export
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestWorkflowExport:
    def test_one_json_file_is_written_per_workflow(self, install_toolkit, in_tmp_cwd):
        toolkit = _FakeWorkflowToolkit(
            documents=[_FakeDocument("wf_0", workflow={"nodes": {"A": 1}}),
                       _FakeDocument("wf_1", workflow={"nodes": {"B": 2}})]
        )
        install_toolkit(toolkit)
        workflow_export(Namespace(projectName=PROJECT, workflows=["wf_0", "wf_1"],
                                  forceOverwrite=False))
        assert json.loads((in_tmp_cwd / "wf_0.json").read_text()) == {"nodes": {"A": 1}}
        assert json.loads((in_tmp_cwd / "wf_1.json").read_text()) == {"nodes": {"B": 2}}

    def test_an_existing_file_is_left_alone_and_the_user_is_told(
        self, install_toolkit, in_tmp_cwd, capsys
    ):
        (in_tmp_cwd / "wf_0.json").write_text("mine")
        install_toolkit(_FakeWorkflowToolkit(documents=[_FakeDocument("wf_0")]))
        workflow_export(Namespace(projectName=PROJECT, workflows=["wf_0"],
                                  forceOverwrite=False))
        assert (in_tmp_cwd / "wf_0.json").read_text() == "mine"
        assert "exists in current directory" in capsys.readouterr().out

    def test_forceoverwrite_replaces_the_existing_file(self, install_toolkit, in_tmp_cwd, capsys):
        (in_tmp_cwd / "wf_0.json").write_text("mine")
        install_toolkit(
            _FakeWorkflowToolkit(documents=[_FakeDocument("wf_0", workflow={"nodes": {}})])
        )
        workflow_export(Namespace(projectName=PROJECT, workflows=["wf_0"],
                                  forceOverwrite=True))
        assert json.loads((in_tmp_cwd / "wf_0.json").read_text()) == {"nodes": {}}
        assert capsys.readouterr().out == ""


# ---------------------------------------------------------------------------
# compareToDisk
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestWorkflowCompareToDisk:
    def test_the_whole_list_is_looked_up_without_stopping_at_the_first_match(
        self, install_toolkit, tmp_path
    ):
        toolkit = _FakeWorkflowToolkit(
            hermesWorkflow=_FakeHermesWorkflow(name="wf_0",
                                               resource=str(tmp_path / "missing.json"))
        )
        install_toolkit(toolkit)
        workflow_compareToDisk(Namespace(projectName=PROJECT, workflows=["wf_0", "wf_1"]))
        name, args, kwargs = toolkit.call("getHermesWorkflowFromDB")
        assert args == (["wf_0", "wf_1"],)
        assert kwargs == dict(returnFirst=False)

    def test_a_workflow_that_is_not_on_the_disk_is_reported(self, install_toolkit, tmp_path, capsys):
        install_toolkit(
            _FakeWorkflowToolkit(
                hermesWorkflow=_FakeHermesWorkflow(name="wf_0",
                                                   resource=str(tmp_path / "missing.json"))
            )
        )
        workflow_compareToDisk(Namespace(projectName=PROJECT, workflows=["wf_0"]))
        printed = capsys.readouterr().out
        assert "does not exist on the disk" in printed
        assert "use export to create it" in printed

    def test_an_unchanged_workflow_reports_that_disk_and_db_are_identical(
        self, install_toolkit, tmp_path, capsys
    ):
        onDisk = tmp_path / "wf_0.json"
        onDisk.write_text(json.dumps({"nodes": {"A": {"value": 1}}}))
        parameters = {"nodes": {"A": {"value": 1}}}
        install_toolkit(
            _FakeWorkflowToolkit(
                hermesWorkflow=_FakeHermesWorkflow(name="wf_0", resource=str(onDisk),
                                                   parametersJSON=parameters),
                localWorkflow=_FakeHermesWorkflow(name="Local", resource=str(onDisk),
                                                  parametersJSON=parameters),
            )
        )
        workflow_compareToDisk(Namespace(projectName=PROJECT, workflows=["wf_0"]))
        printed = capsys.readouterr().out
        assert "Simulation wf_0" in printed
        assert "Disk and DB are identical" in printed

    def test_a_changed_parameter_is_printed_as_a_comparison_table(
        self, install_toolkit, tmp_path, capsys
    ):
        onDisk = tmp_path / "wf_0.json"
        onDisk.write_text(json.dumps({"nodes": {"A": {"value": 2}}}))
        install_toolkit(
            _FakeWorkflowToolkit(
                hermesWorkflow=_FakeHermesWorkflow(
                    name="wf_0", resource=str(onDisk),
                    parametersJSON={"nodes": {"A": {"value": 1}}}),
                localWorkflow=_FakeHermesWorkflow(
                    name="Local", resource=str(onDisk),
                    parametersJSON={"nodes": {"A": {"value": 2}}}),
            )
        )
        workflow_compareToDisk(Namespace(projectName=PROJECT, workflows=["wf_0"]))
        printed = capsys.readouterr().out
        assert "Disk and DB are identical" not in printed
        assert "LocalFile" in printed and "DB" in printed


# ---------------------------------------------------------------------------
# sync
# ---------------------------------------------------------------------------

def _syncArgs(**overrides):
    arguments = dict(projectName=PROJECT, workflows=["wf_0"], force=False, quiet=False)
    arguments.update(overrides)
    return Namespace(**arguments)


def _syncToolkit(onDiskPath, dbParameters, diskParameters, documents=None):
    return _FakeWorkflowToolkit(
        documents=documents if documents is not None else [_FakeDocument("wf_0")],
        hermesWorkflow=_FakeHermesWorkflow(name="wf_0", resource=str(onDiskPath),
                                           parametersJSON=dbParameters),
        localWorkflow=_FakeHermesWorkflow(name="Local", resource=str(onDiskPath),
                                          parametersJSON=diskParameters),
    )


@pytest.mark.unit
class TestWorkflowSyncToDb:
    def test_a_workflow_that_is_not_on_the_disk_is_reported_and_not_synced(
        self, install_toolkit, tmp_path, capsys
    ):
        toolkit = _syncToolkit(tmp_path / "missing.json", {}, {})
        install_toolkit(toolkit)
        workflow_sync_to_db(_syncArgs())
        assert "does not exist on the disk" in capsys.readouterr().out
        assert "updateDocumentWorkflow" not in toolkit.names()

    def test_identical_parameters_are_reported_and_nothing_is_written_to_the_db(
        self, install_toolkit, tmp_path, capsys
    ):
        onDisk = tmp_path / "wf_0.json"
        onDisk.write_text(json.dumps({"nodes": {"A": {"value": 1}}}))
        parameters = {"nodes": {"A": {"value": 1}}}
        toolkit = _syncToolkit(onDisk, parameters, parameters)
        install_toolkit(toolkit)
        workflow_sync_to_db(_syncArgs())
        printed = capsys.readouterr().out
        assert "Simulation wf_0" in printed
        assert "Disk and DB parameters are identical" in printed
        assert "updateDocumentWorkflow" not in toolkit.names()

    def test_a_local_change_is_pushed_to_the_database(self, install_toolkit, tmp_path, capsys):
        onDisk = tmp_path / "wf_0.json"
        onDisk.write_text(json.dumps({"nodes": {"A": {"value": 2}}}))
        document = _FakeDocument("wf_0")
        toolkit = _syncToolkit(onDisk, {"nodes": {"A": {"value": 1}}},
                               {"nodes": {"A": {"value": 2}}}, documents=[document])
        install_toolkit(toolkit)
        workflow_sync_to_db(_syncArgs())
        assert "Found Changes" in capsys.readouterr().out
        name, args, kwargs = toolkit.call("updateDocumentWorkflow")
        assert kwargs["document"] is document
        assert kwargs["workflow"] is toolkit.localWorkflow

    def test_the_force_flag_syncs_even_when_nothing_changed(self, install_toolkit, tmp_path):
        onDisk = tmp_path / "wf_0.json"
        onDisk.write_text(json.dumps({"nodes": {"A": {"value": 1}}}))
        parameters = {"nodes": {"A": {"value": 1}}}
        toolkit = _syncToolkit(onDisk, parameters, parameters)
        install_toolkit(toolkit)
        workflow_sync_to_db(_syncArgs(force=True))
        assert "updateDocumentWorkflow" in toolkit.names()

    def test_the_quiet_flag_suppresses_the_report_but_still_syncs(
        self, install_toolkit, tmp_path, capsys
    ):
        onDisk = tmp_path / "wf_0.json"
        onDisk.write_text(json.dumps({"nodes": {"A": {"value": 2}}}))
        toolkit = _syncToolkit(onDisk, {"nodes": {"A": {"value": 1}}},
                               {"nodes": {"A": {"value": 2}}})
        install_toolkit(toolkit)
        workflow_sync_to_db(_syncArgs(quiet=True))
        assert capsys.readouterr().out == ""
        assert "updateDocumentWorkflow" in toolkit.names()

    def test_an_unknown_workflow_prints_nothing_and_syncs_nothing(
        self, install_toolkit, capsys
    ):
        toolkit = _FakeWorkflowToolkit(documents=[])
        install_toolkit(toolkit)
        workflow_sync_to_db(_syncArgs(workflows=["nosuchflow"]))
        assert capsys.readouterr().out == ""
        assert toolkit.names() == ["getWorkflowListDocumentFromDB"]


# ---------------------------------------------------------------------------
# variations
# ---------------------------------------------------------------------------

def _variationArgs(variationJson, **overrides):
    arguments = dict(projectName=PROJECT, workflow="wf", variationJson=str(variationJson),
                     dry_run=True, overwrite=False, naming_scheme="index")
    arguments.update(overrides)
    return Namespace(**arguments)


def _variationToolkit(workflowDirectory):
    document = _FakeDocument(
        "wf_0",
        workflow={"nodes": {"A": {"value": 1}}},
        groupName="wf",
        resource=str(workflowDirectory / "wf_0"),
    )
    return _FakeWorkflowToolkit(documents=[document])


VARIATION_JSON = [{"$.nodes.A.value": [10, 20]}]


@pytest.mark.unit
class TestCreateWorkflowVariations:
    def test_a_dry_run_writes_one_file_per_variation_without_registering_them(
        self, install_toolkit, tmp_path
    ):
        variations = tmp_path / "variations.json"
        variations.write_text(json.dumps(VARIATION_JSON))
        toolkit = _variationToolkit(tmp_path)
        install_toolkit(toolkit)
        create_workflow_variations(_variationArgs(variations))
        assert json.loads((tmp_path / "wf_0.json").read_text()) == {
            "nodes": {"A": {"value": 10}}
        }
        assert json.loads((tmp_path / "wf_1.json").read_text()) == {
            "nodes": {"A": {"value": 20}}
        }
        assert "addWorkflowFileInGroup" not in toolkit.names()

    def test_a_wet_run_registers_every_variation_file_it_wrote(self, install_toolkit, tmp_path):
        variations = tmp_path / "variations.json"
        variations.write_text(json.dumps(VARIATION_JSON))
        toolkit = _variationToolkit(tmp_path)
        install_toolkit(toolkit)
        create_workflow_variations(_variationArgs(variations, dry_run=False))
        added = [call[1][0] for call in toolkit.calls if call[0] == "addWorkflowFileInGroup"]
        assert added == [str(tmp_path / "wf_0.json"), str(tmp_path / "wf_1.json")]

    def test_an_existing_variation_file_is_kept_when_overwrite_is_not_given(
        self, install_toolkit, tmp_path
    ):
        variations = tmp_path / "variations.json"
        variations.write_text(json.dumps(VARIATION_JSON))
        (tmp_path / "wf_0.json").write_text("hand edited")
        install_toolkit(_variationToolkit(tmp_path))
        create_workflow_variations(_variationArgs(variations))
        assert (tmp_path / "wf_0.json").read_text() == "hand edited"
        assert (tmp_path / "wf_1.json").exists()

    def test_overwrite_replaces_an_existing_variation_file(self, install_toolkit, tmp_path):
        variations = tmp_path / "variations.json"
        variations.write_text(json.dumps(VARIATION_JSON))
        (tmp_path / "wf_0.json").write_text("hand edited")
        install_toolkit(_variationToolkit(tmp_path))
        create_workflow_variations(_variationArgs(variations, overwrite=True))
        assert json.loads((tmp_path / "wf_0.json").read_text()) == {
            "nodes": {"A": {"value": 10}}
        }

    def test_a_wet_run_with_another_naming_scheme_refuses_to_write_anything(
        self, install_toolkit, tmp_path
    ):
        variations = tmp_path / "variations.json"
        variations.write_text(json.dumps(VARIATION_JSON))
        toolkit = _variationToolkit(tmp_path)
        install_toolkit(toolkit)
        create_workflow_variations(
            _variationArgs(variations, dry_run=False, naming_scheme="name")
        )
        assert not (tmp_path / "wf_0.json").exists()
        assert "addWorkflowFileInGroup" not in toolkit.names()

    def test_a_dry_run_with_an_unimplemented_naming_scheme_refuses_to_write_anything(
        self, install_toolkit, tmp_path
    ):
        variations = tmp_path / "variations.json"
        variations.write_text(json.dumps(VARIATION_JSON))
        install_toolkit(_variationToolkit(tmp_path))
        create_workflow_variations(_variationArgs(variations, naming_scheme="name"))
        assert not (tmp_path / "wf_0.json").exists()

    @pytest.mark.xfail(
        strict=True,
        reason="B133: create_workflow_variations only logs an error when the "
               "variation JSON path does not exist and then falls through to "
               "open() it, so the command dies with a raw FileNotFoundError "
               "instead of reporting the bad path. "
               "See the consolidated findings issue.",
    )
    def test_a_missing_variation_file_is_reported_without_crashing(
        self, install_toolkit, tmp_path
    ):
        install_toolkit(_variationToolkit(tmp_path))
        create_workflow_variations(_variationArgs(tmp_path / "nosuchfile.json"))

    def test_a_missing_variation_file_currently_raises_filenotfounderror(
        self, install_toolkit, tmp_path
    ):
        """Characterisation of B133."""
        install_toolkit(_variationToolkit(tmp_path))
        with pytest.raises(FileNotFoundError):
            create_workflow_variations(_variationArgs(tmp_path / "nosuchfile.json"))

    @pytest.mark.xfail(
        strict=True,
        reason="B134: create_workflow_variations logs 'workflow not found.' "
               "when the DB returns no document and then immediately does "
               "workflow_doc_list[0], so an unknown --workflow dies with "
               "IndexError instead of the error message. "
               "See the consolidated findings issue.",
    )
    def test_an_unknown_base_workflow_is_reported_without_crashing(
        self, install_toolkit, tmp_path
    ):
        variations = tmp_path / "variations.json"
        variations.write_text(json.dumps(VARIATION_JSON))
        install_toolkit(_FakeWorkflowToolkit(documents=[]))
        create_workflow_variations(_variationArgs(variations, workflow="nosuchflow"))

    def test_an_unknown_base_workflow_currently_raises_indexerror(
        self, install_toolkit, tmp_path
    ):
        """Characterisation of B134."""
        variations = tmp_path / "variations.json"
        variations.write_text(json.dumps(VARIATION_JSON))
        install_toolkit(_FakeWorkflowToolkit(documents=[]))
        with pytest.raises(IndexError):
            create_workflow_variations(_variationArgs(variations, workflow="nosuchflow"))


# ---------------------------------------------------------------------------
# batch delete
# ---------------------------------------------------------------------------

def _batchArgs(**overrides):
    arguments = dict(projectName=PROJECT, workflow=None, groupName=None, hard=False)
    arguments.update(overrides)
    return Namespace(**arguments)


@pytest.mark.unit
class TestBatchDeleteWorkflows:
    def test_the_named_group_is_deleted_as_given(self, install_toolkit):
        toolkit = _FakeWorkflowToolkit()
        install_toolkit(toolkit)
        batch_delete_workflows(_batchArgs(groupName="myGroup"))
        name, args, kwargs = toolkit.call("deleteWorkflowInGroup")
        assert args == ("myGroup", False)
        assert kwargs["resetCounter"] is True
        assert kwargs["exclude"] == [None]
        assert "getWorkflowListDocumentFromDB" not in toolkit.names()

    def test_the_group_of_the_named_workflow_is_deleted_and_that_workflow_kept(
        self, install_toolkit
    ):
        toolkit = _FakeWorkflowToolkit(
            documents=[_FakeDocument("wf_3", groupName="inferredGroup")]
        )
        install_toolkit(toolkit)
        batch_delete_workflows(_batchArgs(workflow="wf_3"))
        name, args, kwargs = toolkit.call("deleteWorkflowInGroup")
        assert args == ("inferredGroup", False)
        assert kwargs["exclude"] == ["wf_3"]

    def test_an_explicit_group_wins_over_the_group_of_the_named_workflow(self, install_toolkit):
        toolkit = _FakeWorkflowToolkit(
            documents=[_FakeDocument("wf_3", groupName="inferredGroup")]
        )
        install_toolkit(toolkit)
        batch_delete_workflows(_batchArgs(workflow="wf_3", groupName="explicitGroup"))
        assert toolkit.call("deleteWorkflowInGroup")[1] == ("explicitGroup", False)

    def test_the_hard_flag_is_forwarded_as_the_deep_delete_argument(self, install_toolkit):
        toolkit = _FakeWorkflowToolkit()
        install_toolkit(toolkit)
        batch_delete_workflows(_batchArgs(groupName="myGroup", hard=True))
        assert toolkit.call("deleteWorkflowInGroup")[1] == ("myGroup", True)

    def test_an_unknown_workflow_raises_a_value_error(self, install_toolkit):
        toolkit = _FakeWorkflowToolkit(documents=[])
        install_toolkit(toolkit)
        with pytest.raises(ValueError, match="workflow not found"):
            batch_delete_workflows(_batchArgs(workflow="nosuchflow"))
        assert "deleteWorkflowInGroup" not in toolkit.names()

    def test_neither_a_group_nor_a_workflow_raises_a_value_error(self, install_toolkit):
        toolkit = _FakeWorkflowToolkit()
        install_toolkit(toolkit)
        with pytest.raises(ValueError, match="either group name or workflow name"):
            batch_delete_workflows(_batchArgs())
        assert "deleteWorkflowInGroup" not in toolkit.names()


# ---------------------------------------------------------------------------
# list workflows
# ---------------------------------------------------------------------------

def _listArgs(**overrides):
    arguments = dict(projectName=PROJECT, group="wf", nodes=False, parameters=False)
    arguments.update(overrides)
    return Namespace(**arguments)


@pytest.mark.unit
class TestWorkflowList:
    def test_the_workflow_names_of_the_group_are_printed_under_a_title(
        self, install_toolkit, capsys
    ):
        toolkit = _FakeWorkflowToolkit(
            documents=[_FakeDocument("wf_0", groupName="wf")],
            workflows=[{"workflowName": "wf_0", "nodes": ["A"], "parameters": {}},
                       {"workflowName": "wf_1", "nodes": ["A"], "parameters": {}}],
        )
        install_toolkit(toolkit)
        workflow_list(_listArgs())
        printed = capsys.readouterr().out
        assert "The simulations in group *wf*  in project *UNIT_WORKFLOW_PROJECT*" in printed
        assert "* wf_0" in printed and "* wf_1" in printed
        assert "+ A" not in printed
        assert toolkit.call("listWorkflows")[2] == dict(
            workflowGroup="wf", listNodes=False, listParameters=False
        )

    def test_the_group_name_comes_from_the_document_and_not_from_the_argument(
        self, install_toolkit, capsys
    ):
        toolkit = _FakeWorkflowToolkit(
            documents=[_FakeDocument("wf_0", groupName="realGroup")],
            workflows=[{"workflowName": "wf_0", "nodes": [], "parameters": {}}],
        )
        install_toolkit(toolkit)
        workflow_list(_listArgs(group="wf_0"))
        assert "group *realGroup*" in capsys.readouterr().out
        assert toolkit.call("listWorkflows")[2]["workflowGroup"] == "realGroup"

    def test_the_nodes_flag_lists_the_nodes_of_every_workflow(self, install_toolkit, capsys):
        toolkit = _FakeWorkflowToolkit(
            documents=[_FakeDocument("wf_0", groupName="wf")],
            workflows=[{"workflowName": "wf_0", "nodes": ["A", "B"], "parameters": {}}],
        )
        install_toolkit(toolkit)
        workflow_list(_listArgs(nodes=True))
        printed = capsys.readouterr().out
        assert "+ A" in printed and "+ B" in printed
        assert toolkit.call("listWorkflows")[2]["listNodes"] is True

    def test_the_parameters_flag_lists_the_parameters_of_every_node(
        self, install_toolkit, capsys
    ):
        toolkit = _FakeWorkflowToolkit(
            documents=[_FakeDocument("wf_0", groupName="wf")],
            workflows=[{"workflowName": "wf_0", "nodes": ["A"],
                        "parameters": {"A": {"dx": 1, "dy": 2}}}],
        )
        install_toolkit(toolkit)
        workflow_list(_listArgs(parameters=True))
        printed = capsys.readouterr().out
        assert "+ A" in printed
        assert "- dx" in printed and "- dy" in printed
        assert toolkit.call("listWorkflows")[2]["listParameters"] is True

    @pytest.mark.xfail(
        strict=True,
        reason="B132: workflow_list's 'nothing found' branch prints "
               "arguments.object, which the `list workflows` parser never "
               "defines (it declares group/projectName/nodes/parameters), "
               "and then indexes simDocument[0] anyway -- so an unknown "
               "group dies with AttributeError instead of the message. "
               "See the consolidated findings issue.",
    )
    def test_an_unknown_group_is_reported_instead_of_crashing(self, install_toolkit, capsys):
        install_toolkit(_FakeWorkflowToolkit(documents=[]))
        workflow_list(_listArgs(group="nosuchgroup"))
        assert "is not a simulation" in capsys.readouterr().out

    def test_an_unknown_group_currently_raises_attributeerror_on_arguments_object(
        self, install_toolkit
    ):
        """Characterisation of B132."""
        install_toolkit(_FakeWorkflowToolkit(documents=[]))
        with pytest.raises(AttributeError, match="object"):
            workflow_list(_listArgs(group="nosuchgroup"))

    def test_even_with_an_object_attribute_an_unknown_group_still_raises_indexerror(
        self, install_toolkit, capsys
    ):
        """Characterisation of B132: the missing attribute hides a second
        fault -- the handler never returns after reporting, so it indexes
        the empty document list."""
        install_toolkit(_FakeWorkflowToolkit(documents=[]))
        with pytest.raises(IndexError):
            workflow_list(_listArgs(group="nosuchgroup", object="nosuchgroup"))
        assert "is not a simulation" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# list nodes / node parameters
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestWorkflowNodesList:
    def test_a_workflow_in_the_database_is_looked_up_by_name(self, install_toolkit, capsys):
        toolkit = _FakeWorkflowToolkit(
            hermesWorkflow=_FakeHermesWorkflow(nodes={"A": {"dx": 1}, "B": {}})
        )
        asked = install_toolkit(toolkit)
        workflowNodes_list(
            Namespace(projectName=PROJECT, workflowName="wf_0", parameters=False)
        )
        printed = capsys.readouterr().out
        assert "The nodes of the wf_0" in printed
        assert "* A" in printed and "* B" in printed
        assert asked["kwargs"]["projectName"] == PROJECT
        assert toolkit.call("getHermesWorkflowFromDB")[1] == ("wf_0",)

    def test_a_non_file_workflow_without_a_project_name_raises_a_value_error(
        self, install_toolkit
    ):
        install_toolkit(_FakeWorkflowToolkit())
        with pytest.raises(ValueError, match="Must supply a project name"):
            workflowNodes_list(
                Namespace(projectName=None, workflowName="wf_0", parameters=False)
            )

    def test_a_directory_is_looked_up_in_the_database_rather_than_read(
        self, install_toolkit, tmp_path, capsys
    ):
        directory = tmp_path / "wf_0"
        directory.mkdir()
        toolkit = _FakeWorkflowToolkit(hermesWorkflow=_FakeHermesWorkflow(nodes={"A": {}}))
        install_toolkit(toolkit)
        workflowNodes_list(
            Namespace(projectName=PROJECT, workflowName=str(directory), parameters=False)
        )
        assert "* A" in capsys.readouterr().out
        assert "getHermesWorkflowFromDB" in toolkit.names()

    def test_the_parameters_flag_lists_each_node_with_its_parameter_names(
        self, install_toolkit, capsys
    ):
        install_toolkit(
            _FakeWorkflowToolkit(
                hermesWorkflow=_FakeHermesWorkflow(nodes={"A": {"dx": 1, "dy": 2}})
            )
        )
        workflowNodes_list(
            Namespace(projectName=PROJECT, workflowName="wf_0", parameters=True)
        )
        printed = capsys.readouterr().out
        assert "* A" in printed
        assert "+ dx" in printed and "+ dy" in printed

    def test_a_workflow_file_on_the_disk_is_read_from_the_disk(
        self, install_toolkit, fake_hermes, tmp_path, capsys
    ):
        onDisk = tmp_path / "wf_0.json"
        onDisk.write_text(json.dumps({"nodes": {}}))
        fake_hermes["result"] = _FakeHermesWorkflow(nodes={"A": {}, "B": {}})
        toolkit = _FakeWorkflowToolkit()
        install_toolkit(toolkit)
        workflowNodes_list(
            Namespace(projectName=None, workflowName=str(onDisk), parameters=False)
        )
        printed = capsys.readouterr().out
        assert "* A" in printed and "* B" in printed
        assert fake_hermes["Resource_path"] == str(onDisk)
        assert fake_hermes["json"] == {"nodes": {}}
        assert toolkit.names() == []


@pytest.mark.unit
class TestWorkflowNodesListParameters:
    def test_the_parameters_of_a_node_in_the_database_are_printed_as_json(
        self, install_toolkit, capsys
    ):
        install_toolkit(
            _FakeWorkflowToolkit(
                hermesWorkflow=_FakeHermesWorkflow(nodes={"A": {"dx": [1, 2]}})
            )
        )
        workflowNodes_listParameters(
            Namespace(projectName=PROJECT, workflowName="wf_0", nodeName="A")
        )
        printed = capsys.readouterr().out
        assert "The parameters of the node A in the workflow wf_0" in printed
        assert "-\t dx:" in printed
        assert json.dumps([1, 2], indent=4) in printed

    def test_an_unknown_node_raises_a_value_error_listing_the_known_nodes(
        self, install_toolkit
    ):
        install_toolkit(
            _FakeWorkflowToolkit(hermesWorkflow=_FakeHermesWorkflow(nodes={"A": {}, "B": {}}))
        )
        with pytest.raises(ValueError, match="Existing nodes are: A,B"):
            workflowNodes_listParameters(
                Namespace(projectName=PROJECT, workflowName="wf_0", nodeName="Z")
            )

    def test_a_workflow_file_on_the_disk_is_read_from_the_disk(
        self, install_toolkit, fake_hermes, tmp_path, capsys
    ):
        onDisk = tmp_path / "wf_0.json"
        onDisk.write_text(json.dumps({"nodes": {}}))
        fake_hermes["result"] = _FakeHermesWorkflow(nodes={"A": {"dx": 1}})
        toolkit = _FakeWorkflowToolkit()
        install_toolkit(toolkit)
        workflowNodes_listParameters(
            Namespace(projectName=None, workflowName=str(onDisk), nodeName="A")
        )
        assert "-\t dx:  1" in capsys.readouterr().out
        assert toolkit.names() == []


# ---------------------------------------------------------------------------
# compare
# ---------------------------------------------------------------------------

def _comparison():
    return pandas.DataFrame({"wf_0": [1], "wf_1": [2]}, index=["nodes.A.value"])


def _compareArgs(**overrides):
    arguments = dict(projectName=PROJECT, workflows=["wf_0", "wf_1"], longFormat=False,
                     transpose=False, format="pandas", file=None)
    arguments.update(overrides)
    return Namespace(**arguments)


@pytest.mark.unit
class TestWorkflowCompare:
    def test_the_comparison_arguments_are_forwarded_to_the_toolkit(self, install_toolkit, capsys):
        toolkit = _FakeWorkflowToolkit(comparison=_comparison())
        install_toolkit(toolkit)
        workflow_compare(_compareArgs(longFormat=True, transpose=True))
        capsys.readouterr()
        name, args, kwargs = toolkit.call("compareWorkflows")
        assert args == (["wf_0", "wf_1"],)
        assert kwargs == dict(longFormat=True, transpose=True)

    def test_the_default_pandas_format_prints_the_table(self, install_toolkit, capsys):
        install_toolkit(_FakeWorkflowToolkit(comparison=_comparison()))
        workflow_compare(_compareArgs())
        printed = capsys.readouterr().out
        assert "nodes.A.value" in printed and "wf_0" in printed

    def test_the_csv_format_is_written_to_the_requested_file(
        self, install_toolkit, in_tmp_cwd, capsys
    ):
        install_toolkit(_FakeWorkflowToolkit(comparison=_comparison()))
        workflow_compare(_compareArgs(format="csv", file="out"))
        capsys.readouterr()
        assert (in_tmp_cwd / "out.csv").read_text() == _comparison().to_csv()

    def test_the_latex_format_is_written_with_a_tex_extension(
        self, install_toolkit, in_tmp_cwd, capsys
    ):
        install_toolkit(_FakeWorkflowToolkit(comparison=_comparison()))
        workflow_compare(_compareArgs(format="latex", file="out"))
        capsys.readouterr()
        assert (in_tmp_cwd / "out.tex").read_text() == _comparison().to_latex()

    def test_the_json_format_is_written_as_indented_json(
        self, install_toolkit, in_tmp_cwd, capsys
    ):
        install_toolkit(_FakeWorkflowToolkit(comparison=_comparison()))
        workflow_compare(_compareArgs(format="json", file="out"))
        capsys.readouterr()
        written = json.loads((in_tmp_cwd / "out.json").read_text())
        assert written == json.loads(_comparison().to_json())

    def test_a_file_name_that_already_has_an_extension_is_used_as_is(
        self, install_toolkit, in_tmp_cwd, capsys
    ):
        install_toolkit(_FakeWorkflowToolkit(comparison=_comparison()))
        workflow_compare(_compareArgs(format="csv", file="report.txt"))
        capsys.readouterr()
        assert (in_tmp_cwd / "report.txt").exists()
        assert not (in_tmp_cwd / "report.txt.csv").exists()

    def test_an_empty_comparison_is_reported_and_no_file_is_written(
        self, install_toolkit, in_tmp_cwd, capsys
    ):
        install_toolkit(_FakeWorkflowToolkit(comparison=pandas.DataFrame()))
        workflow_compare(_compareArgs(format="csv", file="out"))
        assert "Could not found any workflows to compare" in capsys.readouterr().out
        assert not (in_tmp_cwd / "out.csv").exists()

    @pytest.mark.xfail(
        strict=True,
        reason="B135: workflow_compare keeps the raw DataFrame as `output` "
               "for the default --format pandas, then does "
               "outputFile.write(output) -- so `compare --file x` without an "
               "explicit --format dies with TypeError and leaves an empty "
               "file behind. See the consolidated findings issue.",
    )
    def test_the_pandas_format_can_also_be_written_to_a_file(
        self, install_toolkit, in_tmp_cwd, capsys
    ):
        install_toolkit(_FakeWorkflowToolkit(comparison=_comparison()))
        workflow_compare(_compareArgs(file="out"))
        capsys.readouterr()
        assert (in_tmp_cwd / "out.txt").read_text() != ""

    def test_the_pandas_format_with_a_file_currently_raises_typeerror(
        self, install_toolkit, in_tmp_cwd, capsys
    ):
        """Characterisation of B135."""
        install_toolkit(_FakeWorkflowToolkit(comparison=_comparison()))
        with pytest.raises(TypeError):
            workflow_compare(_compareArgs(file="out"))
        capsys.readouterr()
        assert (in_tmp_cwd / "out.txt").read_text() == ""


# ---------------------------------------------------------------------------
# build / execute
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestWorkflowBuildAndExecute:
    def test_build_hands_the_arguments_to_the_hermes_build_handler(self, fake_hermes):
        arguments = Namespace(workflowName="wf_0")
        workflow_build(arguments)
        assert fake_hermes["build"] is arguments
        assert "execute" not in fake_hermes

    def test_execute_hands_the_arguments_to_the_hermes_execute_handler(self, fake_hermes):
        arguments = Namespace(workflowName="wf_0")
        workflow_execute(arguments)
        assert fake_hermes["execute"] is arguments
        assert "build" not in fake_hermes


def _buildExecuteArgs(**overrides):
    arguments = dict(projectName=PROJECT, workflowName="wf_0.json")
    arguments.update(overrides)
    return Namespace(**arguments)


@pytest.mark.unit
class TestWorkflowBuildExecute:
    def test_the_extension_is_stripped_before_the_database_is_queried(self, install_toolkit):
        toolkit = _FakeWorkflowToolkit(documents=[_FakeDocument("wf_0")])
        install_toolkit(toolkit)
        workflow_buildExecute(_buildExecuteArgs())
        assert toolkit.call("getWorkflowListDocumentFromDB")[1] == ("wf_0",)
        assert toolkit.call("executeWorkflowFromDB")[1] == ("wf_0",)

    def test_without_scheduler_arguments_it_defaults_to_the_local_scheduler(
        self, install_toolkit
    ):
        toolkit = _FakeWorkflowToolkit(documents=[_FakeDocument("wf_0")])
        install_toolkit(toolkit)
        workflow_buildExecute(_buildExecuteArgs())
        assert toolkit.call("executeWorkflowFromDB")[2] == dict(
            scheduler="local", schedulerHost=None, schedulerPort=None, dispatch_id=None
        )

    def test_the_scheduler_options_are_forwarded(self, install_toolkit):
        toolkit = _FakeWorkflowToolkit(documents=[_FakeDocument("wf_0")])
        install_toolkit(toolkit)
        workflow_buildExecute(
            _buildExecuteArgs(scheduler="central", scheduler_host="luigi.example",
                              scheduler_port=8082, dispatch_id="abc123")
        )
        assert toolkit.call("executeWorkflowFromDB")[2] == dict(
            scheduler="central", schedulerHost="luigi.example", schedulerPort=8082,
            dispatch_id="abc123"
        )

    def test_a_workflow_missing_from_the_database_but_present_on_disk_is_added_first(
        self, install_toolkit, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "wf_0.json").write_text(json.dumps({"nodes": {}}))
        toolkit = _FakeWorkflowToolkit(documents=[])
        install_toolkit(toolkit)
        workflow_buildExecute(_buildExecuteArgs())
        assert toolkit.call("addWorkflowFileInGroup")[1] == ("wf_0.json",)
        assert toolkit.call("executeWorkflowFromDB")[1] == ("wf_0",)

    def test_a_workflow_that_is_neither_in_the_database_nor_a_file_raises_a_value_error(
        self, install_toolkit, in_tmp_cwd
    ):
        toolkit = _FakeWorkflowToolkit(documents=[])
        install_toolkit(toolkit)
        with pytest.raises(ValueError, match="is not in the DB and not a file"):
            workflow_buildExecute(_buildExecuteArgs())
        assert "executeWorkflowFromDB" not in toolkit.names()

    def test_it_does_not_ask_the_case_configuration_for_a_missing_project_name(
        self, install_toolkit, in_tmp_cwd
    ):
        """Characterisation: unlike every other handler, workflow_buildExecute
        never calls _confirm_project_name, so `--projectName` is mandatory in
        practice -- a None is handed straight to getToolkit."""
        (in_tmp_cwd / "caseConfiguration.json").write_text(
            json.dumps({"projectName": "FROM_CASE_CONFIG"})
        )
        toolkit = _FakeWorkflowToolkit(documents=[_FakeDocument("wf_0")])
        asked = install_toolkit(toolkit)
        workflow_buildExecute(_buildExecuteArgs(projectName=None))
        assert asked["kwargs"]["projectName"] is None
