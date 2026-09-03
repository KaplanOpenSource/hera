"""hermesWorkflowToolkit: the Luigi command builder, the hermes-object
factories, the template datasource lookups, and the whole DB-backed
workflow lifecycle (add -> retrieve -> update -> list -> compare -> delete)
driven against the in-memory datalayer.

The workflows used here are minimal but real hermes workflows: five nodes
named exactly as ``workflow_Eulerian._requiredNodeList`` demands
(``controlDict``, ``fvSolution``, ``fvSchemes``, ``fileWriter``,
``Parameters``), each a ``general.Parameters`` stub, with only the
``Parameters`` node carrying a value.  A ``solver`` is always declared,
because ``getHermesWorkflowFromJSON`` cannot build a solver-less workflow
at all -- see B136.

Deliberately not covered here:

* ``executeWorkflowFromDB`` -- it shells out to ``python3 -m luigi`` via
  subprocess and writes a generated module to disk; the part of it worth
  unit-testing is ``buildLuigiExecutionCommand``, which is covered below.
* The name helpers ``getworkFlowName`` / ``splitWorkflowName`` /
  ``isInFolder`` -- already covered by
  ``test_hermes_workflow_toolkit_names.py``.
* ``actionModes`` -- a bare enum with no behaviour.
"""
import json
import os

import pandas
import pytest

from hera.simulations.hermesWorkflowToolkit import (
    SCHEDULER_CENTRAL,
    SCHEDULER_LOCAL,
    buildLuigiExecutionCommand,
)

_REQUIRED_NODES = ["controlDict", "fvSolution", "fvSchemes", "fileWriter", "Parameters"]


def _workflowJSON(alpha=1.0, solver="simpleFoam"):
    """A minimal hermes workflow whose only varying parameter is Parameters.alpha."""
    nodes = {
        nodeName: {
            "type": "general.Parameters",
            "Execution": {"input_parameters": {}},
        }
        for nodeName in _REQUIRED_NODES
    }
    nodes["Parameters"]["Execution"]["input_parameters"] = {"alpha": alpha}

    workflowBody = {"nodeList": list(_REQUIRED_NODES), "nodes": nodes}
    if solver is not None:
        workflowBody["solver"] = solver
    return {"workflow": workflowBody}


@pytest.fixture()
def tk(unit_toolkit_factory):
    from hera import toolkitHome

    return unit_toolkit_factory(toolkitHome.SIMULATIONS_WORKFLOWS)


@pytest.fixture()
def twoWorkflowsInGroup(tk):
    """Two workflows in the group 'flow', differing only in Parameters.alpha."""
    first = tk.addWorkflowToGroup(_workflowJSON(1.0), "flow")
    second = tk.addWorkflowToGroup(_workflowJSON(9.0), "flow")
    return first, second


# ---------------------------------------------------------------------------
# buildLuigiExecutionCommand
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestBuildLuigiExecutionCommand:
    def test_the_local_scheduler_is_the_default(self):
        command = buildLuigiExecutionCommand("myModule")
        assert command == "python3 -m luigi --module myModule finalnode_xx_0 --local-scheduler"

    def test_the_local_scheduler_does_not_receive_a_dispatch_id(self):
        """dispatch_id only matters to a central scheduler that has to tell
        concurrent runs apart, so it must not reach the local one."""
        command = buildLuigiExecutionCommand("myModule", dispatch_id="deadbeef")
        assert "--dispatch-id" not in command

    def test_the_central_scheduler_replaces_the_local_flag(self):
        command = buildLuigiExecutionCommand("m", scheduler=SCHEDULER_CENTRAL)
        assert "--local-scheduler" not in command

    def test_the_central_scheduler_carries_the_dispatch_id(self):
        command = buildLuigiExecutionCommand(
            "m", dispatch_id="deadbeef", scheduler=SCHEDULER_CENTRAL
        )
        assert command.endswith("--dispatch-id deadbeef")

    def test_an_omitted_central_address_is_left_to_luigis_own_defaults(self):
        command = buildLuigiExecutionCommand("m", scheduler=SCHEDULER_CENTRAL)
        assert "--scheduler-host" not in command
        assert "--scheduler-port" not in command

    def test_a_central_host_and_port_are_both_passed_through(self):
        command = buildLuigiExecutionCommand(
            "m", scheduler=SCHEDULER_CENTRAL, schedulerHost="luigi.example", schedulerPort=9999
        )
        assert "--scheduler-host luigi.example" in command
        assert "--scheduler-port 9999" in command

    def test_the_target_task_is_overridable(self):
        command = buildLuigiExecutionCommand("m", targetTask="myFinalTask")
        assert "--module m myFinalTask" in command

    def test_an_unrecognised_scheduler_name_falls_back_to_local(self):
        """Only the exact string "central" selects the central scheduler."""
        command = buildLuigiExecutionCommand("m", scheduler="CENTRAL")
        assert "--local-scheduler" in command
        assert buildLuigiExecutionCommand("m", scheduler=SCHEDULER_LOCAL) == command


# ---------------------------------------------------------------------------
# getHermesWorkflowFromJSON
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetHermesWorkflowFromJSON:
    def test_the_solver_field_selects_the_matching_openfoam_workflow_class(self, tk):
        workflowObject = tk.getHermesWorkflowFromJSON(_workflowJSON())
        from hera.simulations.openFoam.OFWorkflow import workflow_simpleFoam

        assert isinstance(workflowObject, workflow_simpleFoam)

    def test_the_name_and_resource_are_handed_to_the_workflow_object(self, tk, tmp_path):
        resource = str(tmp_path / "here.json")
        workflowObject = tk.getHermesWorkflowFromJSON(
            _workflowJSON(), name="myFlow", resource=resource
        )
        assert workflowObject.name == "myFlow"
        assert workflowObject.Resource_path == resource

    def test_a_json_string_is_accepted_as_well_as_a_dict(self, tk):
        workflowObject = tk.getHermesWorkflowFromJSON(json.dumps(_workflowJSON(3.0)))
        assert workflowObject.parametersJSON["Parameters"] == {"alpha": 3.0}

    def test_a_json_file_path_is_accepted_as_well_as_a_dict(self, tk, tmp_path):
        path = tmp_path / "wf.json"
        path.write_text(json.dumps(_workflowJSON(4.0)))
        workflowObject = tk.getHermesWorkflowFromJSON(str(path))
        assert workflowObject.parametersJSON["Parameters"] == {"alpha": 4.0}

    def test_an_unknown_solver_is_rejected_by_name(self, tk):
        with pytest.raises(ValueError, match="noSuchSolver"):
            tk.getHermesWorkflowFromJSON(_workflowJSON(solver="noSuchSolver"))

    @pytest.mark.xfail(
        strict=True,
        reason="B136: a workflow with no 'solver' field should fall back to "
               "the generic hermes workflow class, but the fallback is "
               "pydoc.locate('hermes.workflow'), which resolves to the "
               "hermes.workflow *package module* (the class lives at "
               "hermes.workflow.workflow.workflow), so the call raises "
               "TypeError: 'module' object is not callable. "
               "See the consolidated findings issue.",
    )
    def test_a_workflow_without_a_solver_should_use_the_generic_hermes_class(self, tk):
        workflowObject = tk.getHermesWorkflowFromJSON(_workflowJSON(solver=None))
        assert workflowObject.parametersJSON["Parameters"] == {"alpha": 1.0}

    def test_a_workflow_without_a_solver_currently_raises_a_typeerror(self, tk):
        """Characterisation of B136."""
        with pytest.raises(TypeError, match="not callable"):
            tk.getHermesWorkflowFromJSON(_workflowJSON(solver=None))


# ---------------------------------------------------------------------------
# The template datasource lookups
# ---------------------------------------------------------------------------

@pytest.fixture()
def templateDatasources(tk, tmp_path):
    """One 'Flow' template datasource and one 'Node' template datasource."""
    from hera.datalayer import datatypes

    path = tmp_path / "template.json"
    path.write_text(json.dumps({"template": True}))
    tk.addDataSource(
        "myFlowTemplate", str(path), datatypes.JSON_DICT,
        component="Flow", solver="simpleFoam", name="myFlowTemplate",
    )
    tk.addDataSource(
        "myNodeTemplate", str(path), datatypes.JSON_DICT,
        component="Node", solver="simpleFoam", name="myNodeTemplate",
    )
    return path


@pytest.mark.unit
class TestListHermesSolverTemplates:
    def test_a_project_with_no_templates_gives_an_empty_frame(self, tk):
        result = tk.listHermesSolverTemplates("simpleFoam")
        assert isinstance(result, pandas.DataFrame)
        assert result.empty

    def test_the_templates_of_a_solver_are_indexed_by_their_name(self, tk, templateDatasources):
        result = tk.listHermesSolverTemplates("simpleFoam")
        assert result.index.name == "name"
        assert sorted(result.index) == ["myFlowTemplate", "myNodeTemplate"]

    def test_a_different_solver_matches_none_of_them(self, tk, templateDatasources):
        assert tk.listHermesSolverTemplates("buoyantReactingFoam").empty


@pytest.mark.unit
class TestComponentFilteredTemplateLookups:
    """getHermesFlowTemplate / getHermesNodeTemplate / listHermesNodesTemplates
    all filter on the datasource's ``component`` field -- and all pass the
    filter as ``desc__component``, which the datalayer prefixes with ``desc``
    a second time. See B137."""

    def test_the_component_field_is_queryable_when_it_is_not_double_prefixed(
        self, tk, templateDatasources
    ):
        """Baseline: the datasources really are there, and really do carry
        component=Flow / component=Node -- so a failure of the three methods
        below is about the query they build, not about missing data."""
        assert len(tk.getDataSourceDocumentsList(component="Flow")) == 1
        assert len(tk.getDataSourceDocumentsList(component="Node")) == 1

    @pytest.mark.xfail(
        strict=True,
        reason="B137: the component filter is passed as desc__component=..., "
               "but collection.getDocuments feeds every extra kwarg through "
               "dictToMongoQuery(prefix='desc'), producing the query field "
               "desc__desc__component -- which matches nothing. So "
               "getHermesFlowTemplate always returns None. "
               "See the consolidated findings issue.",
    )
    def test_a_registered_flow_template_should_be_retrievable(self, tk, templateDatasources):
        assert tk.getHermesFlowTemplate("myFlowTemplate") == {"template": True}

    def test_a_registered_flow_template_is_currently_never_found(self, tk, templateDatasources):
        """Characterisation of B137."""
        assert tk.getHermesFlowTemplate("myFlowTemplate") is None

    @pytest.mark.xfail(
        strict=True,
        reason="B137: getHermesNodeTemplate filters on desc__component, which "
               "the datalayer re-prefixes into desc__desc__component and so "
               "matches nothing. See the consolidated findings issue.",
    )
    def test_a_registered_node_template_should_be_retrievable(self, tk, templateDatasources):
        assert tk.getHermesNodeTemplate("myNodeTemplate") == {"template": True}

    def test_a_registered_node_template_is_currently_never_found(self, tk, templateDatasources):
        """Characterisation of B137."""
        assert tk.getHermesNodeTemplate("myNodeTemplate") is None

    @pytest.mark.xfail(
        strict=True,
        reason="B137: listHermesNodesTemplates filters on desc__component, "
               "re-prefixed into desc__desc__component, so the listing is "
               "always empty. (Its body would then fail a second time: it "
               "reads doc.desc['desc'], a nesting addDataSource never "
               "creates.) See the consolidated findings issue.",
    )
    def test_a_registered_node_template_should_appear_in_the_node_listing(
        self, tk, templateDatasources
    ):
        assert "myNodeTemplate" in list(tk.listHermesNodesTemplates().index)

    def test_the_node_listing_is_currently_always_empty(self, tk, templateDatasources):
        """Characterisation of B137."""
        assert tk.listHermesNodesTemplates().empty


# ---------------------------------------------------------------------------
# addWorkflowToGroup / addWorkflowFileInGroup
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestAddWorkflowToGroup:
    def test_the_first_workflow_of_a_group_is_numbered_zero(self, tk):
        doc = tk.addWorkflowToGroup(_workflowJSON(), "flow")
        assert doc.desc["workflowName"] == "flow_0000"
        assert doc.desc["groupName"] == "flow"
        assert doc.desc["groupID"] == 0

    def test_the_stored_document_carries_the_solver_workflow_and_parameters(self, tk):
        doc = tk.addWorkflowToGroup(_workflowJSON(2.5), "flow")
        assert doc.desc["solver"] == "simpleFoam"
        assert doc.desc["parameters"]["Parameters"] == {"alpha": 2.5}
        assert doc.desc["workflow"]["workflow"]["nodeList"] == _REQUIRED_NODES

    def test_it_lands_in_the_simulations_collection_as_a_hermes_workflow(self, tk):
        doc = tk.addWorkflowToGroup(_workflowJSON(), "flow")
        assert doc.type == tk.DOCTYPE_WORKFLOW
        assert len(tk.getSimulationsDocuments(type=tk.DOCTYPE_WORKFLOW)) == 1

    def test_a_second_distinct_workflow_gets_the_next_number(self, tk, twoWorkflowsInGroup):
        first, second = twoWorkflowsInGroup
        assert [first.desc["workflowName"], second.desc["workflowName"]] == [
            "flow_0000", "flow_0001",
        ]

    def test_re_adding_an_identical_workflow_returns_the_existing_document(self, tk):
        first = tk.addWorkflowToGroup(_workflowJSON(1.0), "flow")
        again = tk.addWorkflowToGroup(_workflowJSON(1.0), "flow")
        assert again.desc["workflowName"] == first.desc["workflowName"]
        assert len(tk.getWorkflowDocumentsInGroup("flow")) == 1

    def test_an_identical_workflow_is_recognised_even_from_another_group(self, tk):
        """The idempotence check is on the parameters alone, so the requested
        group name is ignored once a matching workflow already exists."""
        first = tk.addWorkflowToGroup(_workflowJSON(1.0), "flow")
        again = tk.addWorkflowToGroup(_workflowJSON(1.0), "otherGroup")
        assert again.desc["groupName"] == first.desc["groupName"] == "flow"

    def test_the_default_resource_is_the_workflow_name_inside_the_files_directory(
        self, tk, unit_files_directory
    ):
        doc = tk.addWorkflowToGroup(_workflowJSON(), "flow")
        assert doc.resource == os.path.join(unit_files_directory, "flow_0000.json")

    def test_an_explicit_resource_is_stored_as_an_absolute_path(self, tk, tmp_path):
        target = tmp_path / "somewhere" / "case"
        doc = tk.addWorkflowToGroup(_workflowJSON(), "flow", resource=str(target))
        assert doc.resource == str(target)

    def test_writing_to_file_produces_the_workflow_json_on_disk(self, tk, unit_files_directory):
        doc = tk.addWorkflowToGroup(_workflowJSON(6.0), "flow", writeWorkflowToFile=True)
        with open(doc.resource) as handle:
            written = json.load(handle)
        assert written["workflow"]["nodes"]["Parameters"]["Execution"]["input_parameters"] == {
            "alpha": 6.0
        }


@pytest.mark.unit
class TestAddWorkflowFileInGroup:
    def test_a_file_inside_the_files_directory_is_registered_under_its_base_name(
        self, tk, unit_files_directory
    ):
        path = os.path.join(unit_files_directory, "myflow_0003.json")
        with open(path, "w") as handle:
            json.dump(_workflowJSON(5.0), handle)

        doc = tk.addWorkflowFileInGroup(path)
        assert doc.desc["groupName"] == "myflow"
        assert doc.resource == path

    def test_the_group_counter_and_not_the_file_name_decides_the_workflow_id(
        self, tk, unit_files_directory
    ):
        """The file is called myflow_0003.json, but only its group prefix is
        used -- the id comes from the (fresh) group counter."""
        path = os.path.join(unit_files_directory, "myflow_0003.json")
        with open(path, "w") as handle:
            json.dump(_workflowJSON(5.0), handle)

        doc = tk.addWorkflowFileInGroup(path)
        assert doc.desc["workflowName"] == "myflow_0000"

    def test_adding_the_same_file_twice_does_not_duplicate_it(self, tk, unit_files_directory):
        path = os.path.join(unit_files_directory, "myflow.json")
        with open(path, "w") as handle:
            json.dump(_workflowJSON(5.0), handle)

        tk.addWorkflowFileInGroup(path)
        tk.addWorkflowFileInGroup(path)
        assert len(tk.getWorkflowDocumentsInGroup("myflow")) == 1

    def test_a_file_outside_the_files_directory_is_refused(self, tk, tmp_path):
        outside = tmp_path / "elsewhere.json"
        outside.write_text(json.dumps(_workflowJSON()))
        with pytest.raises(ValueError, match="is not in"):
            tk.addWorkflowFileInGroup(str(outside))


# ---------------------------------------------------------------------------
# Document retrieval
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetWorkflowDocumentByName:
    def test_an_existing_workflow_is_found_by_its_full_name(self, tk, twoWorkflowsInGroup):
        doc = tk.getWorkflowDocumentByName("flow_0001")
        assert doc.desc["parameters"]["Parameters"] == {"alpha": 9.0}

    def test_an_unknown_name_gives_none(self, tk, twoWorkflowsInGroup):
        assert tk.getWorkflowDocumentByName("flow_9999") is None

    def test_a_group_name_is_not_accepted_as_a_workflow_name(self, tk, twoWorkflowsInGroup):
        """Unlike getWorkflowDocumentFromDB, this method matches the
        workflowName field only."""
        assert tk.getWorkflowDocumentByName("flow") is None

    def test_additional_criteria_narrow_the_search(self, tk, twoWorkflowsInGroup):
        assert tk.getWorkflowDocumentByName("flow_0000", solver="simpleFoam") is not None
        assert tk.getWorkflowDocumentByName("flow_0000", solver="otherFoam") is None

    def test_nothing_is_found_in_the_cache_collection(self, tk, twoWorkflowsInGroup):
        """The workflows were added to Simulations, so the Cache kind is empty."""
        assert tk.getWorkflowDocumentByName("flow_0000", dockind="Cache") is None


@pytest.mark.unit
class TestGetWorkflowDocumentFromDB:
    def test_a_workflow_name_matches_exactly_that_workflow(self, tk, twoWorkflowsInGroup):
        docList = tk.getWorkflowDocumentFromDB("flow_0000")
        assert [doc.desc["workflowName"] for doc in docList] == ["flow_0000"]

    def test_a_group_name_matches_every_workflow_in_the_group(self, tk, twoWorkflowsInGroup):
        docList = tk.getWorkflowDocumentFromDB("flow")
        assert sorted(doc.desc["workflowName"] for doc in docList) == ["flow_0000", "flow_0001"]

    def test_a_resource_path_matches_the_workflow_stored_there(self, tk, twoWorkflowsInGroup):
        first, _ = twoWorkflowsInGroup
        docList = tk.getWorkflowDocumentFromDB(first.resource)
        assert [doc.desc["workflowName"] for doc in docList] == ["flow_0000"]

    def test_a_workflow_json_file_matches_by_its_parameter_values(
        self, tk, twoWorkflowsInGroup, tmp_path
    ):
        path = tmp_path / "lookup.json"
        path.write_text(json.dumps(_workflowJSON(9.0)))
        docList = tk.getWorkflowDocumentFromDB(str(path))
        assert [doc.desc["workflowName"] for doc in docList] == ["flow_0001"]

    def test_a_workflow_dict_matches_by_its_parameter_values(self, tk, twoWorkflowsInGroup):
        docList = tk.getWorkflowDocumentFromDB(_workflowJSON(9.0))
        assert [doc.desc["workflowName"] for doc in docList] == ["flow_0001"]

    def test_a_hermes_workflow_object_matches_by_its_parameter_values(
        self, tk, twoWorkflowsInGroup
    ):
        workflowObject = tk.getHermesWorkflowFromJSON(_workflowJSON(1.0))
        docList = tk.getWorkflowDocumentFromDB(workflowObject)
        assert [doc.desc["workflowName"] for doc in docList] == ["flow_0000"]

    def test_a_workflow_dict_that_matches_nothing_gives_an_empty_result(
        self, tk, twoWorkflowsInGroup
    ):
        assert len(tk.getWorkflowDocumentFromDB(_workflowJSON(-1.0))) == 0

    def test_a_string_that_is_neither_a_name_nor_loadable_json_gives_an_empty_result(
        self, tk, twoWorkflowsInGroup
    ):
        assert len(tk.getWorkflowDocumentFromDB("definitely not a workflow")) == 0

    def test_a_directory_path_gives_an_empty_result(self, tk, unit_files_directory):
        """loadJSON raises IsADirectoryError, which is swallowed."""
        assert len(tk.getWorkflowDocumentFromDB(unit_files_directory)) == 0

    def test_an_unsupported_input_type_gives_an_empty_result(self, tk, twoWorkflowsInGroup):
        assert len(tk.getWorkflowDocumentFromDB(17)) == 0

    def test_the_cache_kind_does_not_see_the_simulation_documents(self, tk, twoWorkflowsInGroup):
        assert len(tk.getWorkflowDocumentFromDB("flow_0000", dockind="Cache")) == 0


@pytest.mark.unit
class TestGetWorkflowListDocumentFromDB:
    def test_a_single_identifier_behaves_like_getworkflowdocumentfromdb(
        self, tk, twoWorkflowsInGroup
    ):
        docList = tk.getWorkflowListDocumentFromDB("flow")
        assert len(docList) == 2

    def test_a_list_of_identifiers_concatenates_their_results(self, tk, twoWorkflowsInGroup):
        docList = tk.getWorkflowListDocumentFromDB(["flow_0000", "flow_0001"])
        assert sorted(doc.desc["workflowName"] for doc in docList) == ["flow_0000", "flow_0001"]

    def test_a_list_entry_that_matches_nothing_contributes_nothing(self, tk, twoWorkflowsInGroup):
        docList = tk.getWorkflowListDocumentFromDB(["flow_0000", "flow_9999"])
        assert [doc.desc["workflowName"] for doc in docList] == ["flow_0000"]


@pytest.mark.unit
class TestGetWorkflowListOfSolvers:
    def test_every_workflow_of_the_solver_is_returned(self, tk, twoWorkflowsInGroup):
        assert len(tk.getWorkflowListOfSolvers("simpleFoam")) == 2

    def test_another_solver_matches_nothing(self, tk, twoWorkflowsInGroup):
        assert len(tk.getWorkflowListOfSolvers("buoyantReactingFoam")) == 0

    def test_additional_criteria_are_applied_on_top_of_the_solver(self, tk, twoWorkflowsInGroup):
        docList = tk.getWorkflowListOfSolvers("simpleFoam", groupName="flow")
        assert len(docList) == 2
        assert len(tk.getWorkflowListOfSolvers("simpleFoam", groupName="other")) == 0


@pytest.mark.unit
class TestGetWorkflowDocumentsInGroup:
    def test_adding_a_workflow_makes_it_appear_in_its_group(self, tk):
        tk.addWorkflowToGroup(_workflowJSON(), "flow")
        docList = tk.getWorkflowDocumentsInGroup("flow")
        assert [doc.desc["workflowName"] for doc in docList] == ["flow_0000"]

    def test_an_unknown_group_is_empty(self, tk, twoWorkflowsInGroup):
        assert len(tk.getWorkflowDocumentsInGroup("noSuchGroup")) == 0

    def test_groups_do_not_leak_into_each_other(self, tk):
        tk.addWorkflowToGroup(_workflowJSON(1.0), "flowA")
        tk.addWorkflowToGroup(_workflowJSON(2.0), "flowB")
        assert len(tk.getWorkflowDocumentsInGroup("flowA")) == 1
        assert len(tk.getWorkflowDocumentsInGroup("flowB")) == 1


# ---------------------------------------------------------------------------
# Document -> hermes object
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetHemresWorkflowFromDocument:
    def test_a_single_document_becomes_a_workflow_named_after_it(self, tk, twoWorkflowsInGroup):
        first, _ = twoWorkflowsInGroup
        workflowObject = tk.getHemresWorkflowFromDocument(first)
        assert workflowObject.name == "flow_0000"
        assert workflowObject.parametersJSON["Parameters"] == {"alpha": 1.0}

    def test_the_documents_resource_becomes_the_workflows_resource_path(
        self, tk, twoWorkflowsInGroup
    ):
        first, _ = twoWorkflowsInGroup
        workflowObject = tk.getHemresWorkflowFromDocument([first])
        assert workflowObject.Resource_path == first.resource

    def test_returning_all_gives_one_workflow_object_per_document(self, tk, twoWorkflowsInGroup):
        first, second = twoWorkflowsInGroup
        workflowObjects = tk.getHemresWorkflowFromDocument([first, second], returnFirst=False)
        assert [w.name for w in workflowObjects] == ["flow_0000", "flow_0001"]

    def test_an_empty_document_list_logs_and_then_raises_indexerror(self, tk):
        """The emptiness is detected and logged, but not turned into a
        graceful return -- the very next statement indexes docList[0]."""
        with pytest.raises(IndexError):
            tk.getHemresWorkflowFromDocument([])


@pytest.mark.unit
class TestGetHermesWorkflowFromDB:
    def test_a_list_of_names_yields_the_first_matching_workflow(self, tk, twoWorkflowsInGroup):
        workflowObject = tk.getHermesWorkflowFromDB(["flow_0001"])
        assert workflowObject.parametersJSON["Parameters"] == {"alpha": 9.0}

    def test_returning_all_yields_every_workflow_of_the_group(self, tk, twoWorkflowsInGroup):
        workflowObjects = tk.getHermesWorkflowFromDB(["flow"], returnFirst=False)
        assert sorted(w.name for w in workflowObjects) == ["flow_0000", "flow_0001"]

    def test_nothing_found_gives_none(self, tk, twoWorkflowsInGroup):
        assert tk.getHermesWorkflowFromDB(["flow_9999"]) is None

    @pytest.mark.xfail(
        strict=True,
        reason="B138: for a non-list identifier getWorkflowListDocumentFromDB "
               "returns the mongoengine QuerySet from getWorkflowDocumentFromDB "
               "unchanged, and getHemresWorkflowFromDocument then tests "
               "isinstance(documentList, list) -- false for a QuerySet -- so it "
               "wraps the whole QuerySet in a one-element list and reads .desc "
               "off it: AttributeError. Only the list-input path works. "
               "See the consolidated findings issue.",
    )
    def test_a_plain_workflow_name_should_also_yield_a_workflow(self, tk, twoWorkflowsInGroup):
        workflowObject = tk.getHermesWorkflowFromDB("flow_0001")
        assert workflowObject.parametersJSON["Parameters"] == {"alpha": 9.0}

    def test_a_plain_workflow_name_currently_raises_an_attributeerror(
        self, tk, twoWorkflowsInGroup
    ):
        """Characterisation of B138."""
        with pytest.raises(AttributeError, match="QuerySet"):
            tk.getHermesWorkflowFromDB("flow_0001")


# ---------------------------------------------------------------------------
# updateDocumentWorkflow
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestUpdateDocumentWorkflow:
    def test_the_new_parameters_are_persisted_to_the_database(self, tk, twoWorkflowsInGroup):
        first, _ = twoWorkflowsInGroup
        replacement = tk.getHermesWorkflowFromJSON(_workflowJSON(42.0))
        tk.updateDocumentWorkflow(first, replacement)

        reread = tk.getWorkflowDocumentByName("flow_0000")
        assert reread.desc["parameters"]["Parameters"] == {"alpha": 42.0}

    def test_the_new_workflow_json_is_persisted_too(self, tk, twoWorkflowsInGroup):
        first, _ = twoWorkflowsInGroup
        replacement = tk.getHermesWorkflowFromJSON(_workflowJSON(42.0))
        tk.updateDocumentWorkflow(first, replacement)

        reread = tk.getWorkflowDocumentByName("flow_0000")
        stored = reread.desc["workflow"]["workflow"]["nodes"]["Parameters"]
        assert stored["Execution"]["input_parameters"] == {"alpha": 42.0}

    def test_the_updated_document_is_then_findable_by_its_new_parameters(
        self, tk, twoWorkflowsInGroup
    ):
        first, _ = twoWorkflowsInGroup
        replacement = tk.getHermesWorkflowFromJSON(_workflowJSON(42.0))
        tk.updateDocumentWorkflow(first, replacement)

        docList = tk.getWorkflowDocumentFromDB(_workflowJSON(42.0))
        assert [doc.desc["workflowName"] for doc in docList] == ["flow_0000"]

    def test_the_name_and_group_of_the_document_are_left_alone(self, tk, twoWorkflowsInGroup):
        first, _ = twoWorkflowsInGroup
        tk.updateDocumentWorkflow(first, tk.getHermesWorkflowFromJSON(_workflowJSON(42.0)))
        reread = tk.getWorkflowDocumentByName("flow_0000")
        assert reread.desc["groupName"] == "flow"
        assert reread.desc["groupID"] == 0


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestCompareWorkflowObj:
    def test_the_differing_parameter_is_reported_per_workflow(self, tk):
        first = tk.getHermesWorkflowFromJSON(_workflowJSON(1.0), name="w1")
        second = tk.getHermesWorkflowFromJSON(_workflowJSON(9.0), name="w2")
        table = tk.compareWorkflowObj([first, second])
        assert list(table.columns) == ["w1", "w2"]
        assert table.loc["Parameters_alpha"].tolist() == [1.0, 9.0]

    def test_the_parameter_path_has_its_dots_replaced_by_underscores(self, tk):
        """So the result can be handed to DataFrame.query, which rejects dots."""
        first = tk.getHermesWorkflowFromJSON(_workflowJSON(1.0), name="w1")
        second = tk.getHermesWorkflowFromJSON(_workflowJSON(9.0), name="w2")
        assert list(tk.compareWorkflowObj([first, second]).index) == ["Parameters_alpha"]

    def test_identical_workflows_have_nothing_to_report(self, tk):
        first = tk.getHermesWorkflowFromJSON(_workflowJSON(1.0), name="w1")
        second = tk.getHermesWorkflowFromJSON(_workflowJSON(1.0), name="w2")
        assert tk.compareWorkflowObj([first, second]).empty

    @pytest.mark.xfail(
        strict=True,
        reason="B139: compareWorkflowObj always passes "
               "changeDotToUnderscore=True, and compareDataframeConfigurations "
               "combines that with longFormat by calling .replace('.','_') on "
               "ret.T.columns -- which in long format is the integer row index, "
               "so every longFormat=True comparison dies with "
               "AttributeError: 'int' object has no attribute 'replace'. "
               "See the consolidated findings issue.",
    )
    def test_the_long_format_should_give_one_row_per_workflow(self, tk):
        first = tk.getHermesWorkflowFromJSON(_workflowJSON(1.0), name="w1")
        second = tk.getHermesWorkflowFromJSON(_workflowJSON(9.0), name="w2")
        table = tk.compareWorkflowObj([first, second], longFormat=True)
        assert sorted(table["datasetName"]) == ["w1", "w2"]

    def test_the_long_format_currently_raises_an_attributeerror(self, tk):
        """Characterisation of B139."""
        first = tk.getHermesWorkflowFromJSON(_workflowJSON(1.0), name="w1")
        second = tk.getHermesWorkflowFromJSON(_workflowJSON(9.0), name="w2")
        with pytest.raises(AttributeError, match="'int' object has no attribute 'replace'"):
            tk.compareWorkflowObj([first, second], longFormat=True)


@pytest.mark.unit
class TestCompareWorkflows:
    @pytest.fixture()
    def twoFiles(self, tmp_path):
        first = tmp_path / "a.json"
        first.write_text(json.dumps(_workflowJSON(1.0)))
        second = tmp_path / "b.json"
        second.write_text(json.dumps(_workflowJSON(9.0)))
        return str(first), str(second)

    def test_two_workflow_files_are_compared_column_per_file(self, tk, twoFiles):
        table = tk.compareWorkflows(list(twoFiles))
        assert list(table.columns) == list(twoFiles)
        assert table.loc["Parameters_alpha"].tolist() == [1.0, 9.0]

    def test_transposing_puts_the_workflows_in_the_rows(self, tk, twoFiles):
        table = tk.compareWorkflows(list(twoFiles), transpose=True)
        assert list(table.index) == list(twoFiles)
        assert list(table.columns) == ["Parameters_alpha"]

    def test_an_identifier_that_is_neither_in_the_db_nor_a_file_leaves_nothing_to_compare(
        self, tk
    ):
        """The failure is logged and the name skipped, so the comparison ends
        up with no workflows at all and pandas refuses to concatenate."""
        with pytest.raises(ValueError, match="No objects to concatenate"):
            tk.compareWorkflows(["noSuchWorkflow"])

    @pytest.mark.xfail(
        strict=True,
        reason="B140: every DB-backed comparison rebuilds the hermes objects "
               "with workflow(..., resource=simulationDoc['resource']), but "
               "hermes.workflow.__init__ takes Resource_path, not resource -- so "
               "compareWorkflows, compareWorkflowInGroup and workflowTable all "
               "raise TypeError: unexpected keyword argument 'resource' as soon "
               "as a matching document exists. "
               "See the consolidated findings issue.",
    )
    def test_a_group_name_should_compare_the_workflows_stored_in_that_group(
        self, tk, twoWorkflowsInGroup
    ):
        table = tk.compareWorkflows("flow")
        assert sorted(table.columns) == ["flow_0000", "flow_0001"]

    def test_a_group_name_currently_raises_a_typeerror(self, tk, twoWorkflowsInGroup):
        """Characterisation of B140."""
        with pytest.raises(TypeError, match="unexpected keyword argument 'resource'"):
            tk.compareWorkflows("flow")


@pytest.mark.unit
class TestCompareWorkflowInGroupAndWorkflowTable:
    def test_an_empty_group_has_nothing_to_compare(self, tk):
        assert tk.compareWorkflowInGroup("noSuchGroup") is None

    def test_an_empty_group_gives_the_same_answer_through_workflowtable(self, tk):
        assert tk.workflowTable("noSuchGroup") is None

    @pytest.mark.xfail(
        strict=True,
        reason="B140: compareWorkflowInGroup rebuilds each stored workflow "
               "with workflow(..., resource=...) while hermes.workflow.__init__ "
               "expects Resource_path, so any non-empty group raises TypeError. "
               "See the consolidated findings issue.",
    )
    def test_a_populated_group_should_be_comparable(self, tk, twoWorkflowsInGroup):
        table = tk.compareWorkflowInGroup("flow")
        assert sorted(table.index) == ["flow_0000", "flow_0001"]

    def test_a_populated_group_currently_raises_a_typeerror(self, tk, twoWorkflowsInGroup):
        """Characterisation of B140."""
        with pytest.raises(TypeError, match="unexpected keyword argument 'resource'"):
            tk.compareWorkflowInGroup("flow")

    def test_workflowtable_inherits_the_same_failure(self, tk, twoWorkflowsInGroup):
        """Characterisation of B140: workflowTable is a thin delegation to
        compareWorkflowInGroup, so it cannot behave any better."""
        with pytest.raises(TypeError, match="unexpected keyword argument 'resource'"):
            tk.workflowTable("flow")


# ---------------------------------------------------------------------------
# Listing
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestListWorkflows:
    def test_by_default_only_the_workflow_names_are_listed(self, tk, twoWorkflowsInGroup):
        listing = tk.listWorkflows("flow")
        assert sorted(entry["workflowName"] for entry in listing) == ["flow_0000", "flow_0001"]
        assert all(list(entry) == ["workflowName"] for entry in listing)

    def test_an_unknown_group_lists_nothing(self, tk, twoWorkflowsInGroup):
        assert tk.listWorkflows("noSuchGroup") == []

    def test_listing_nodes_adds_the_node_list_of_each_workflow(self, tk, twoWorkflowsInGroup):
        listing = tk.listWorkflows("flow", listNodes=True)
        assert all(entry["nodes"] == _REQUIRED_NODES for entry in listing)

    def test_listing_parameters_adds_the_stored_parameters(self, tk, twoWorkflowsInGroup):
        listing = tk.listWorkflows("flow", listParameters=True)
        byName = {entry["workflowName"]: entry["parameters"] for entry in listing}
        assert byName["flow_0000"]["Parameters"] == {"alpha": 1.0}
        assert byName["flow_0001"]["Parameters"] == {"alpha": 9.0}


@pytest.mark.unit
class TestListGroups:
    def test_the_solver_the_group_and_every_workflow_are_printed(
        self, tk, twoWorkflowsInGroup, capsys
    ):
        tk.listGroups()
        printed = capsys.readouterr().out
        for expected in ["simpleFoam", "* flow", "+ flow_0000", "+ flow_0001"]:
            assert expected in printed

    def test_it_returns_nothing_it_only_prints(self, tk, twoWorkflowsInGroup):
        assert tk.listGroups() is None

    def test_the_workflow_names_can_be_suppressed(self, tk, twoWorkflowsInGroup, capsys):
        tk.listGroups(workflowName=False)
        printed = capsys.readouterr().out
        assert "* flow" in printed
        assert "flow_0000" not in printed

    def test_filtering_by_another_solver_prints_nothing(self, tk, twoWorkflowsInGroup, capsys):
        tk.listGroups(solver="buoyantReactingFoam")
        assert capsys.readouterr().out == ""

    def test_a_project_with_no_workflows_prints_nothing(self, tk, capsys):
        tk.listGroups()
        assert capsys.readouterr().out == ""


# ---------------------------------------------------------------------------
# Deletion and naming
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestDeleteWorkflowInGroup:
    def test_every_workflow_of_the_group_is_removed(self, tk, twoWorkflowsInGroup):
        tk.deleteWorkflowInGroup("flow")
        assert len(tk.getWorkflowDocumentsInGroup("flow")) == 0

    def test_other_groups_are_untouched(self, tk):
        tk.addWorkflowToGroup(_workflowJSON(1.0), "flowA")
        tk.addWorkflowToGroup(_workflowJSON(2.0), "flowB")
        tk.deleteWorkflowInGroup("flowA")
        assert len(tk.getWorkflowDocumentsInGroup("flowB")) == 1

    def test_excluded_workflows_survive_the_deletion(self, tk, twoWorkflowsInGroup):
        tk.deleteWorkflowInGroup("flow", exclude=["flow_0000"])
        remaining = tk.getWorkflowDocumentsInGroup("flow")
        assert [doc.desc["workflowName"] for doc in remaining] == ["flow_0000"]

    def test_a_shallow_delete_leaves_the_resource_on_disk(self, tk, tmp_path):
        caseDirectory = tmp_path / "case"
        caseDirectory.mkdir()
        tk.addWorkflowToGroup(_workflowJSON(), "flow", resource=str(caseDirectory))
        tk.deleteWorkflowInGroup("flow")
        assert caseDirectory.exists()

    def test_a_deep_delete_removes_the_resource_directory_too(self, tk, tmp_path):
        caseDirectory = tmp_path / "case"
        caseDirectory.mkdir()
        (caseDirectory / "inside.txt").write_text("x")
        tk.addWorkflowToGroup(_workflowJSON(), "flow", resource=str(caseDirectory))
        tk.deleteWorkflowInGroup("flow", deepDelete=True)
        assert not caseDirectory.exists()

    def test_resetting_the_counter_restarts_the_group_numbering_at_one(
        self, tk, twoWorkflowsInGroup
    ):
        """As the docstring says: the reset takes the counter back to 1, not
        to 0 -- the zeroth name is only ever handed out to a brand new group."""
        tk.deleteWorkflowInGroup("flow", resetCounter=True)
        doc = tk.addWorkflowToGroup(_workflowJSON(1.0), "flow")
        assert doc.desc["workflowName"] == "flow_0001"

    def test_keeping_the_counter_continues_the_group_numbering(self, tk, twoWorkflowsInGroup):
        tk.deleteWorkflowInGroup("flow", resetCounter=False)
        doc = tk.addWorkflowToGroup(_workflowJSON(1.0), "flow")
        assert doc.desc["workflowName"] == "flow_0002"


@pytest.mark.unit
class TestFindAvailableName:
    def test_a_brand_new_group_starts_at_zero(self, tk):
        newID, name = tk.findAvailableName("brandNew")
        assert (newID, name) == (0, "brandNew_0000")

    def test_successive_calls_hand_out_successive_names(self, tk):
        tk.findAvailableName("brandNew")
        assert tk.findAvailableName("brandNew") == (1, "brandNew_0001")

    def test_the_name_is_built_the_same_way_as_getworkflowname(self, tk):
        newID, name = tk.findAvailableName("brandNew")
        assert name == tk.getworkFlowName("brandNew", newID)

    @pytest.mark.xfail(
        strict=True,
        reason="B141: findAvailableName draws from the counter "
               "'simulations_<group>', while addWorkflowToGroup draws from the "
               "counter '<group>'. The two never agree, so the 'available' name "
               "it returns can be -- and for a group populated the normal way "
               "always is -- the name of a workflow that already exists. "
               "See the consolidated findings issue.",
    )
    def test_the_available_name_should_not_already_be_taken(self, tk, twoWorkflowsInGroup):
        _, name = tk.findAvailableName("flow")
        assert tk.getWorkflowDocumentByName(name) is None

    def test_the_available_name_currently_collides_with_an_existing_workflow(
        self, tk, twoWorkflowsInGroup
    ):
        """Characterisation of B141."""
        newID, name = tk.findAvailableName("flow")
        assert (newID, name) == (0, "flow_0000")
        assert tk.getWorkflowDocumentByName("flow_0000") is not None
