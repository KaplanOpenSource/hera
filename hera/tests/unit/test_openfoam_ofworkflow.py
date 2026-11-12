"""OFWorkflow: the hermes-workflow specializations for OpenFOAM.

Covers every accessor and mutator of ``abstractWorkflow`` (the group/groupID
properties and their setters, the node accessors ``controlDict`` /
``fvSolution`` / ``fvScheme`` / ``fileWriter`` / ``parameters`` /
``buildAllRun`` / ``defineNewBoundaryConditions``, and ``__setitem__``), of
``workflow_Eulerian`` (its required-node validation, ``blockMesh`` /
``snappyHexMesh``, the IC/BC dispatcher and the two blockMesh geometry
setters), of ``workflow_Lagrangian`` and
``workflow_StochasticLagrangianSolver`` (``dispersionName``,
``originalFlowFieldName``, ``dispersionFlowFieldName`` and its setter,
``dispersionDuration``), and of the five thin solver subclasses
(``workflow_simpleFoam``, ``workflow_buoyantReactingFoam``,
``workflow_scalarTransportFoam``, ``workflow_indoorFOAMBoussinesq``,
``workflow_homogenousWindLogProfile``).

The workflows here are real ``hermes.workflow`` objects built from a minimal
but genuine workflow JSON: nodes named exactly as
``abstractWorkflow._requiredNodeList`` demands (``controlDict``,
``fvSolution``, ``fvSchemes``, ``fileWriter``, ``Parameters``), each a
``general.Parameters`` stub, plus whichever optional nodes a given test
needs.  Nothing is stubbed -- the real hermes package is installed, and the
constructor really does build the task network.  The hera-document-backed
tests use a document produced by the real workflow toolkit against the
in-memory datalayer, so ``document.save()`` in the setters is exercised for
real.

Bugs pinned here: B145 (``controlDict`` reads the wrong node name), B146
(``fvScheme`` reads the wrong node name), B147 (``__setitem__`` calls a
misspelt, name-mangled super method), B148 and B149 (``set_blockMesh_blockHeight``
moves the bottom face instead of the top one, and then derives a zero cell
count from the vertex it has just overwritten) and B150 (``workflowGroupID``'s
setter rejects the integer type its own getter yields).

Deliberately not covered here:

* ``foam_snappyhexmesh_addobject`` beyond the fact that it is an
  unimplemented stub -- its body is a docstring and nothing else, so there
  is no behaviour to assert.
* Anything that actually runs OpenFOAM or luigi.  ``build()``, ``write()``
  and the task-network internals belong to hermes, not to hera, and are
  outside this module's surface.
* ``workflowDescirption`` beyond its default -- it is a bare class
  attribute that no code in the repository ever assigns.
"""
import pytest

from hera.simulations.openFoam.OFWorkflow import (
    abstractWorkflow,
    workflow_buoyantReactingFoam,
    workflow_Eulerian,
    workflow_homogenousWindLogProfile,
    workflow_indoorFOAMBoussinesq,
    workflow_Lagrangian,
    workflow_scalarTransportFoam,
    workflow_simpleFoam,
    workflow_StochasticLagrangianSolver,
)

# Exactly the nodes abstractWorkflow.__init__ records as required, and so
# exactly the nodes workflow_Eulerian.__init__ insists on finding.
_REQUIRED_NODES = ["controlDict", "fvSolution", "fvSchemes", "fileWriter", "Parameters"]


def _workflowJSON(extraNodes=(), parameters=None, solver="simpleFoam"):
    """A minimal but real hermes workflow JSON.

    Every node is a ``general.Parameters`` stub; only ``Parameters`` may
    carry values.  A fresh dict is built on every call because the hermes
    constructor mutates the JSON it is handed (it appends ``finalnode_xx``).
    """
    nodeNames = list(_REQUIRED_NODES) + [name for name in extraNodes]
    nodes = {
        nodeName: {"type": "general.Parameters", "Execution": {"input_parameters": {}}}
        for nodeName in nodeNames
    }
    if parameters is not None:
        nodes["Parameters"]["Execution"]["input_parameters"] = dict(parameters)
    return {"workflow": {"nodeList": nodeNames, "nodes": nodes, "solver": solver}}


@pytest.fixture()
def workingDirectory(tmp_path):
    """A throwaway working directory, so no Resource_path points at the repo."""
    return str(tmp_path)


@pytest.fixture()
def eulerian(workingDirectory):
    """A valid Eulerian workflow carrying blockMesh and buildAllRun as well."""
    return workflow_Eulerian(
        _workflowJSON(extraNodes=["blockMesh", "buildAllRun", "defineNewBoundaryConditions"]),
        name="eulerian",
        WD_path=workingDirectory,
    )


@pytest.fixture()
def workflowDocument(unit_toolkit_factory):
    """A real hera workflow document, with a groupName and a groupID in desc."""
    from hera import toolkitHome

    toolkit = unit_toolkit_factory(toolkitHome.SIMULATIONS_WORKFLOWS)
    document = toolkit.addWorkflowToGroup(_workflowJSON(parameters={"alpha": 1.0}), "flow")
    return toolkit, document


# ---------------------------------------------------------------------------
# abstractWorkflow: construction
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestAbstractWorkflowConstruction:
    def test_a_workflow_json_becomes_a_real_hermes_workflow(self, workingDirectory):
        workflowObject = abstractWorkflow(
            _workflowJSON(), name="abstract", WD_path=workingDirectory
        )
        import hermes

        assert isinstance(workflowObject, hermes.workflow)
        assert workflowObject.name == "abstract"

    def test_the_required_node_list_is_recorded_on_the_instance(self, workingDirectory):
        workflowObject = abstractWorkflow(_workflowJSON(), WD_path=workingDirectory)
        assert workflowObject._requiredNodeList == _REQUIRED_NODES

    def test_the_required_node_list_is_not_shared_between_instances(self, workingDirectory):
        """It is declared as a class attribute but assigned per instance, so
        appending to one workflow's list must not touch another's."""
        first = abstractWorkflow(_workflowJSON(), WD_path=workingDirectory)
        second = abstractWorkflow(_workflowJSON(), WD_path=workingDirectory)
        first._requiredNodeList.append("blockMesh")
        assert second._requiredNodeList == _REQUIRED_NODES
        assert abstractWorkflow._requiredNodeList is None

    def test_the_hera_document_defaults_to_none(self, workingDirectory):
        workflowObject = abstractWorkflow(_workflowJSON(), WD_path=workingDirectory)
        assert workflowObject.workflowHeraDocument is None

    def test_the_hera_document_is_kept_as_handed_over(
        self, workingDirectory, workflowDocument
    ):
        _, document = workflowDocument
        workflowObject = abstractWorkflow(
            _workflowJSON(), workflowHeraDocument=document, WD_path=workingDirectory
        )
        assert workflowObject.workflowHeraDocument is document

    def test_the_workflow_description_defaults_to_none(self, workingDirectory):
        workflowObject = abstractWorkflow(_workflowJSON(), WD_path=workingDirectory)
        assert workflowObject.workflowDescirption is None

    def test_the_abstract_class_does_not_validate_the_required_nodes(
        self, workingDirectory
    ):
        """Only workflow_Eulerian enforces the list; the abstract base merely
        records it."""
        incomplete = _workflowJSON()
        del incomplete["workflow"]["nodes"]["controlDict"]
        incomplete["workflow"]["nodeList"].remove("controlDict")
        assert abstractWorkflow(incomplete, WD_path=workingDirectory) is not None


# ---------------------------------------------------------------------------
# abstractWorkflow: the node accessors
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestNodeAccessors:
    def test_fvsolution_reaches_the_fvsolution_node(self, eulerian):
        assert eulerian.fvSolution.name == "fvSolution"

    def test_filewriter_reaches_the_filewriter_node(self, eulerian):
        assert eulerian.fileWriter.name == "fileWriter"

    def test_parameters_reaches_the_parameters_node(self, workingDirectory):
        workflowObject = abstractWorkflow(
            _workflowJSON(parameters={"alpha": 3.0}), WD_path=workingDirectory
        )
        assert workflowObject.parameters.name == "Parameters"
        assert workflowObject.parameters.parameters == {"alpha": 3.0}

    def test_buildallrun_reaches_the_buildallrun_node_when_it_exists(self, eulerian):
        assert eulerian.buildAllRun.name == "buildAllRun"

    def test_definenewboundaryconditions_reaches_its_node_when_it_exists(self, eulerian):
        assert eulerian.defineNewBoundaryConditions.name == "defineNewBoundaryConditions"

    def test_an_optional_node_that_is_absent_gives_none_rather_than_raising(
        self, workingDirectory
    ):
        """hermes' __getitem__ answers None for an unknown node, so the
        optional-node accessors are safe on a bare workflow."""
        workflowObject = abstractWorkflow(_workflowJSON(), WD_path=workingDirectory)
        assert workflowObject.buildAllRun is None
        assert workflowObject.defineNewBoundaryConditions is None

    def test_an_accessor_returns_a_live_view_of_the_node(self, eulerian):
        """The node object writes straight through to the workflow JSON, so a
        value set through the accessor is visible in the JSON itself."""
        eulerian.fvSolution["solverTolerance"] = 1e-6
        stored = eulerian.nodes["fvSolution"]["Execution"]["input_parameters"]
        assert stored == {"solverTolerance": 1e-6}


@pytest.mark.unit
class TestControlDictAccessor:
    """B145: the accessor and the required-node list disagree on the casing."""

    def test_the_required_node_is_spelt_with_a_lowercase_c(self, eulerian):
        """Baseline: the node really is there under the name the class
        validates against, so the failure below is about the accessor's key,
        not about a missing node."""
        assert "controlDict" in eulerian.nodes
        assert eulerian["controlDict"].name == "controlDict"

    @pytest.mark.xfail(
        strict=True,
        reason="B145: abstractWorkflow.controlDict returns self['ControlDict'] "
               "(capital C), but _requiredNodeList -- which workflow_Eulerian "
               "enforces at construction -- names the node 'controlDict', as do "
               "all the flow templates under hera/doc. hermes' __getitem__ "
               "answers None for an unknown node instead of raising, so the "
               "accessor silently yields None on every workflow that has just "
               "passed the class's own validation. "
               "See the consolidated findings issue.",
    )
    def test_the_control_dict_node_should_be_reachable_through_the_accessor(self, eulerian):
        assert eulerian.controlDict.name == "controlDict"

    def test_the_control_dict_accessor_is_currently_always_none(self, eulerian):
        """Characterisation of B145."""
        assert eulerian.controlDict is None

    def test_the_accessor_only_answers_for_a_node_named_with_a_capital_c(
        self, workingDirectory
    ):
        """Characterisation of B145: the accessor is not broken outright, it
        reads a name that only the Lagrangian dispersion templates use."""
        workflowObject = abstractWorkflow(
            _workflowJSON(extraNodes=["ControlDict"]), WD_path=workingDirectory
        )
        assert workflowObject.controlDict.name == "ControlDict"


@pytest.mark.unit
class TestFvSchemeAccessor:
    """B146: the accessor reads a node name that exists nowhere."""

    def test_the_required_node_is_spelt_in_the_plural(self, eulerian):
        """Baseline: the node is present under the validated name."""
        assert eulerian["fvSchemes"].name == "fvSchemes"

    @pytest.mark.xfail(
        strict=True,
        reason="B146: abstractWorkflow.fvScheme returns self['fvScheme'], "
               "singular, while the node is named 'fvSchemes' in "
               "_requiredNodeList and in every workflow template in the "
               "repository -- no node named 'fvScheme' exists anywhere. Since "
               "hermes' __getitem__ answers None for an unknown node, the "
               "accessor can never return anything but None. "
               "See the consolidated findings issue.",
    )
    def test_the_fvschemes_node_should_be_reachable_through_the_accessor(self, eulerian):
        assert eulerian.fvScheme.name == "fvSchemes"

    def test_the_fvscheme_accessor_is_currently_always_none(self, eulerian):
        """Characterisation of B146."""
        assert eulerian.fvScheme is None


# ---------------------------------------------------------------------------
# abstractWorkflow: workflowGroup / workflowGroupID
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestWorkflowGroupWithoutADocument:
    def test_the_group_name_of_a_document_less_workflow_is_none(self, workingDirectory):
        workflowObject = abstractWorkflow(_workflowJSON(), WD_path=workingDirectory)
        assert workflowObject.workflowGroup is None

    def test_the_group_id_of_a_document_less_workflow_is_none(self, workingDirectory):
        workflowObject = abstractWorkflow(_workflowJSON(), WD_path=workingDirectory)
        assert workflowObject.workflowGroupID is None

    def test_setting_the_group_name_without_a_document_is_a_silent_no_op(
        self, workingDirectory
    ):
        """The type check runs first, but the write is guarded by the document,
        so a document-less workflow accepts the assignment and forgets it."""
        workflowObject = abstractWorkflow(_workflowJSON(), WD_path=workingDirectory)
        workflowObject.workflowGroup = "someGroup"
        assert workflowObject.workflowGroup is None

    def test_setting_the_group_id_without_a_document_is_a_silent_no_op(
        self, workingDirectory
    ):
        workflowObject = abstractWorkflow(_workflowJSON(), WD_path=workingDirectory)
        workflowObject.workflowGroupID = "someID"
        assert workflowObject.workflowGroupID is None

    def test_a_non_string_group_name_is_rejected_even_without_a_document(
        self, workingDirectory
    ):
        """The validation precedes the document guard, so it always applies."""
        workflowObject = abstractWorkflow(_workflowJSON(), WD_path=workingDirectory)
        with pytest.raises(ValueError, match="Group name must be a string"):
            workflowObject.workflowGroup = 17


@pytest.mark.unit
class TestWorkflowGroupWithADocument:
    @pytest.fixture()
    def documentBacked(self, workflowDocument, workingDirectory):
        toolkit, document = workflowDocument
        workflowObject = workflow_simpleFoam(
            _workflowJSON(), workflowHeraDocument=document, name="flow_0000",
            WD_path=workingDirectory,
        )
        return toolkit, document, workflowObject

    def test_the_group_name_comes_from_the_documents_desc(self, documentBacked):
        _, document, workflowObject = documentBacked
        assert workflowObject.workflowGroup == document.desc["groupName"] == "flow"

    def test_the_group_id_comes_from_the_documents_desc(self, documentBacked):
        _, document, workflowObject = documentBacked
        assert workflowObject.workflowGroupID == document.desc["groupID"] == 0

    def test_a_new_group_name_round_trips_through_the_getter(self, documentBacked):
        _, _, workflowObject = documentBacked
        workflowObject.workflowGroup = "renamed"
        assert workflowObject.workflowGroup == "renamed"

    def test_a_new_group_name_is_persisted_to_the_database(self, documentBacked):
        toolkit, document, workflowObject = documentBacked
        workflowObject.workflowGroup = "renamed"
        reread = toolkit.getWorkflowDocumentByName(document.desc["workflowName"])
        assert reread.desc["groupName"] == "renamed"

    def test_a_new_group_id_round_trips_through_the_getter(self, documentBacked):
        _, _, workflowObject = documentBacked
        workflowObject.workflowGroupID = "0007"
        assert workflowObject.workflowGroupID == "0007"

    def test_a_new_group_id_is_persisted_to_the_database(self, documentBacked):
        toolkit, document, workflowObject = documentBacked
        workflowObject.workflowGroupID = "0007"
        reread = toolkit.getWorkflowDocumentByName(document.desc["workflowName"])
        assert reread.desc["groupID"] == "0007"

    def test_the_group_name_setter_leaves_the_group_id_alone(self, documentBacked):
        _, _, workflowObject = documentBacked
        workflowObject.workflowGroup = "renamed"
        assert workflowObject.workflowGroupID == 0

    def test_a_non_string_group_name_is_refused_and_nothing_is_written(
        self, documentBacked
    ):
        _, _, workflowObject = documentBacked
        with pytest.raises(ValueError, match="Group name must be a string"):
            workflowObject.workflowGroup = ["flow"]
        assert workflowObject.workflowGroup == "flow"

    @pytest.mark.xfail(
        strict=True,
        reason="B150: the workflowGroupID getter yields whatever "
               "addWorkflowToGroup stored, which is an int (the group counter), "
               "but its setter demands isinstance(value,str) and raises "
               "ValueError('Group name must be a string' -- the group-name "
               "message, copy-pasted) for anything else. So the property cannot "
               "round-trip its own value: `wf.workflowGroupID = "
               "wf.workflowGroupID` raises, and writing the id back needs a "
               "string, which then makes the stored type disagree with every "
               "other document in the group. "
               "See the consolidated findings issue.",
    )
    def test_the_group_id_should_accept_the_type_its_getter_yields(self, documentBacked):
        _, _, workflowObject = documentBacked
        workflowObject.workflowGroupID = workflowObject.workflowGroupID
        assert workflowObject.workflowGroupID == 0

    def test_the_group_id_setter_currently_refuses_an_integer(self, documentBacked):
        """Characterisation of B150."""
        _, _, workflowObject = documentBacked
        assert isinstance(workflowObject.workflowGroupID, int)
        with pytest.raises(ValueError, match="Group name must be a string"):
            workflowObject.workflowGroupID = 1


# ---------------------------------------------------------------------------
# abstractWorkflow.__setitem__
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestSetItem:
    """B147: adding a node through the mapping interface cannot work."""

    @staticmethod
    def _newNodeValue():
        return {
            "buildAllRun": {"name": "myCommand", "couldRunInParallel": False},
            "fileWriter": {"myFile": "content"},
            "node": {"type": "general.Parameters", "Execution": {"input_parameters": {}}},
        }

    @pytest.mark.xfail(
        strict=True,
        reason="B147: abstractWorkflow.__setitem__ ends with "
               "super().__setitem(key=key,value=value) -- the trailing "
               "'__setitem' is a typo for '__setitem__', and because it starts "
               "with two underscores inside a class body it is name-mangled to "
               "_abstractWorkflow__setitem, which exists on nothing. Every node "
               "assignment therefore raises AttributeError: 'super' object has "
               "no attribute '_abstractWorkflow__setitem'. "
               "See the consolidated findings issue.",
    )
    def test_assigning_a_node_should_add_it_to_the_workflow(self, eulerian):
        eulerian["myNewNode"] = self._newNodeValue()
        assert "myNewNode" in eulerian.nodes

    def test_assigning_a_node_currently_raises_an_attributeerror(self, eulerian):
        """Characterisation of B147."""
        with pytest.raises(AttributeError, match="_abstractWorkflow__setitem"):
            eulerian["myNewNode"] = self._newNodeValue()

    def test_the_failed_assignment_has_already_mutated_buildallrun_and_filewriter(
        self, eulerian
    ):
        """Characterisation of B147: the two side effects happen before the
        failing call, so a caught AttributeError leaves the workflow with an
        execution entry and a file-writer entry for a node that was never
        added."""
        with pytest.raises(AttributeError):
            eulerian["myNewNode"] = self._newNodeValue()

        assert eulerian.buildAllRun["myNewNode"] == {
            "name": "myCommand", "couldRunInParallel": False,
        }
        assert eulerian.fileWriter["myNewNode"] == {"myFile": "content"}
        assert "myNewNode" not in eulerian.nodes

    def test_a_value_without_a_buildallrun_entry_fails_on_that_key_first(self, eulerian):
        with pytest.raises(KeyError, match="buildAllRun"):
            eulerian["myNewNode"] = {"fileWriter": {}}

    def test_a_workflow_without_a_buildallrun_node_cannot_be_assigned_to_at_all(
        self, workingDirectory
    ):
        """buildAllRun is None for a workflow that lacks the node, and None
        does not support item assignment."""
        workflowObject = abstractWorkflow(_workflowJSON(), WD_path=workingDirectory)
        with pytest.raises(TypeError, match="does not support item assignment"):
            workflowObject["myNewNode"] = self._newNodeValue()


# ---------------------------------------------------------------------------
# workflow_Eulerian: validation and node accessors
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestEulerianValidation:
    @pytest.mark.parametrize("missingNode", _REQUIRED_NODES)
    def test_every_required_node_is_insisted_upon_by_name(
        self, workingDirectory, missingNode
    ):
        workflowBody = _workflowJSON()
        del workflowBody["workflow"]["nodes"][missingNode]
        workflowBody["workflow"]["nodeList"].remove(missingNode)
        with pytest.raises(ValueError, match=f"The node {missingNode} does not exist"):
            workflow_Eulerian(workflowBody, WD_path=workingDirectory)

    def test_the_rejection_says_it_is_not_a_flow_workflow(self, workingDirectory):
        workflowBody = _workflowJSON()
        del workflowBody["workflow"]["nodes"]["fvSolution"]
        with pytest.raises(ValueError, match="Not a flow workflow"):
            workflow_Eulerian(workflowBody, WD_path=workingDirectory)

    def test_the_validation_reads_the_nodes_and_not_the_nodelist(self, workingDirectory):
        """A node present in 'nodes' but missing from 'nodeList' still counts as
        existing, because the check is `node not in self.workflowJSON['nodes']`."""
        workflowBody = _workflowJSON()
        workflowBody["workflow"]["nodeList"].remove("fvSolution")
        assert workflow_Eulerian(workflowBody, WD_path=workingDirectory) is not None

    def test_a_complete_workflow_is_accepted(self, eulerian):
        assert isinstance(eulerian, workflow_Eulerian)
        assert isinstance(eulerian, abstractWorkflow)


@pytest.mark.unit
class TestEulerianMeshNodeAccessors:
    def test_blockmesh_reaches_the_blockmesh_node(self, eulerian):
        assert eulerian.blockMesh.name == "blockMesh"

    def test_a_missing_blockmesh_node_gives_none(self, workingDirectory):
        """blockMesh is not in the required list, so a valid Eulerian workflow
        may have none -- and then the accessor answers None."""
        workflowObject = workflow_Eulerian(_workflowJSON(), WD_path=workingDirectory)
        assert workflowObject.blockMesh is None

    def test_snappyhexmesh_is_none_when_the_node_is_absent(self, eulerian):
        assert eulerian.snappyHexMesh is None

    def test_snappyhexmesh_reaches_the_node_when_it_is_present(self, workingDirectory):
        workflowObject = workflow_Eulerian(
            _workflowJSON(extraNodes=["snappyHexMesh"]), WD_path=workingDirectory
        )
        assert workflowObject.snappyHexMesh.name == "snappyHexMesh"


# ---------------------------------------------------------------------------
# workflow_Eulerian: the IC/BC dispatcher
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestICandBC:
    def test_the_node_handler_writes_the_fields_into_the_boundary_node(self, eulerian):
        eulerian.ICandBCHandler_node({"data": {"U": [1, 0, 0]}})
        assert eulerian.defineNewBoundaryConditions["fields"] == {"U": [1, 0, 0]}

    def test_the_dispatcher_routes_by_the_type_field(self, eulerian):
        eulerian.prepareICandBC({"ICandBC": [{"type": "node", "data": {"p": 0}}]})
        assert eulerian.defineNewBoundaryConditions["fields"] == {"p": 0}

    def test_the_dispatcher_applies_every_entry_in_order(self, eulerian):
        """All the entries here target the same key, so the last one wins --
        which is exactly how the sequencing shows up."""
        eulerian.prepareICandBC(
            {"ICandBC": [{"type": "node", "data": {"p": 0}},
                         {"type": "node", "data": {"p": 9}}]}
        )
        assert eulerian.defineNewBoundaryConditions["fields"] == {"p": 9}

    def test_an_empty_icandbc_list_does_nothing(self, eulerian):
        eulerian.prepareICandBC({"ICandBC": []})
        assert list(eulerian.defineNewBoundaryConditions.keys()) == []

    def test_the_dispatcher_returns_nothing(self, eulerian):
        assert eulerian.prepareICandBC({"ICandBC": [{"type": "node", "data": {}}]}) is None

    def test_an_unknown_handler_type_surfaces_as_a_missing_attribute(self, eulerian):
        """The dispatch is a getattr on a name built from the type, so an
        unsupported type is reported as a missing method rather than as a
        validation error."""
        with pytest.raises(AttributeError, match="ICandBCHandler_noSuchType"):
            eulerian.prepareICandBC({"ICandBC": [{"type": "noSuchType"}]})

    def test_the_node_handler_needs_the_boundary_node_to_exist(self, workingDirectory):
        workflowObject = workflow_Eulerian(_workflowJSON(), WD_path=workingDirectory)
        with pytest.raises(TypeError, match="does not support item assignment"):
            workflowObject.ICandBCHandler_node({"data": {}})


# ---------------------------------------------------------------------------
# workflow_Eulerian: the blockMesh geometry setters
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestSetBlockMeshBlockBoundaries:
    @pytest.fixture()
    def bounded(self, eulerian):
        eulerian.set_blockMesh_blockBoundaries(
            minx=0, maxx=10, miny=0, maxy=20, minz=0, maxz=30, dx=1, dy=2, dz=3
        )
        return eulerian

    def test_eight_vertices_are_written(self, bounded):
        assert len(bounded.blockMesh["vertices"]) == 8

    def test_the_bottom_face_is_the_four_corners_at_the_minimum_height(self, bounded):
        """OpenFOAM's blockMesh convention: vertices 0-3 are the base, walked
        counter-clockwise."""
        assert bounded.blockMesh["vertices"][:4] == [
            [0, 0, 0], [10, 0, 0], [10, 20, 0], [0, 20, 0],
        ]

    def test_the_top_face_repeats_those_corners_at_the_maximum_height(self, bounded):
        assert bounded.blockMesh["vertices"][4:] == [
            [0, 0, 30], [10, 0, 30], [10, 20, 30], [0, 20, 30],
        ]

    def test_the_hex_is_the_eight_vertex_indices_in_order(self, bounded):
        assert bounded.blockMesh["hex"] == [0, 1, 2, 3, 4, 5, 6, 7]

    def test_the_cell_count_is_the_span_divided_by_the_spacing(self, bounded):
        assert bounded.blockMesh["cellCount"] == [10, 10, 10]

    def test_the_cell_counts_are_integers(self, bounded):
        assert all(isinstance(count, int) for count in bounded.blockMesh["cellCount"])

    def test_a_span_that_does_not_divide_evenly_is_rounded_up(self, eulerian):
        """Ceiling division, so the block is always fully covered."""
        eulerian.set_blockMesh_blockBoundaries(
            minx=0, maxx=10, miny=0, maxy=10, minz=0, maxz=10, dx=3, dy=3, dz=3
        )
        assert eulerian.blockMesh["cellCount"] == [4, 4, 4]

    def test_it_returns_nothing(self, eulerian):
        result = eulerian.set_blockMesh_blockBoundaries(0, 1, 0, 1, 0, 1, 1, 1, 1)
        assert result is None

    def test_it_needs_a_blockmesh_node(self, workingDirectory):
        workflowObject = workflow_Eulerian(_workflowJSON(), WD_path=workingDirectory)
        with pytest.raises(TypeError, match="does not support item assignment"):
            workflowObject.set_blockMesh_blockBoundaries(0, 1, 0, 1, 0, 1, 1, 1, 1)


@pytest.mark.unit
class TestSetBlockMeshBlockHeight:
    """B148 and B149: raising the block moves the wrong face, and then computes
    the cell count from the vertex it has just overwritten."""

    @pytest.fixture()
    def raised(self, eulerian):
        eulerian.set_blockMesh_blockBoundaries(
            minx=0, maxx=10, miny=0, maxy=20, minz=0, maxz=30, dx=1, dy=2, dz=3
        )
        eulerian.set_blockMesh_blockHeight(Z=100, dz=5)
        return eulerian

    @pytest.mark.xfail(
        strict=True,
        reason="B148: set_blockMesh_blockHeight is meant to lift the top face of "
               "the block to Z, but it writes "
               "`for i in range(len(verticsList[4:])): verticsList[i][2] = Z`. "
               "The slice is only used for its length (4), while the assignment "
               "indexes i itself -- so vertices 0-3, the *bottom* face, are "
               "moved to Z and the top face keeps its old height. The block ends "
               "up inverted. The index should be i+4. "
               "See the consolidated findings issue.",
    )
    def test_raising_the_block_should_move_the_top_face(self, raised):
        assert [vertex[2] for vertex in raised.blockMesh["vertices"]] == [
            0, 0, 0, 0, 100, 100, 100, 100,
        ]

    def test_raising_the_block_currently_moves_the_bottom_face(self, raised):
        """Characterisation of B148."""
        assert [vertex[2] for vertex in raised.blockMesh["vertices"]] == [
            100, 100, 100, 100, 30, 30, 30, 30,
        ]

    def test_the_x_and_y_coordinates_are_left_untouched(self, raised):
        assert [vertex[:2] for vertex in raised.blockMesh["vertices"]] == [
            [0, 0], [10, 0], [10, 20], [0, 20],
            [0, 0], [10, 0], [10, 20], [0, 20],
        ]

    @pytest.mark.xfail(
        strict=True,
        reason="B149: after moving the vertices, set_blockMesh_blockHeight reads "
               "minZ = verticsList[0][2] -- but vertex 0 is one of the four it "
               "has just overwritten with Z (see B148), so minZ == Z and "
               "cellCount[2] = (Z-minZ)/dz is always exactly 0. A zero cell "
               "count in z makes the mesh unbuildable. It is also a float, "
               "unlike the ints set_blockMesh_blockBoundaries writes. "
               "See the consolidated findings issue.",
    )
    def test_raising_the_block_should_recount_the_cells_in_z(self, raised):
        assert raised.blockMesh["cellCount"][2] == 20

    def test_the_z_cell_count_currently_collapses_to_zero(self, raised):
        """Characterisation of B149."""
        assert raised.blockMesh["cellCount"][2] == 0.0

    def test_the_x_and_y_cell_counts_are_left_untouched(self, raised):
        assert raised.blockMesh["cellCount"][:2] == [10, 10]

    def test_it_returns_nothing(self, eulerian):
        eulerian.set_blockMesh_blockBoundaries(0, 1, 0, 1, 0, 1, 1, 1, 1)
        assert eulerian.set_blockMesh_blockHeight(Z=5, dz=1) is None

    def test_it_needs_the_boundaries_to_have_been_set_first(self, eulerian):
        """It reads blockMesh['vertices'], which only exists once
        set_blockMesh_blockBoundaries has written it."""
        with pytest.raises(KeyError, match="vertices"):
            eulerian.set_blockMesh_blockHeight(Z=5, dz=1)


@pytest.mark.unit
class TestFoamSnappyHexMeshAddObject:
    def test_it_is_an_unimplemented_stub(self, eulerian):
        """Its body is a docstring and nothing else, so it accepts anything and
        answers None without touching the workflow."""
        assert eulerian.foam_snappyhexmesh_addobject("someData", "someFile.obj") is None
        assert eulerian.snappyHexMesh is None


# ---------------------------------------------------------------------------
# workflow_Lagrangian
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestLagrangianWorkflow:
    def test_it_is_an_abstract_workflow_and_not_an_eulerian_one(self, workingDirectory):
        workflowObject = workflow_Lagrangian(
            _workflowJSON(), name="lagrangian", WD_path=workingDirectory
        )
        assert isinstance(workflowObject, abstractWorkflow)
        assert not isinstance(workflowObject, workflow_Eulerian)

    def test_it_does_not_impose_the_eulerian_required_nodes(self, workingDirectory):
        """The list is still recorded, but only workflow_Eulerian enforces it,
        so a Lagrangian workflow without a controlDict node is accepted."""
        workflowBody = _workflowJSON()
        del workflowBody["workflow"]["nodes"]["controlDict"]
        workflowBody["workflow"]["nodeList"].remove("controlDict")

        workflowObject = workflow_Lagrangian(workflowBody, WD_path=workingDirectory)
        assert workflowObject._requiredNodeList == _REQUIRED_NODES
        assert "controlDict" not in workflowObject.nodes

    def test_the_hera_document_is_forwarded_to_the_base_class(
        self, workingDirectory, workflowDocument
    ):
        _, document = workflowDocument
        workflowObject = workflow_Lagrangian(
            _workflowJSON(), workflowHeraDocument=document, WD_path=workingDirectory
        )
        assert workflowObject.workflowGroup == "flow"


# ---------------------------------------------------------------------------
# workflow_StochasticLagrangianSolver
# ---------------------------------------------------------------------------

_DISPERSION_PARAMETERS = {
    "originalFlowField": "myFlow",
    "dispersionFlowField": "myDispersion",
    "dispersionDuration": 3600,
}


@pytest.mark.unit
class TestStochasticLagrangianSolver:
    @pytest.fixture()
    def dispersion(self, workingDirectory):
        return workflow_StochasticLagrangianSolver(
            _workflowJSON(parameters=_DISPERSION_PARAMETERS),
            name="myDispersionRun",
            WD_path=workingDirectory,
        )

    def test_a_missing_dispersion_flow_field_is_refused_at_construction(
        self, workingDirectory
    ):
        parameters = dict(_DISPERSION_PARAMETERS)
        del parameters["dispersionFlowField"]
        with pytest.raises(ValueError, match="must have a dispersionFlowField"):
            workflow_StochasticLagrangianSolver(
                _workflowJSON(parameters=parameters), WD_path=workingDirectory
            )

    def test_it_is_a_lagrangian_workflow(self, dispersion):
        assert isinstance(dispersion, workflow_Lagrangian)

    def test_the_dispersion_name_is_the_workflow_name(self, dispersion):
        assert dispersion.dispersionName == "myDispersionRun"

    def test_an_unnamed_workflow_has_no_dispersion_name(self, workingDirectory):
        workflowObject = workflow_StochasticLagrangianSolver(
            _workflowJSON(parameters=_DISPERSION_PARAMETERS), WD_path=workingDirectory
        )
        assert workflowObject.dispersionName is None

    def test_the_original_flow_field_name_comes_from_the_parameters(self, dispersion):
        assert dispersion.originalFlowFieldName == "myFlow"

    def test_the_dispersion_flow_field_name_joins_both_fields_with_dff(self, dispersion):
        assert dispersion.dispersionFlowFieldName == "myFlow_DFF_myDispersion"

    def test_the_dispersion_duration_comes_from_the_parameters(self, dispersion):
        assert dispersion.dispersionDuration == 3600

    def test_a_missing_duration_is_only_noticed_when_it_is_read(self, workingDirectory):
        """Only dispersionFlowField is validated at construction."""
        workflowObject = workflow_StochasticLagrangianSolver(
            _workflowJSON(parameters={"dispersionFlowField": "d"}),
            WD_path=workingDirectory,
        )
        with pytest.raises(KeyError, match="dispersionDuration"):
            workflowObject.dispersionDuration

    def test_setting_the_dispersion_flow_field_name_writes_the_parameter(self, dispersion):
        dispersion.dispersionFlowFieldName = "otherDispersion"
        assert dispersion.parameters.parameters["dispersionFlowField"] == "otherDispersion"

    def test_the_setter_coerces_its_value_to_a_string(self, dispersion):
        dispersion.dispersionFlowFieldName = 17
        assert dispersion.parameters.parameters["dispersionFlowField"] == "17"

    def test_the_setter_writes_the_component_and_not_the_composite_name(self, dispersion):
        """Asymmetric by design: the setter names the dispersion half, while the
        getter prefixes the original flow field and the _DFF_ marker, so the
        round-trip is deliberately not the identity."""
        dispersion.dispersionFlowFieldName = "otherDispersion"
        assert dispersion.dispersionFlowFieldName == "myFlow_DFF_otherDispersion"

    def test_the_setter_leaves_the_original_flow_field_alone(self, dispersion):
        dispersion.dispersionFlowFieldName = "otherDispersion"
        assert dispersion.originalFlowFieldName == "myFlow"

    def test_the_composite_name_follows_a_change_of_original_flow_field(self, dispersion):
        dispersion.parameters["originalFlowField"] = "anotherFlow"
        assert dispersion.dispersionFlowFieldName == "anotherFlow_DFF_myDispersion"


# ---------------------------------------------------------------------------
# The thin solver subclasses
# ---------------------------------------------------------------------------

_SOLVER_CLASSES = [
    workflow_simpleFoam,
    workflow_buoyantReactingFoam,
    workflow_scalarTransportFoam,
    workflow_indoorFOAMBoussinesq,
    workflow_homogenousWindLogProfile,
]


@pytest.mark.unit
@pytest.mark.parametrize(
    "solverClass", _SOLVER_CLASSES, ids=[cls.__name__ for cls in _SOLVER_CLASSES]
)
class TestSolverSubclasses:
    """All five add nothing but a name -- their __init__ is a bare super() call,
    so each must behave exactly like workflow_Eulerian."""

    def test_it_is_an_eulerian_workflow(self, solverClass, workingDirectory):
        workflowObject = solverClass(_workflowJSON(), WD_path=workingDirectory)
        assert isinstance(workflowObject, workflow_Eulerian)

    def test_it_inherits_the_eulerian_required_node_validation(
        self, solverClass, workingDirectory
    ):
        workflowBody = _workflowJSON()
        del workflowBody["workflow"]["nodes"]["fileWriter"]
        with pytest.raises(ValueError, match="Not a flow workflow"):
            solverClass(workflowBody, WD_path=workingDirectory)

    def test_the_name_and_the_hera_document_reach_the_base_class(
        self, solverClass, workingDirectory, workflowDocument
    ):
        _, document = workflowDocument
        workflowObject = solverClass(
            _workflowJSON(), workflowHeraDocument=document, name="named",
            WD_path=workingDirectory,
        )
        assert workflowObject.name == "named"
        assert workflowObject.workflowGroup == "flow"

    def test_it_exposes_the_node_accessors(self, solverClass, workingDirectory):
        workflowObject = solverClass(
            _workflowJSON(extraNodes=["blockMesh"]), WD_path=workingDirectory
        )
        assert workflowObject.fvSolution.name == "fvSolution"
        assert workflowObject.blockMesh.name == "blockMesh"
