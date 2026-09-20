"""openFoam/postProcess/VTKPipeline.py: the pipeline/filter tree, its JSON
form, the DB-backed cache bookkeeping of ``registeredVTKPipeLine``, and the
argument forwarding into ParaView.

What is covered *deeply* (real logic, no ParaView involved at all)
-----------------------------------------------------------------
``VTKPipeLine`` (construction, ``newVTKPipelineFilter``, ``addFilter``,
``addFilterFromObj``, ``addExistingFilter``, ``__setitem__``,
``__getitem__``, ``registerPipeline``, ``toJSON``, ``allFilterNames``),
the whole of ``VTKFilter`` (including ``set_param`` / ``enforce_param`` and
the ``fullName`` property, which is broken), every ``vtkFilter_*``
subclass and its setters -- all of these only build Python dicts, lists
and strings, so the assertions below are on real values.

On ``registeredVTKPipeLine``: ``__init__`` (both the DB branch and the
directory branch, plus the rejection), ``_resolveRequestedFilters``,
``_parseTimeList``, ``_filterCachedTimesteps``, ``_updateCacheDB``,
``_buildFilterQuery``, ``getFilterOutputFileExt``,
``getFilterOutputFilePath`` and ``clearCache`` are all driven against the
real in-memory datalayer (mongomock), so the cache documents, counters,
queries and file deletions asserted here are genuine.

What is covered only *shallowly*, and why
-----------------------------------------
``paraview`` is not installed and cannot be (it ships as a binary, not on
PyPI).  ``hera/tests/unit/_stubs.py`` registers ``paraview.simple`` as an
empty namespace module, so *nothing* in it exists until a test puts it
there.  That is actually the useful seam: the tests below install
recording stand-ins (``_RecordingProxy``) and assert on the recorded call
log -- which filter type was instantiated, with which ``Input`` and
``guiName``, which properties were set in which order, and that
``UpdatePipeline`` was called.  Those calls are the entire observable
behaviour of ``_buildFilterLayer``.

Consequently:

* ``_buildFilterLayer`` -- covered for structure (names, nesting order,
  parameter application order, nested ``a.b`` property traversal).  The
  geometric *effect* of those calls is unobservable here and is not
  asserted.
* ``_buildAndExecuteParaViewPipeline`` -- argument forwarding only: which
  reader is created, whether ``CellArrays`` is narrowed, exactly which
  ``filtersDict`` reaches ``writeCase``, and that every ParaView source is
  deleted afterwards.  ``paraviewOpenFOAM.initializeReader`` /
  ``writeCase`` are patched on the class, because the real ones dereference
  a ``pvsimple`` name that does not even exist in ``pvOpenFOAMBase`` when
  ``vtkmodules`` is absent.
* ``getData`` -- the orchestration (time resolution, cache hit/miss, DB
  persistence, reload) is real; only step 4 (the ParaView execution) is
  stubbed, with the stub writing a real parquet file so that the final
  ``document.getData()`` reload is exercised for real.  The *contents* of
  the returned frame are whatever the stub wrote; no physics is asserted.
* ``getRegularData`` / ``getNonRegularData`` -- pure forwarding, tested as
  forwarding.
* ``registeredVTKPipeLine.__init__`` with ``serverName`` set is not
  covered: it reaches ``pvsimple.Connect`` inside ``pvOpenFOAMBase``, whose
  ``pvsimple`` import is skipped entirely when ``vtkmodules`` is missing,
  so the only thing observable would be a stub artefact (``NameError``).

Deliberately not duplicated: ``OFToolkit.getVTKPipelineCacheDocuments`` /
``clearVTKPipelineCache`` live in ``toolkit.py`` and are covered by
``test_openfoam_toolkit_vtk_cache.py`` (which pins B101).
"""
import copy
import os

import pandas
import pytest

from hera.simulations.openFoam import (
    CASETYPE_DECOMPOSED,
    CASETYPE_RECONSTRUCTED,
    TYPE_VTK_FILTER,
)
from hera.simulations.openFoam.postProcess import VTKPipeline as pipelineModule
from hera.simulations.openFoam.postProcess.VTKPipeline import (
    VTKFilter,
    VTKPipeLine,
    registeredVTKPipeLine,
    vtkFilter_CellCenters,
    vtkFilter_DescriptiveStatistics,
    vtkFilter_ExtractBlock,
    vtkFilter_IntegrateVariables,
    vtkFilter_PlotOverLine,
    vtkFilter_Slice,
)
from hera.simulations.openFoam.postProcess.pvOpenFOAMBase import paraviewOpenFOAM
from hera.utils import dictToMongoQuery

ALL_FILTER_TYPES = [
    "CellCenters",
    "DescriptiveStatistics",
    "ExtractBlock",
    "IntegrateVariables",
    "PlotOverLine",
    "Slice",
]


# ---------------------------------------------------------------------------
# Recording stand-ins for the (absent) ParaView proxies
# ---------------------------------------------------------------------------

class _RecordingProxy:
    """A stand-in for a ParaView proxy whose only behaviour is its log.

    Property *writes* are recorded and thrown away; property *reads* always
    return a nested sub-proxy.  That mirrors the one piece of real ParaView
    behaviour the pipeline code depends on -- assigning ``SliceType =
    "Plane"`` swaps in a Plane sub-proxy, so the following
    ``SliceType.Origin`` assignment lands on that sub-proxy -- without
    pretending to implement anything else.
    """

    def __init__(self, log, path):
        d = self.__dict__
        d["_log"] = log
        d["_path"] = path
        d["_children"] = {}

    def __setattr__(self, name, value):
        self.__dict__["_log"].append(("set", self.__dict__["_path"], name, value))

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        children = self.__dict__["_children"]
        if name not in children:
            children[name] = _RecordingProxy(
                self.__dict__["_log"], f"{self.__dict__['_path']}.{name}"
            )
        return children[name]

    def UpdatePipeline(self):
        self.__dict__["_log"].append(("update", self.__dict__["_path"]))

    def __repr__(self):
        return f"<proxy {self.__dict__['_path']}>"


@pytest.fixture()
def paraviewLog(monkeypatch):
    """Install recording filter constructors on the stubbed paraview.simple.

    Returns the shared call log.  Entries are one of:
        ("create", filterType, guiName, inputProxy)
        ("set", proxyPath, propertyName, value)
        ("update", proxyPath)
        ("delete", proxyName)
    """
    log = []

    def _makeConstructor(filterType):
        def _construct(**kwargs):
            proxy = _RecordingProxy(log, kwargs.get("guiName"))
            log.append(("create", filterType, kwargs.get("guiName"), kwargs.get("Input")))
            return proxy

        return _construct

    for filterType in ALL_FILTER_TYPES:
        monkeypatch.setattr(
            pipelineModule.pvsimple, filterType, _makeConstructor(filterType), raising=False
        )

    monkeypatch.setattr(pipelineModule.pvsimple, "GetSources", lambda: {}, raising=False)
    monkeypatch.setattr(
        pipelineModule.pvsimple,
        "Delete",
        lambda proxy: log.append(("delete", proxy)),
        raising=False,
    )
    return log


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def of(unit_toolkit_factory):
    from hera import toolkitHome

    return unit_toolkit_factory(toolkitHome.SIMULATIONS_OPENFOAM)


@pytest.fixture()
def caseDirectory(tmp_path):
    """A minimal on-disk case directory named ``<group>_<index>``."""
    directory = tmp_path / "flow_0001"
    (directory / "0").mkdir(parents=True)
    (directory / "10").mkdir()
    return directory


@pytest.fixture()
def pipeline(of):
    """A pipeline with one write-enabled Slice under a non-writing ExtractBlock."""
    pl = VTKPipeLine(of)
    block = pl.addFilter("ExtractBlock", VTKPipeLine.FILTER_EXTRACTBLOCK, write=False)
    sliceFilter = block.addFilter("slice", VTKPipeLine.FILTER_SLICE, write=True)
    sliceFilter.setPlaneOrigin([0, 0, 10])
    sliceFilter.setPlaneNormal([0, 0, 1])
    return pl


@pytest.fixture()
def registered(pipeline, caseDirectory):
    """``pipeline`` bound to ``caseDirectory`` through the directory branch."""
    return pipeline.registerPipeline(str(caseDirectory))


# ---------------------------------------------------------------------------
# VTKPipeLine: construction and the filter factory
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestVTKPipeLineConstruction:
    def test_a_new_pipeline_starts_with_no_filters(self, of):
        assert VTKPipeLine(of).filters == {}

    def test_it_adopts_an_existing_filter_dict_by_reference(self, of):
        existing = {"a": "placeholder"}
        assert VTKPipeLine(of, vtkPipeline=existing).filters is existing

    def test_two_pipelines_do_not_share_the_default_filter_dict(self, of):
        first = VTKPipeLine(of)
        second = VTKPipeLine(of)
        first.addFilter("s", VTKPipeLine.FILTER_SLICE)
        assert second.filters == {}

    def test_it_keeps_the_datalayer_it_was_given(self, of):
        assert VTKPipeLine(of).datalayer is of

    def test_the_filter_name_constants_are_all_real_filter_types(self):
        constants = [
            VTKPipeLine.FILTER_CELLCENTERS,
            VTKPipeLine.FILTER_SLICE,
            VTKPipeLine.FILTER_PLOTOVERLINE,
            VTKPipeLine.FILTER_EXTRACTBLOCK,
            VTKPipeLine.FILTER_INTEGRATEVARIABLES,
        ]
        assert all(name in ALL_FILTER_TYPES for name in constants)


@pytest.mark.unit
class TestNewVTKPipelineFilter:
    @pytest.mark.parametrize("filterType", ALL_FILTER_TYPES)
    def test_it_builds_the_class_matching_each_known_filter_type(self, filterType):
        built = VTKPipeLine.newVTKPipelineFilter(name="n", filterType=filterType)
        assert type(built).__name__ == f"vtkFilter_{filterType}"
        assert built.filterType == filterType
        assert built.name == "n"

    def test_it_rejects_an_unknown_filter_type_and_lists_the_known_ones(self):
        with pytest.raises(ValueError) as err:
            VTKPipeLine.newVTKPipelineFilter(name="n", filterType="NoSuchFilter")
        message = str(err.value)
        assert "NoSuchFilter is not Known" in message
        assert all(known in message for known in ALL_FILTER_TYPES)

    def test_the_known_types_are_discovered_from_the_modules_vtkFilter_classes(self):
        discovered = sorted(
            name.split("_")[1] for name in dir(pipelineModule) if name.startswith("vtkFilter")
        )
        assert discovered == ALL_FILTER_TYPES

    def test_it_defaults_the_parameter_list_to_empty(self):
        assert VTKPipeLine.newVTKPipelineFilter(name="n", filterType="CellCenters").params == []

    def test_two_filters_built_with_the_default_do_not_share_a_parameter_list(self):
        first = VTKPipeLine.newVTKPipelineFilter(name="a", filterType="Slice")
        second = VTKPipeLine.newVTKPipelineFilter(name="b", filterType="Slice")
        first.set_param("Q", 1)
        assert second.params == []

    def test_it_forwards_the_write_flag(self):
        assert VTKPipeLine.newVTKPipelineFilter(name="n", filterType="Slice", write=False).write is False

    def test_it_forwards_an_explicit_parameter_list(self):
        params = [("Resolution", 7)]
        built = VTKPipeLine.newVTKPipelineFilter(
            name="n", filterType="PlotOverLine", params=params
        )
        assert built.params == [("Resolution", 7)]

    def test_a_type_that_passes_the_name_check_but_cannot_be_located_is_refused(self, monkeypatch):
        """The second guard: a name in dir(module) that pydoc cannot resolve.

        Unreachable through the public API -- the name list is derived from
        the very module pydoc is asked to look in -- so the only way to
        exercise the guard is to make the lookup fail.
        """
        monkeypatch.setattr(pipelineModule.pydoc, "locate", lambda path: None)
        with pytest.raises(RuntimeError, match="Slice does not exist"):
            VTKPipeLine.newVTKPipelineFilter(name="n", filterType="Slice")


@pytest.mark.unit
class TestVTKPipeLineAddFilter:
    def test_add_filter_registers_the_new_filter_under_its_name(self, of):
        pl = VTKPipeLine(of)
        added = pl.addFilter("mySlice", VTKPipeLine.FILTER_SLICE)
        assert pl.filters == {"mySlice": added}
        assert pl["mySlice"] is added

    def test_add_filter_returns_the_filter_so_it_can_be_configured(self, of):
        pl = VTKPipeLine(of)
        assert isinstance(pl.addFilter("s", VTKPipeLine.FILTER_SLICE), vtkFilter_Slice)

    def test_add_filter_forwards_an_explicit_parameter_list(self, of):
        pl = VTKPipeLine(of)
        added = pl.addFilter("s", VTKPipeLine.FILTER_SLICE, params=[("SliceType", "Plane")])
        assert added.params == [("SliceType", "Plane")]

    def test_add_filter_from_obj_keys_the_filter_by_its_own_name(self, of):
        pl = VTKPipeLine(of)
        existing = vtkFilter_CellCenters(name="cc", write=True)
        pl.addFilterFromObj(existing)
        assert pl["cc"] is existing

    def test_the_deprecated_add_existing_filter_still_adds_the_filter(self, of):
        pl = VTKPipeLine(of)
        existing = vtkFilter_CellCenters(name="cc", write=True)
        pl.addExistingFilter(existing)
        assert pl["cc"] is existing

    def test_setitem_can_key_a_filter_under_a_name_of_its_own(self, of):
        pl = VTKPipeLine(of)
        existing = vtkFilter_CellCenters(name="realName", write=True)
        pl["otherKey"] = existing
        assert pl["otherKey"] is existing

    def test_a_second_filter_with_the_same_name_replaces_the_first(self, of):
        pl = VTKPipeLine(of)
        pl.addFilter("s", VTKPipeLine.FILTER_SLICE)
        second = pl.addFilter("s", VTKPipeLine.FILTER_CELLCENTERS)
        assert pl.filters == {"s": second}


@pytest.mark.unit
class TestVTKPipeLineLookup:
    def test_it_resolves_a_dotted_path_down_the_downstream_tree(self, pipeline):
        assert pipeline["ExtractBlock.slice"] is pipeline["ExtractBlock"]["slice"]

    def test_an_unknown_root_filter_raises_a_named_key_error(self, pipeline):
        with pytest.raises(KeyError, match="The filter nope is not found in the current pipeline"):
            pipeline["nope"]

    def test_an_unknown_leaf_of_a_known_path_raises_a_named_key_error(self, pipeline):
        with pytest.raises(KeyError, match="not found"):
            pipeline["ExtractBlock.nope"]

    def test_a_non_string_key_is_not_treated_as_a_path(self, pipeline):
        with pytest.raises(AttributeError):
            pipeline[17]


@pytest.mark.unit
class TestVTKPipeLineToJSON:
    def test_the_json_is_wrapped_in_a_filters_key(self, pipeline):
        assert list(pipeline.toJSON().keys()) == ["filters"]

    def test_it_records_type_write_params_and_downstream_for_each_filter(self, pipeline):
        block = pipeline.toJSON()["filters"]["ExtractBlock"]
        assert block["filterType"] == "ExtractBlock"
        assert block["write"] is False
        assert block["params"] == [("Selectors", [])]
        assert list(block["downstream"]) == ["slice"]

    def test_nested_filters_appear_under_their_fathers_downstream(self, pipeline):
        nested = pipeline.toJSON()["filters"]["ExtractBlock"]["downstream"]["slice"]
        assert nested["filterType"] == "Slice"
        assert nested["write"] is True
        assert nested["params"] == [
            ("SliceType", "Plane"),
            ("SliceType.Origin", [0, 0, 10]),
            ("SliceType.Normal", [0, 0, 1]),
        ]

    def test_an_empty_pipeline_serialises_to_an_empty_filter_map(self, of):
        assert VTKPipeLine(of).toJSON() == dict(filters={})

    def test_the_json_is_keyed_by_the_filter_name_not_the_dict_key(self, of):
        """Characterisation: toJSON ignores the key ``__setitem__`` was given."""
        pl = VTKPipeLine(of)
        pl["otherKey"] = vtkFilter_CellCenters(name="realName", write=True)
        assert list(pl.toJSON()["filters"]) == ["realName"]


@pytest.mark.unit
class TestAllFilterNames:
    def test_it_returns_dotted_paths_depth_first(self, pipeline):
        assert pipeline.allFilterNames() == ["ExtractBlock", "ExtractBlock.slice"]

    def test_write_only_drops_the_filters_not_marked_for_output(self, pipeline):
        assert pipeline.allFilterNames(writeOnly=True) == ["ExtractBlock.slice"]

    def test_an_empty_pipeline_has_no_filter_names(self, of):
        assert VTKPipeLine(of).allFilterNames() == []

    def test_siblings_are_returned_in_insertion_order(self, of):
        pl = VTKPipeLine(of)
        pl.addFilter("first", VTKPipeLine.FILTER_CELLCENTERS)
        pl.addFilter("second", VTKPipeLine.FILTER_SLICE)
        assert pl.allFilterNames() == ["first", "second"]

    def test_it_descends_more_than_one_level(self, of):
        pl = VTKPipeLine(of)
        block = pl.addFilter("a", VTKPipeLine.FILTER_EXTRACTBLOCK)
        middle = block.addFilter("b", VTKPipeLine.FILTER_EXTRACTBLOCK)
        middle.addFilter("c", VTKPipeLine.FILTER_CELLCENTERS)
        assert pl.allFilterNames() == ["a", "a.b", "a.b.c"]


# ---------------------------------------------------------------------------
# VTKFilter
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestVTKFilter:
    def test_a_new_filter_has_an_empty_downstream_of_its_own(self):
        first = VTKFilter(name="a", filterType="Slice", write=True, params=[])
        second = VTKFilter(name="b", filterType="Slice", write=True, params=[])
        first["x"] = "placeholder"
        assert second.downstream == {}

    def test_json_coerces_a_truthy_write_flag_to_a_boolean(self):
        assert VTKFilter("a", "Slice", "yes", []).toJSON()["a"]["write"] is True

    def test_json_coerces_a_missing_write_flag_to_false(self):
        assert VTKFilter("a", "Slice", None, []).toJSON()["a"]["write"] is False

    def test_json_of_a_leaf_has_an_empty_downstream_map(self):
        assert VTKFilter("a", "Slice", True, []).toJSON()["a"]["downstream"] == {}

    def test_add_filter_defaults_write_to_none_which_serialises_as_false(self):
        father = VTKFilter("a", "ExtractBlock", True, [])
        child = father.addFilter("b", "CellCenters")
        assert child.write is None
        assert father.toJSON()["a"]["downstream"]["b"]["write"] is False

    def test_add_filter_forwards_an_explicit_parameter_list_to_the_child(self):
        father = VTKFilter("a", "ExtractBlock", True, [])
        child = father.addFilter("b", "Slice", params=[("SliceType", "Plane")])
        assert child.params == [("SliceType", "Plane")]

    def test_add_filter_from_obj_keys_the_child_by_its_own_name(self):
        father = VTKFilter("a", "ExtractBlock", True, [])
        child = vtkFilter_CellCenters(name="cc", write=True)
        father.addFilterFromObj(child)
        assert father.downstream == {"cc": child}

    def test_it_resolves_a_dotted_downstream_path(self):
        father = VTKFilter("a", "ExtractBlock", True, [])
        middle = father.addFilter("b", "ExtractBlock")
        leaf = middle.addFilter("c", "CellCenters")
        assert father["b.c"] is leaf

    def test_an_unknown_downstream_name_raises_a_named_key_error(self):
        father = VTKFilter("a", "ExtractBlock", True, [])
        with pytest.raises(KeyError, match="Filter nope not found"):
            father["nope"]

    def test_setitem_adds_a_downstream_filter_under_an_arbitrary_key(self):
        father = VTKFilter("a", "ExtractBlock", True, [])
        father["k"] = "placeholder"
        assert father.downstream == {"k": "placeholder"}


@pytest.mark.unit
class TestSetParam:
    def test_it_appends_a_parameter_that_was_not_set_before(self):
        f = VTKFilter("a", "Slice", True, [])
        f.set_param("Resolution", 10)
        assert f.params == [("Resolution", 10)]

    def test_it_replaces_a_parameter_in_place_keeping_the_order(self):
        f = VTKFilter("a", "Slice", True, [("first", 1), ("second", 2)])
        f.set_param("first", 99)
        assert f.params == [("first", 99), ("second", 2)]

    def test_setting_the_same_value_twice_does_not_duplicate_it(self):
        f = VTKFilter("a", "Slice", True, [])
        f.set_param("Resolution", 10)
        f.set_param("Resolution", 10)
        assert f.params == [("Resolution", 10)]

    def test_it_rebinds_the_list_rather_than_mutating_the_one_it_was_given(self):
        given = [("first", 1)]
        f = VTKFilter("a", "Slice", True, given)
        f.set_param("first", 99)
        assert given == [("first", 1)]
        assert f.params == [("first", 99)]

    def test_a_parameter_stored_as_a_list_pair_is_still_replaced_by_name(self):
        """Characterisation: only ``param[0]`` is inspected, so list pairs match."""
        f = VTKFilter("a", "Slice", True, [["first", 1]])
        f.set_param("first", 99)
        assert f.params == [("first", 99)]


@pytest.mark.unit
class TestEnforceParam:
    def test_it_sets_the_parameter_when_nothing_set_it_yet(self):
        f = VTKFilter("a", "Slice", True, [])
        f.enforce_param("SliceType", "Plane")
        assert f.params == [("SliceType", "Plane")]

    def test_it_accepts_the_parameter_already_being_the_enforced_value(self):
        f = VTKFilter("a", "Slice", True, [("SliceType", "Plane")])
        f.enforce_param("SliceType", "Plane")
        assert f.params == [("SliceType", "Plane")]

    def test_it_refuses_to_overwrite_a_conflicting_value(self):
        f = VTKFilter("a", "Slice", True, [("SliceType", "Sphere")])
        with pytest.raises(RuntimeError, match="SliceType must not be different than Plane"):
            f.enforce_param("SliceType", "Plane")

    def test_the_conflicting_value_is_left_untouched_after_the_refusal(self):
        f = VTKFilter("a", "Slice", True, [("SliceType", "Sphere")])
        with pytest.raises(RuntimeError):
            f.enforce_param("SliceType", "Plane")
        assert f.params == [("SliceType", "Sphere")]

    def test_an_already_enforced_value_stored_as_a_list_pair_is_rejected(self):
        """Characterisation: the membership test compares a tuple against the
        stored pairs, so a pair stored as a list (the shape a JSON round-trip
        produces) is reported as a conflict even when the value is identical.
        """
        f = VTKFilter("a", "Slice", True, [["SliceType", "Plane"]])
        with pytest.raises(RuntimeError, match="must not be different"):
            f.enforce_param("SliceType", "Plane")


@pytest.mark.unit
class TestVTKFilterFullName:
    @pytest.mark.xfail(
        strict=True,
        reason="B316: VTKFilter.fullName binds `fatherName` only inside "
               "`if filter.father is not None`, then returns it "
               "unconditionally.  `father` is a class attribute that is None "
               "and is never assigned anywhere in the module (neither "
               "__init__, addFilter, addFilterFromObj nor __setitem__ set "
               "it), so the guard is always false and the property always "
               "raises UnboundLocalError.  The dead branch is broken too: it "
               "calls `filter.father.fullName()`, invoking the property's "
               "string result as if it were a method.  See the consolidated "
               "findings issue.",
    )
    def test_a_root_filter_reports_its_own_name_as_its_full_name(self):
        assert VTKFilter("a", "Slice", True, []).fullName == "a"

    def test_full_name_currently_raises_for_every_filter(self):
        """Characterisation of B316."""
        f = VTKFilter("a", "Slice", True, [])
        assert f.father is None
        with pytest.raises(UnboundLocalError):
            f.fullName

    def test_the_dead_father_branch_calls_the_property_as_a_method(self, monkeypatch):
        """Characterisation of B316: even a filter that *does* have a father
        cannot report a full name, because the branch invokes
        ``father.fullName()`` -- the property already evaluated to a string.
        """

        class _StubFather:
            fullName = "root"

        child = VTKFilter("child", "Slice", True, [])
        monkeypatch.setattr(child, "father", _StubFather())
        with pytest.raises(TypeError, match="not callable"):
            child.fullName

    def test_no_code_path_in_the_module_ever_assigns_a_father(self):
        """Characterisation of B316: the guard cannot become true."""
        father = VTKFilter("a", "ExtractBlock", True, [])
        child = father.addFilter("b", "CellCenters")
        father.addFilterFromObj(vtkFilter_CellCenters(name="c", write=True))
        assert child.father is None
        assert father.downstream["c"].father is None


# ---------------------------------------------------------------------------
# The concrete filter subclasses
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestSliceFilter:
    def test_it_declares_the_slice_filter_type(self):
        assert vtkFilter_Slice(name="s", write=True).filterType == "Slice"

    def test_setting_the_origin_enforces_a_plane_slice_first(self):
        s = vtkFilter_Slice(name="s", write=True)
        s.setPlaneOrigin([1, 2, 3])
        assert s.params == [("SliceType", "Plane"), ("SliceType.Origin", [1, 2, 3])]

    def test_setting_the_normal_enforces_a_plane_slice_first(self):
        s = vtkFilter_Slice(name="s", write=True)
        s.setPlaneNormal([0, 0, 1])
        assert s.params == [("SliceType", "Plane"), ("SliceType.Normal", [0, 0, 1])]

    def test_origin_and_normal_together_enforce_the_plane_only_once(self):
        s = vtkFilter_Slice(name="s", write=True)
        s.setPlaneOrigin([1, 2, 3])
        s.setPlaneNormal([0, 0, 1])
        assert [name for name, _ in s.params].count("SliceType") == 1

    def test_the_origin_can_be_changed_without_reordering_the_parameters(self):
        s = vtkFilter_Slice(name="s", write=True)
        s.setPlaneOrigin([1, 2, 3])
        s.setPlaneNormal([0, 0, 1])
        s.setPlaneOrigin([9, 9, 9])
        assert s.params == [
            ("SliceType", "Plane"),
            ("SliceType.Origin", [9, 9, 9]),
            ("SliceType.Normal", [0, 0, 1]),
        ]

    def test_a_conflicting_slice_type_makes_the_origin_setter_refuse(self):
        s = vtkFilter_Slice(name="s", write=True, params=[("SliceType", "Sphere")])
        with pytest.raises(RuntimeError, match="SliceType"):
            s.setPlaneOrigin([1, 2, 3])


@pytest.mark.unit
class TestPlotOverLineFilter:
    def test_it_declares_the_plot_over_line_filter_type(self):
        assert vtkFilter_PlotOverLine(name="p", write=True).filterType == "PlotOverLine"

    def test_the_sampling_pattern_constants_are_the_paraview_strings(self):
        assert vtkFilter_PlotOverLine.SAMPLE_UNIFORMLY == "Sample Uniformly"
        assert vtkFilter_PlotOverLine.SAMPLE_AT_CELL_BOUNDARIES == "Sample At Cell Boundaries"
        assert vtkFilter_PlotOverLine.SAMPLE_AT_SEGMENT_CENTERS == "Sample At Segment Centers"

    def test_setting_the_sample_pattern_records_it_as_a_parameter(self):
        p = vtkFilter_PlotOverLine(name="p", write=True)
        p.setSamplePattern(vtkFilter_PlotOverLine.SAMPLE_AT_CELL_BOUNDARIES)
        assert p.params == [("SamplingPattern", "Sample At Cell Boundaries")]

    def test_a_uniform_resolution_enforces_uniform_sampling_first(self):
        p = vtkFilter_PlotOverLine(name="p", write=True)
        p.setUniformSampleResolution(120)
        assert p.params == [("SamplingPattern", "Sample Uniformly"), ("Resolution", 120)]

    def test_a_uniform_resolution_is_refused_under_another_sampling_pattern(self):
        p = vtkFilter_PlotOverLine(name="p", write=True)
        p.setSamplePattern(vtkFilter_PlotOverLine.SAMPLE_AT_SEGMENT_CENTERS)
        with pytest.raises(RuntimeError, match="SamplingPattern"):
            p.setUniformSampleResolution(120)

    def test_the_end_points_are_recorded_in_order(self):
        p = vtkFilter_PlotOverLine(name="p", write=True)
        p.setPoints([0, 0, 0], [1, 1, 1])
        assert p.params == [("Point1", [0, 0, 0]), ("Point2", [1, 1, 1])]


@pytest.mark.unit
class TestCellCentersFilter:
    def test_it_declares_the_cell_centers_filter_type_and_no_parameters(self):
        cc = vtkFilter_CellCenters(name="cc", write=True)
        assert (cc.filterType, cc.params) == ("CellCenters", [])


@pytest.mark.unit
class TestExtractBlockFilter:
    def test_an_empty_patch_list_still_produces_a_selectors_parameter(self):
        eb = vtkFilter_ExtractBlock(name="eb", write=True)
        assert eb.params == [("Selectors", [])]

    def test_each_patch_becomes_a_root_boundary_selector(self):
        eb = vtkFilter_ExtractBlock(name="eb", write=True, patchList=["ground", "top"])
        assert eb.params == [("Selectors", ["/Root/boundary/ground", "/Root/boundary/top"])]

    def test_the_internal_mesh_selector_is_appended_last_when_requested(self):
        eb = vtkFilter_ExtractBlock(name="eb", write=True, patchList=["ground"], internalMesh=True)
        assert eb.params == [("Selectors", ["/Root/boundary/ground", "/Root/internalMesh"])]

    def test_the_selectors_parameter_is_placed_before_the_caller_supplied_ones(self):
        eb = vtkFilter_ExtractBlock(name="eb", write=True, params=[("Other", 1)])
        assert eb.params == [("Selectors", []), ("Other", 1)]

    def test_two_filters_built_from_the_defaults_do_not_share_a_patch_list(self):
        first = vtkFilter_ExtractBlock(name="a", write=True)
        second = vtkFilter_ExtractBlock(name="b", write=True)
        first.setRegionsToExtract(patchList=["ground"], internalMesh=False)
        assert second.params == [("Selectors", [])]

    def test_resetting_the_regions_replaces_the_selectors_in_place(self):
        eb = vtkFilter_ExtractBlock(name="eb", write=True, patchList=["ground"])
        eb.setRegionsToExtract(patchList=["top"], internalMesh=False)
        assert eb.params == [("Selectors", ["/Root/boundary/top"])]

    def test_the_setter_includes_the_internal_mesh_by_default(self):
        """Characterisation: the setter defaults internalMesh to True while the
        constructor defaults it to False."""
        eb = vtkFilter_ExtractBlock(name="eb", write=True)
        eb.setRegionsToExtract()
        assert eb.params == [("Selectors", ["/Root/internalMesh"])]


@pytest.mark.unit
class TestDescriptiveStatisticsFilter:
    def test_it_declares_the_descriptive_statistics_filter_type(self):
        ds = vtkFilter_DescriptiveStatistics(name="ds", write=True)
        assert (ds.filterType, ds.params) == ("DescriptiveStatistics", [])

    def test_the_variables_of_interest_are_recorded_as_one_parameter(self):
        ds = vtkFilter_DescriptiveStatistics(name="ds", write=True)
        ds.setVariablesOfInterest(["U_x", "T"])
        assert ds.params == [("VariablesofInterest", ["U_x", "T"])]

    def test_the_variables_of_interest_can_be_replaced(self):
        ds = vtkFilter_DescriptiveStatistics(name="ds", write=True)
        ds.setVariablesOfInterest(["U_x"])
        ds.setVariablesOfInterest(["T"])
        assert ds.params == [("VariablesofInterest", ["T"])]


@pytest.mark.unit
class TestIntegrateVariablesFilter:
    def test_it_declares_the_integrate_variables_filter_type(self):
        iv = vtkFilter_IntegrateVariables(name="iv", write=True, params=[])
        assert iv.filterType == "IntegrateVariables"

    @pytest.mark.xfail(
        strict=True,
        reason="B317: vtkFilter_IntegrateVariables.__init__ takes **kwargs "
               "and forwards them to VTKFilter.__init__, whose `params` is a "
               "required positional argument -- but unlike every one of its "
               "five sibling filter classes (Slice, PlotOverLine, "
               "CellCenters, DescriptiveStatistics, and ExtractBlock via an "
               "explicit params=None default) it omits the "
               "`kwargs.setdefault(\"params\", [])` line, so constructing it "
               "the same way as its siblings raises TypeError.  See the "
               "consolidated findings issue.",
    )
    def test_it_can_be_built_without_an_explicit_parameter_list(self):
        assert vtkFilter_IntegrateVariables(name="iv", write=True).params == []

    def test_it_currently_demands_an_explicit_parameter_list(self):
        """Characterisation of B317."""
        with pytest.raises(TypeError, match="params"):
            vtkFilter_IntegrateVariables(name="iv", write=True)

    def test_the_factory_hides_the_missing_default_by_always_passing_params(self):
        """Characterisation of B317: newVTKPipelineFilter always supplies one."""
        assert VTKPipeLine.newVTKPipelineFilter(name="iv", filterType="IntegrateVariables").params == []


# ---------------------------------------------------------------------------
# registeredVTKPipeLine: binding a pipeline to a case
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestRegisterPipeline:
    def test_it_returns_a_registered_pipeline_bound_to_this_pipeline(self, pipeline, caseDirectory):
        reg = pipeline.registerPipeline(str(caseDirectory))
        assert isinstance(reg, registeredVTKPipeLine)
        assert reg.vtkpipeline is pipeline
        assert reg.datalayer is pipeline.datalayer

    def test_the_case_type_defaults_to_the_decomposed_case(self, registered):
        assert registered.pvOFBase.caseType == CASETYPE_DECOMPOSED

    def test_it_forwards_the_case_type_to_the_paraview_base(self, pipeline, caseDirectory):
        reg = pipeline.registerPipeline(str(caseDirectory), caseType=CASETYPE_RECONSTRUCTED)
        assert reg.pvOFBase.caseType == CASETYPE_RECONSTRUCTED

    def test_the_paraview_base_is_pointed_at_the_case_path(self, registered, caseDirectory):
        assert registered.pvOFBase.casePath == str(caseDirectory)


@pytest.mark.unit
class TestRegisteredPipelineFromDirectory:
    def test_the_case_path_is_the_directory_itself(self, registered, caseDirectory):
        assert registered.casePath == str(caseDirectory)

    def test_no_simulation_document_is_attached(self, registered):
        assert registered.simulationDocument is None

    def test_the_workflow_name_is_the_directory_basename(self, registered):
        assert registered.simulationParams["workflowName"] == "flow_0001"

    def test_the_group_name_is_the_part_before_the_first_underscore(self, registered):
        assert registered.simulationParams["groupName"] == "flow"

    def test_a_directory_without_an_underscore_gets_an_empty_group_name(self, pipeline, tmp_path):
        directory = tmp_path / "plaincase"
        directory.mkdir()
        reg = pipeline.registerPipeline(str(directory))
        assert reg.simulationParams["groupName"] == ""

    def test_a_pseudo_document_has_no_workflow_parameters(self, registered):
        assert registered.simulationParams["workflowParameters"] == {}

    def test_the_default_timestep_block_size_is_fifty(self, registered):
        assert registered.tsBlockNum == 50

    def test_a_path_that_is_neither_in_the_db_nor_a_directory_is_refused(self, pipeline, tmp_path):
        missing = str(tmp_path / "nope")
        with pytest.raises(ValueError, match="is not in the DB, and does not represent a valid case"):
            pipeline.registerPipeline(missing)


@pytest.mark.unit
class TestRegisteredPipelineFromDatabase:
    @pytest.fixture()
    def workflowDocument(self, of, caseDirectory):
        return of.addSimulationsDocument(
            resource=str(caseDirectory),
            dataFormat="string",
            type=of.DOCTYPE_WORKFLOW,
            desc=dict(
                workflowName="myFlow",
                groupName="myGroup",
                parameters={"Parameters": {"alpha": 1.0}},
            ),
        )

    def test_it_finds_the_case_by_workflow_name(self, pipeline, workflowDocument, caseDirectory):
        reg = pipeline.registerPipeline("myFlow")
        assert reg.casePath == str(caseDirectory)

    def test_it_keeps_the_simulation_document_it_found(self, pipeline, workflowDocument):
        reg = pipeline.registerPipeline("myFlow")
        assert reg.simulationDocument is not None
        assert reg.simulationDocument.desc["workflowName"] == "myFlow"

    def test_the_simulation_parameters_come_from_the_document_description(
        self, pipeline, workflowDocument
    ):
        reg = pipeline.registerPipeline("myFlow")
        assert reg.simulationParams == dict(
            workflowName="myFlow",
            groupName="myGroup",
            workflowParameters={"Parameters": {"alpha": 1.0}},
        )

    def test_the_database_branch_wins_over_an_identically_named_directory(
        self, of, pipeline, caseDirectory, tmp_path, monkeypatch
    ):
        (tmp_path / "myFlow").mkdir()
        monkeypatch.chdir(tmp_path)
        of.addSimulationsDocument(
            resource=str(caseDirectory),
            dataFormat="string",
            type=of.DOCTYPE_WORKFLOW,
            desc=dict(workflowName="myFlow", groupName="g", parameters={}),
        )
        reg = pipeline.registerPipeline("myFlow")
        assert reg.casePath == str(caseDirectory)
        assert reg.simulationDocument is not None

    def test_a_workflow_document_missing_a_group_name_is_not_handled(self, of, pipeline, caseDirectory):
        """Characterisation: the three desc keys are read without a guard."""
        of.addSimulationsDocument(
            resource=str(caseDirectory),
            dataFormat="string",
            type=of.DOCTYPE_WORKFLOW,
            desc=dict(workflowName="bare"),
        )
        with pytest.raises(KeyError, match="groupName"):
            pipeline.registerPipeline("bare")


# ---------------------------------------------------------------------------
# registeredVTKPipeLine: the small pure helpers
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestResolveRequestedFilters:
    def test_none_means_every_write_enabled_filter(self, registered):
        assert registered._resolveRequestedFilters(None) == ["ExtractBlock.slice"]

    def test_a_single_name_is_wrapped_in_a_list(self, registered):
        assert registered._resolveRequestedFilters("a.b") == ["a.b"]

    def test_a_list_of_names_is_kept(self, registered):
        assert registered._resolveRequestedFilters(["a", "b"]) == ["a", "b"]


@pytest.mark.unit
class TestParseTimeList:
    CASE_TIMES = [0.0, 10.0, 20.0, 30.0]

    def test_none_means_every_timestep_in_the_case(self, registered):
        assert registered._parseTimeList(None, self.CASE_TIMES) is self.CASE_TIMES

    def test_an_explicit_list_is_passed_through_unchanged(self, registered):
        given = [10.0, 30.0]
        assert registered._parseTimeList(given, self.CASE_TIMES) is given

    def test_a_full_range_selects_the_inclusive_span(self, registered):
        got = registered._parseTimeList("10:20", self.CASE_TIMES)
        assert list(got) == [10.0, 20.0]

    def test_an_open_start_falls_back_to_the_first_case_time(self, registered):
        got = registered._parseTimeList(":20", self.CASE_TIMES)
        assert list(got) == [0.0, 10.0, 20.0]

    def test_an_open_end_falls_back_to_the_last_case_time(self, registered):
        got = registered._parseTimeList("20:", self.CASE_TIMES)
        assert list(got) == [20.0, 30.0]

    def test_a_fully_open_range_selects_everything(self, registered):
        got = registered._parseTimeList(":", self.CASE_TIMES)
        assert list(got) == self.CASE_TIMES

    def test_a_range_that_matches_nothing_yields_no_timesteps(self, registered):
        assert list(registered._parseTimeList("100:200", self.CASE_TIMES)) == []

    def test_a_bare_number_is_read_as_an_open_ended_range(self, registered):
        """Characterisation: "20" is not a single timestep, it is "20:"."""
        assert list(registered._parseTimeList("20", self.CASE_TIMES)) == [20.0, 30.0]

    def test_a_three_part_range_is_not_validated(self, registered):
        """Characterisation: only two bounds exist, so a third part overflows."""
        with pytest.raises(IndexError):
            registered._parseTimeList("0:10:20", self.CASE_TIMES)

    def test_a_non_numeric_bound_is_not_validated(self, registered):
        """Characterisation: bounds go straight into float()."""
        with pytest.raises(ValueError):
            registered._parseTimeList("early:late", self.CASE_TIMES)


@pytest.mark.unit
class TestFilterOutputFileNaming:
    def test_a_regular_mesh_is_written_as_zarr(self, registered):
        assert registered.getFilterOutputFileExt(True) == "zarr"

    def test_a_non_regular_mesh_is_written_as_parquet(self, registered):
        assert registered.getFilterOutputFileExt(False) == "parquet"

    def test_the_first_path_for_a_filter_gets_counter_zero(self, registered, caseDirectory):
        path = registered.getFilterOutputFilePath("slice", "parquet")
        assert path == os.path.join(str(caseDirectory), "vtkpipelinedata", "slice_0.parquet")

    def test_each_request_advances_the_per_filter_counter(self, registered):
        first = registered.getFilterOutputFilePath("slice", "parquet")
        second = registered.getFilterOutputFilePath("slice", "parquet")
        assert os.path.basename(first) == "slice_0.parquet"
        assert os.path.basename(second) == "slice_1.parquet"

    def test_counters_are_kept_separately_per_filter(self, registered):
        registered.getFilterOutputFilePath("first", "parquet")
        second = registered.getFilterOutputFilePath("second", "parquet")
        assert os.path.basename(second) == "second_0.parquet"

    def test_a_dotted_filter_path_becomes_an_underscored_file_name(self, registered):
        path = registered.getFilterOutputFilePath("ExtractBlock.slice", "parquet")
        assert os.path.basename(path) == "ExtractBlock_slice_0.parquet"

    def test_the_extension_is_whatever_the_caller_asked_for(self, registered):
        path = registered.getFilterOutputFilePath("slice", "zarr")
        assert os.path.basename(path) == "slice_0.zarr"

    def test_asking_for_an_existing_path_before_any_exists_returns_nothing(self, registered):
        assert registered.getFilterOutputFilePath("slice", "parquet", generate_new=False) is None

    def test_asking_for_an_existing_path_reuses_the_current_counter(self, registered):
        generated = registered.getFilterOutputFilePath("slice", "parquet")
        reused = registered.getFilterOutputFilePath("slice", "parquet", generate_new=False)
        assert reused == generated

    def test_the_path_is_always_absolute_even_from_a_relative_case_path(
        self, of, pipeline, tmp_path, monkeypatch
    ):
        (tmp_path / "relcase").mkdir()
        monkeypatch.chdir(tmp_path)
        reg = pipeline.registerPipeline("relcase")
        assert os.path.isabs(reg.getFilterOutputFilePath("slice", "parquet"))


@pytest.mark.unit
class TestBuildFilterQuery:
    def test_the_query_carries_the_simulation_the_pipeline_and_the_filter(self, registered):
        qry = registered._buildFilterQuery(filterName="ExtractBlock.slice")
        assert qry["filterName"] == "ExtractBlock.slice"
        assert qry["simulation"]["workflowName"] == "flow_0001"
        assert qry["pipeline"] == registered.vtkpipeline.toJSON()

    def test_the_regular_mesh_flag_is_omitted_when_it_is_not_given(self, registered):
        qry = registered._buildFilterQuery(filterName="f")
        assert "regularMesh" not in qry["simulation"]

    def test_the_regular_mesh_flag_is_recorded_under_the_simulation(self, registered):
        qry = registered._buildFilterQuery(filterName="f", regularMesh=True)
        assert qry["simulation"]["regularMesh"] is True

    def test_the_query_is_usable_as_a_mongo_query(self, registered):
        mongoQuery = dictToMongoQuery(registered._buildFilterQuery(filterName="f"))
        assert mongoQuery["filterName"] == "f"
        assert mongoQuery["simulation__workflowName"] == "flow_0001"

    @pytest.mark.xfail(
        strict=True,
        reason="B318: _buildFilterQuery builds `dict(simulation=self."
               "simulationParams, ...)`, which stores the very same dict "
               "object rather than a copy, and then writes into it "
               "(`qry['simulation']['regularMesh'] = regularMesh`).  Building "
               "a query therefore mutates the object's own simulationParams. "
               "_updateCacheDB compounds it by writing "
               "`recordData['simulation']['timeList'] = timeList` into the "
               "same aliased dict, so after one getData() every later cache "
               "query silently carries a regularMesh and a timeList "
               "constraint that was never asked for.  See the consolidated "
               "findings issue.",
    )
    def test_building_a_query_leaves_the_simulation_parameters_alone(self, registered):
        registered._buildFilterQuery(filterName="f", regularMesh=True)
        assert registered.simulationParams == dict(
            workflowName="flow_0001", groupName="flow", workflowParameters={}
        )

    def test_the_regular_mesh_flag_currently_leaks_into_later_queries(self, registered):
        """Characterisation of B318."""
        registered._buildFilterQuery(filterName="f", regularMesh=True)
        assert registered.simulationParams["regularMesh"] is True
        leaked = registered._buildFilterQuery(filterName="f")
        assert leaked["simulation"]["regularMesh"] is True

    def test_every_query_shares_the_one_simulation_parameters_dict(self, registered):
        """Characterisation of B318."""
        first = registered._buildFilterQuery(filterName="a")
        second = registered._buildFilterQuery(filterName="b")
        assert first["simulation"] is second["simulation"] is registered.simulationParams


# ---------------------------------------------------------------------------
# registeredVTKPipeLine: the DB cache bookkeeping
# ---------------------------------------------------------------------------

def _cacheDescriptionFor(registered, filterName, regularMesh, timeList):
    """A cache-document desc exactly as ``_updateCacheDB`` would write it."""
    desc = copy.deepcopy(registered._buildFilterQuery(filterName=filterName, regularMesh=regularMesh))
    desc["simulation"]["timeList"] = list(timeList)
    return desc


@pytest.mark.unit
class TestFilterCachedTimesteps:
    def test_a_filter_with_no_cache_document_is_scheduled_for_computation(self, registered):
        timeList, toProcess, outputs, docs = registered._filterCachedTimesteps(
            ["f"], [0.0, 10.0], False, "parquet", False
        )
        assert toProcess == ["f"]
        assert list(timeList) == [0.0, 10.0]
        assert docs == {}

    def test_a_filter_with_no_cache_document_gets_a_fresh_output_path(self, registered):
        _, _, outputs, _ = registered._filterCachedTimesteps(
            ["f"], [0.0], False, "parquet", False
        )
        assert os.path.basename(outputs["f"]) == "f_0.parquet"

    def test_a_fully_cached_filter_is_not_recomputed(self, of, registered):
        of.addCacheDocument(
            resource="/somewhere/f_0.parquet",
            dataFormat="parquet",
            type=TYPE_VTK_FILTER,
            desc=_cacheDescriptionFor(registered, "f", False, [0.0, 10.0]),
        )
        timeList, toProcess, outputs, docs = registered._filterCachedTimesteps(
            ["f"], [0.0, 10.0], False, "parquet", False
        )
        assert toProcess == []
        assert list(timeList) == []
        assert outputs["f"] == "/somewhere/f_0.parquet"
        assert list(docs) == ["f"]

    def test_a_partially_cached_filter_is_recomputed_for_the_missing_times_only(self, of, registered):
        of.addCacheDocument(
            resource="/somewhere/f_0.parquet",
            dataFormat="parquet",
            type=TYPE_VTK_FILTER,
            desc=_cacheDescriptionFor(registered, "f", False, [0.0]),
        )
        timeList, toProcess, _, _ = registered._filterCachedTimesteps(
            ["f"], [0.0, 10.0], False, "parquet", False
        )
        assert toProcess == ["f"]
        assert list(timeList) == [10.0]

    def test_overwrite_recomputes_a_fully_cached_filter(self, of, registered):
        of.addCacheDocument(
            resource="/somewhere/f_0.parquet",
            dataFormat="parquet",
            type=TYPE_VTK_FILTER,
            desc=_cacheDescriptionFor(registered, "f", False, [0.0, 10.0]),
        )
        _, toProcess, _, _ = registered._filterCachedTimesteps(
            ["f"], [0.0, 10.0], False, "parquet", True
        )
        assert toProcess == ["f"]

    def test_a_cache_document_for_the_other_mesh_kind_is_not_reused(self, of, registered):
        of.addCacheDocument(
            resource="/somewhere/f_0.zarr",
            dataFormat="zarr_xarray",
            type=TYPE_VTK_FILTER,
            desc=_cacheDescriptionFor(registered, "f", True, [0.0, 10.0]),
        )
        _, toProcess, outputs, docs = registered._filterCachedTimesteps(
            ["f"], [0.0, 10.0], False, "parquet", False
        )
        assert toProcess == ["f"]
        assert docs == {}
        assert outputs["f"].endswith("f_0.parquet")

    def test_a_cache_document_for_another_filter_is_not_reused(self, of, registered):
        of.addCacheDocument(
            resource="/somewhere/other_0.parquet",
            dataFormat="parquet",
            type=TYPE_VTK_FILTER,
            desc=_cacheDescriptionFor(registered, "other", False, [0.0, 10.0]),
        )
        _, toProcess, _, docs = registered._filterCachedTimesteps(
            ["f"], [0.0, 10.0], False, "parquet", False
        )
        assert toProcess == ["f"]
        assert docs == {}

    @pytest.mark.xfail(
        strict=True,
        reason="B319: _filterCachedTimesteps rebinds the single `timeList` "
               "local inside its per-filter loop "
               "(`timeList = [ts for ts in timeList if ts not in dbTimeList]`) "
               "and returns that one list for every filter.  A filter whose "
               "cache is complete therefore empties the timestep list for all "
               "the *other* filters iterated after it, so a filter with no "
               "cache at all is appended to filtersToProcess (len(docList)==0) "
               "and then computed over zero timesteps -- after which "
               "_updateCacheDB records a cache document claiming an empty "
               "timeList for an output file that holds nothing.  Two filters "
               "with disjoint cached times make it order dependent as well: "
               "whichever is visited first is the one recomputed.  See the "
               "consolidated findings issue.",
    )
    def test_a_cached_filter_does_not_steal_the_timesteps_of_an_uncached_one(self, of, registered):
        of.addCacheDocument(
            resource="/somewhere/cached_0.parquet",
            dataFormat="parquet",
            type=TYPE_VTK_FILTER,
            desc=_cacheDescriptionFor(registered, "cached", False, [0.0, 10.0]),
        )
        timeList, toProcess, _, _ = registered._filterCachedTimesteps(
            ["cached", "fresh"], [0.0, 10.0], False, "parquet", False
        )
        assert toProcess == ["fresh"]
        assert list(timeList) == [0.0, 10.0]

    def test_the_timestep_list_is_currently_shared_across_filters(self, of, registered):
        """Characterisation of B319."""
        of.addCacheDocument(
            resource="/somewhere/cached_0.parquet",
            dataFormat="parquet",
            type=TYPE_VTK_FILTER,
            desc=_cacheDescriptionFor(registered, "cached", False, [0.0, 10.0]),
        )
        timeList, toProcess, _, _ = registered._filterCachedTimesteps(
            ["cached", "fresh"], [0.0, 10.0], False, "parquet", False
        )
        assert toProcess == ["fresh"]
        assert list(timeList) == []

    def test_which_fully_cached_filter_is_recomputed_depends_on_the_order(self, of, registered):
        """Characterisation of B319: with two filters cached over disjoint
        times, each one's subtraction leaves the other's times behind, so the
        filter visited first is recomputed even though its cache is complete.
        """
        for filterName, cachedTimes in (("a", [0.0]), ("b", [10.0])):
            of.addCacheDocument(
                resource=f"/somewhere/{filterName}_0.parquet",
                dataFormat="parquet",
                type=TYPE_VTK_FILTER,
                desc=_cacheDescriptionFor(registered, filterName, False, cachedTimes),
            )

        _, firstOrder, _, _ = registered._filterCachedTimesteps(
            ["a", "b"], [0.0, 10.0], False, "parquet", False
        )
        _, secondOrder, _, _ = registered._filterCachedTimesteps(
            ["b", "a"], [0.0, 10.0], False, "parquet", False
        )
        assert firstOrder == ["a"]
        assert secondOrder == ["b"]


@pytest.mark.unit
class TestUpdateCacheDB:
    def test_a_first_computation_creates_a_cache_document(self, of, registered, tmp_path):
        output = tmp_path / "slice_0.parquet"
        registered._updateCacheDB(["slice"], [0.0, 10.0], {}, {"slice": str(output)}, False)
        docList = of.getCacheDocuments(type=TYPE_VTK_FILTER)
        assert len(docList) == 1
        assert docList[0].desc["filterName"] == "slice"
        assert docList[0].desc["simulation"]["timeList"] == [0.0, 10.0]

    def test_the_document_resource_is_the_absolute_output_path(self, of, registered, monkeypatch, tmp_path):
        monkeypatch.chdir(tmp_path)
        registered._updateCacheDB(["slice"], [0.0], {}, {"slice": "relative.parquet"}, False)
        resource = of.getCacheDocuments(type=TYPE_VTK_FILTER)[0].resource
        assert resource == os.path.join(str(tmp_path), "relative.parquet")

    def test_a_non_regular_mesh_is_recorded_as_parquet(self, of, registered, tmp_path):
        registered._updateCacheDB(["slice"], [0.0], {}, {"slice": str(tmp_path / "o")}, False)
        assert of.getCacheDocuments(type=TYPE_VTK_FILTER)[0].dataFormat == of.datatypes.PARQUET

    def test_a_regular_mesh_is_recorded_as_zarr(self, of, registered, tmp_path):
        registered._updateCacheDB(["slice"], [0.0], {}, {"slice": str(tmp_path / "o")}, True)
        assert of.getCacheDocuments(type=TYPE_VTK_FILTER)[0].dataFormat == of.datatypes.ZARR_XARRAY

    def test_the_new_document_is_published_back_to_the_caller(self, registered, tmp_path):
        documents = {}
        registered._updateCacheDB(["slice"], [0.0], documents, {"slice": str(tmp_path / "o")}, False)
        assert list(documents) == ["slice"]

    def test_an_existing_document_gets_the_new_timesteps_merged_in_sorted(self, of, registered, tmp_path):
        registered._updateCacheDB(["slice"], [10.0], {}, {"slice": str(tmp_path / "o")}, False)
        existing = {"slice": of.getCacheDocuments(type=TYPE_VTK_FILTER)[0]}
        registered._updateCacheDB(["slice"], [0.0, 20.0], existing, {}, False)
        assert existing["slice"].desc["simulation"]["timeList"] == [0.0, 10.0, 20.0]

    def test_merging_into_an_existing_document_creates_no_second_document(self, of, registered, tmp_path):
        registered._updateCacheDB(["slice"], [10.0], {}, {"slice": str(tmp_path / "o")}, False)
        existing = {"slice": of.getCacheDocuments(type=TYPE_VTK_FILTER)[0]}
        registered._updateCacheDB(["slice"], [0.0], existing, {}, False)
        assert len(of.getCacheDocuments(type=TYPE_VTK_FILTER)) == 1

    def test_nothing_is_written_when_no_filter_needs_processing(self, of, registered):
        registered._updateCacheDB([], [0.0], {}, {}, False)
        assert list(of.getCacheDocuments(type=TYPE_VTK_FILTER)) == []


@pytest.mark.unit
class TestClearCache:
    def test_it_asks_for_every_filter_when_no_name_is_given(self, registered, monkeypatch):
        seen = []

        def _record(self, **kwargs):
            seen.append(kwargs)
            return []

        monkeypatch.setattr(type(registered.datalayer), "deleteSimulationsDocuments", _record)
        registered.clearCache()
        assert [kw["filterName"] for kw in seen] == ["ExtractBlock", "ExtractBlock.slice"]

    def test_a_single_filter_name_limits_the_deletion(self, registered, monkeypatch):
        seen = []
        monkeypatch.setattr(
            type(registered.datalayer),
            "deleteSimulationsDocuments",
            lambda self, **kwargs: seen.append(kwargs) or [],
        )
        registered.clearCache(filterName="ExtractBlock.slice")
        assert [kw["filterName"] for kw in seen] == ["ExtractBlock.slice"]

    def test_it_always_scopes_the_deletion_to_vtk_filter_documents(self, registered, monkeypatch):
        seen = []
        monkeypatch.setattr(
            type(registered.datalayer),
            "deleteSimulationsDocuments",
            lambda self, **kwargs: seen.append(kwargs) or [],
        )
        registered.clearCache(filterName="f")
        assert seen[0]["type"] == TYPE_VTK_FILTER

    def test_the_regular_mesh_flag_is_added_to_the_deletion_query(self, registered, monkeypatch):
        seen = []
        monkeypatch.setattr(
            type(registered.datalayer),
            "deleteSimulationsDocuments",
            lambda self, **kwargs: seen.append(kwargs) or [],
        )
        registered.clearCache(regularMesh=True, filterName="f")
        assert seen[0]["simulation__regularMesh"] is True

    def test_the_regular_mesh_flag_is_absent_when_it_is_not_given(self, registered, monkeypatch):
        seen = []
        monkeypatch.setattr(
            type(registered.datalayer),
            "deleteSimulationsDocuments",
            lambda self, **kwargs: seen.append(kwargs) or [],
        )
        registered.clearCache(filterName="f")
        assert "simulation__regularMesh" not in seen[0]

    def test_a_deleted_documents_output_file_is_removed_from_disk(self, registered, monkeypatch, tmp_path):
        output = tmp_path / "slice_0.parquet"
        output.write_text("payload")
        doc = dict(
            resource=str(output),
            desc=dict(workflowName="flow_0001", pipeline=dict(filters={})),
        )
        monkeypatch.setattr(
            type(registered.datalayer), "deleteSimulationsDocuments", lambda self, **kwargs: [doc]
        )
        registered.clearCache(filterName="f")
        assert not output.exists()

    def test_a_deleted_documents_output_directory_is_removed_recursively(
        self, registered, monkeypatch, tmp_path
    ):
        output = tmp_path / "slice_0.zarr"
        (output / "inner").mkdir(parents=True)
        doc = dict(
            resource=str(output),
            desc=dict(workflowName="flow_0001", pipeline=dict(filters={})),
        )
        monkeypatch.setattr(
            type(registered.datalayer), "deleteSimulationsDocuments", lambda self, **kwargs: [doc]
        )
        registered.clearCache(filterName="f")
        assert not output.exists()

    def test_an_already_missing_output_file_is_tolerated(self, registered, monkeypatch, tmp_path):
        doc = dict(
            resource=str(tmp_path / "gone.parquet"),
            desc=dict(workflowName="flow_0001", pipeline=dict(filters={})),
        )
        monkeypatch.setattr(
            type(registered.datalayer), "deleteSimulationsDocuments", lambda self, **kwargs: [doc]
        )
        registered.clearCache(filterName="f")

    @pytest.mark.xfail(
        strict=True,
        reason="B320: clearCache deletes through "
               "`self.datalayer.deleteSimulationsDocuments(...)` -- the "
               "Simulations collection -- but the cache it is meant to clear "
               "is written by _updateCacheDB with `addCacheDocument` and read "
               "back by _filterCachedTimesteps with `getCacheDocuments`, i.e. "
               "the Cache collection.  The query it builds is correct (the "
               "test asserts the same query does match through "
               "getCacheDocuments); only the collection is wrong, so "
               "clearCache is a no-op against real pipeline output and leaves "
               "both the documents and the files on disk.  See the "
               "consolidated findings issue.",
    )
    def test_it_removes_the_cache_documents_the_pipeline_itself_wrote(self, of, registered, tmp_path):
        output = tmp_path / "slice_0.parquet"
        output.write_text("payload")
        desc = _cacheDescriptionFor(registered, "ExtractBlock.slice", False, [0.0])
        of.addCacheDocument(
            resource=str(output), dataFormat="parquet", type=TYPE_VTK_FILTER, desc=desc
        )
        registered.clearCache(regularMesh=False, filterName="ExtractBlock.slice")
        assert list(of.getCacheDocuments(type=TYPE_VTK_FILTER)) == []

    def test_the_cache_documents_currently_survive_the_clear(self, of, registered, tmp_path):
        """Characterisation of B320."""
        output = tmp_path / "slice_0.parquet"
        output.write_text("payload")
        desc = _cacheDescriptionFor(registered, "ExtractBlock.slice", False, [0.0])
        of.addCacheDocument(
            resource=str(output), dataFormat="parquet", type=TYPE_VTK_FILTER, desc=desc
        )
        query = dictToMongoQuery(
            registered._buildFilterQuery(filterName="ExtractBlock.slice", regularMesh=False)
        )
        assert len(of.getCacheDocuments(type=TYPE_VTK_FILTER, **query)) == 1

        registered.clearCache(regularMesh=False, filterName="ExtractBlock.slice")

        assert len(of.getCacheDocuments(type=TYPE_VTK_FILTER)) == 1
        assert output.exists()

    @pytest.mark.xfail(
        strict=True,
        reason="B321: clearCache's per-document debug line reads "
               "`doc['desc']['workflowName']`, but the desc that this very "
               "class writes (_updateCacheDB -> _buildFilterQuery) nests the "
               "name as desc['simulation']['workflowName'] and has no "
               "top-level workflowName at all.  The f-string is built eagerly, "
               "regardless of log level, so the moment the delete query "
               "actually matches one of the pipeline's own cache documents the "
               "loop dies with KeyError before removing anything from disk.  "
               "Sibling of B101 in toolkit.clearVTKPipelineCache.  See the "
               "consolidated findings issue.",
    )
    def test_it_removes_the_output_file_of_a_document_it_wrote_itself(
        self, registered, monkeypatch, tmp_path
    ):
        output = tmp_path / "slice_0.parquet"
        output.write_text("payload")
        doc = dict(
            resource=str(output),
            desc=_cacheDescriptionFor(registered, "ExtractBlock.slice", False, [0.0]),
        )
        monkeypatch.setattr(
            type(registered.datalayer), "deleteSimulationsDocuments", lambda self, **kwargs: [doc]
        )
        registered.clearCache(filterName="ExtractBlock.slice")
        assert not output.exists()

    def test_a_document_written_by_this_class_currently_raises_on_clear(
        self, registered, monkeypatch, tmp_path
    ):
        """Characterisation of B321."""
        output = tmp_path / "slice_0.parquet"
        output.write_text("payload")
        doc = dict(
            resource=str(output),
            desc=_cacheDescriptionFor(registered, "ExtractBlock.slice", False, [0.0]),
        )
        assert "workflowName" not in doc["desc"]
        assert "workflowName" in doc["desc"]["simulation"]
        monkeypatch.setattr(
            type(registered.datalayer), "deleteSimulationsDocuments", lambda self, **kwargs: [doc]
        )
        with pytest.raises(KeyError, match="workflowName"):
            registered.clearCache(filterName="ExtractBlock.slice")
        assert output.exists()


# ---------------------------------------------------------------------------
# _buildFilterLayer: the ParaView filter tree (structure only)
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestBuildFilterLayer:
    def test_a_none_structure_creates_nothing(self, registered, paraviewLog):
        assert registered._buildFilterLayer(None, "READER", None) == []
        assert paraviewLog == []

    def test_an_empty_structure_creates_nothing(self, registered, paraviewLog):
        assert registered._buildFilterLayer(None, "READER", {}) == []
        assert paraviewLog == []

    def test_it_returns_the_full_names_of_the_filters_it_created(self, registered, paraviewLog):
        names = registered._buildFilterLayer(
            None, "READER", registered.vtkpipeline.toJSON()["filters"]
        )
        assert names == ["ExtractBlock", "ExtractBlock.slice"]

    def test_each_filter_is_instantiated_by_its_paraview_type_name(self, registered, paraviewLog):
        registered._buildFilterLayer(None, "READER", registered.vtkpipeline.toJSON()["filters"])
        created = [(entry[1], entry[2]) for entry in paraviewLog if entry[0] == "create"]
        assert created == [("ExtractBlock", "ExtractBlock"), ("Slice", "ExtractBlock.slice")]

    def test_a_root_filter_is_wired_to_the_father_it_was_given(self, registered, paraviewLog):
        registered._buildFilterLayer(None, "READER", registered.vtkpipeline.toJSON()["filters"])
        rootCreate = [entry for entry in paraviewLog if entry[0] == "create"][0]
        assert rootCreate[3] == "READER"

    def test_a_downstream_filter_is_wired_to_its_father_proxy(self, registered, paraviewLog):
        registered._buildFilterLayer(None, "READER", registered.vtkpipeline.toJSON()["filters"])
        creates = [entry for entry in paraviewLog if entry[0] == "create"]
        fatherProxy = creates[1][3]
        assert isinstance(fatherProxy, _RecordingProxy)
        assert repr(fatherProxy) == "<proxy ExtractBlock>"

    def test_a_father_name_prefixes_the_gui_names(self, registered, paraviewLog):
        names = registered._buildFilterLayer(
            "root", "READER", registered.vtkpipeline.toJSON()["filters"]
        )
        assert names == ["root.ExtractBlock", "root.ExtractBlock.slice"]

    def test_the_parameters_are_applied_in_the_declared_order(self, registered, paraviewLog):
        registered._buildFilterLayer(None, "READER", registered.vtkpipeline.toJSON()["filters"])
        sets = [entry for entry in paraviewLog if entry[0] == "set"]
        assert sets == [
            ("set", "ExtractBlock", "Selectors", []),
            ("set", "ExtractBlock.slice", "SliceType", "Plane"),
            ("set", "ExtractBlock.slice.SliceType", "Origin", [0, 0, 10]),
            ("set", "ExtractBlock.slice.SliceType", "Normal", [0, 0, 1]),
        ]

    def test_a_dotted_parameter_name_is_set_on_the_nested_sub_proxy(self, registered, paraviewLog):
        structure = {
            "s": dict(filterType="Slice", write=True, params=[("a.b.c", 5)], downstream={})
        }
        registered._buildFilterLayer(None, "READER", structure)
        assert ("set", "s.a.b", "c", 5) in paraviewLog

    def test_every_filter_is_updated_after_its_parameters_are_set(self, registered, paraviewLog):
        registered._buildFilterLayer(None, "READER", registered.vtkpipeline.toJSON()["filters"])
        assert [entry for entry in paraviewLog if entry[0] == "update"] == [
            ("update", "ExtractBlock"),
            ("update", "ExtractBlock.slice"),
        ]

    def test_a_father_is_created_and_updated_before_its_child(self, registered, paraviewLog):
        registered._buildFilterLayer(None, "READER", registered.vtkpipeline.toJSON()["filters"])
        order = [(entry[0], entry[1]) for entry in paraviewLog if entry[0] in ("create", "update")]
        assert order == [
            ("create", "ExtractBlock"),
            ("update", "ExtractBlock"),
            ("create", "Slice"),
            ("update", "ExtractBlock.slice"),
        ]

    def test_siblings_are_built_in_declaration_order(self, registered, paraviewLog):
        structure = {
            "first": dict(filterType="CellCenters", write=True, params=[], downstream={}),
            "second": dict(filterType="Slice", write=True, params=[], downstream={}),
        }
        assert registered._buildFilterLayer(None, "READER", structure) == ["first", "second"]

    def test_a_filter_type_that_paraview_does_not_expose_is_not_caught(self, registered, paraviewLog):
        """Characterisation: the type name goes straight into getattr."""
        structure = {"x": dict(filterType="NoSuchParaViewFilter", write=True, params=[])}
        with pytest.raises(AttributeError):
            registered._buildFilterLayer(None, "READER", structure)

    def test_a_structure_without_a_params_key_is_not_caught(self, registered, paraviewLog):
        """Characterisation: 'params' and 'filterType' are read without a guard."""
        with pytest.raises(KeyError, match="params"):
            registered._buildFilterLayer(None, "READER", {"x": dict(filterType="Slice")})


# ---------------------------------------------------------------------------
# _buildAndExecuteParaViewPipeline: argument forwarding only
# ---------------------------------------------------------------------------

@pytest.fixture()
def paraviewExecution(monkeypatch, paraviewLog):
    """Patch the two ParaView entry points of paraviewOpenFOAM and record them.

    ``initializeReader``/``writeCase`` dereference a ``pvsimple`` name that
    is never bound in ``pvOpenFOAMBase`` when ``vtkmodules`` is missing, so
    the real ones cannot run here at all.  The recorder is the observable.
    """
    calls = {"readers": [], "writeCase": []}

    def _initializeReader(self, readerName="reader"):
        calls["readers"].append(readerName)
        return _RecordingProxy(paraviewLog, f"__{readerName}__")

    def _writeCase(self, **kwargs):
        calls["writeCase"].append(kwargs)

    monkeypatch.setattr(paraviewOpenFOAM, "initializeReader", _initializeReader)
    monkeypatch.setattr(paraviewOpenFOAM, "writeCase", _writeCase)
    return calls


@pytest.mark.unit
class TestBuildAndExecuteParaViewPipeline:
    def test_it_creates_the_reader_named_reader(self, registered, paraviewExecution):
        registered._buildAndExecuteParaViewPipeline(
            ["ExtractBlock.slice"], {"ExtractBlock.slice": "/out.parquet"}, [0.0], None, False, False
        )
        assert paraviewExecution["readers"] == ["reader"]

    def test_the_root_filters_are_wired_to_the_reader(self, registered, paraviewExecution, paraviewLog):
        registered._buildAndExecuteParaViewPipeline(
            ["ExtractBlock.slice"], {"ExtractBlock.slice": "/out.parquet"}, [0.0], None, False, False
        )
        rootCreate = [entry for entry in paraviewLog if entry[0] == "create"][0]
        assert repr(rootCreate[3]) == "<proxy __reader__>"

    def test_the_reader_is_narrowed_to_the_requested_field_names(
        self, registered, paraviewExecution, paraviewLog
    ):
        registered._buildAndExecuteParaViewPipeline(
            ["ExtractBlock.slice"], {"ExtractBlock.slice": "/out.parquet"}, [0.0], ["U", "T"], False, False
        )
        assert ("set", "__reader__", "CellArrays", ["U", "T"]) in paraviewLog

    def test_the_reader_is_left_alone_when_no_field_names_are_given(
        self, registered, paraviewExecution, paraviewLog
    ):
        registered._buildAndExecuteParaViewPipeline(
            ["ExtractBlock.slice"], {"ExtractBlock.slice": "/out.parquet"}, [0.0], None, False, False
        )
        assert [entry for entry in paraviewLog if entry[1] == "__reader__"] == []

    def test_only_the_filters_to_process_are_handed_to_the_writer(self, registered, paraviewExecution):
        outputs = {"ExtractBlock": "/block.parquet", "ExtractBlock.slice": "/slice.parquet"}
        registered._buildAndExecuteParaViewPipeline(
            ["ExtractBlock.slice"], outputs, [0.0], None, False, False
        )
        assert paraviewExecution["writeCase"][0]["filtersDict"] == {
            "ExtractBlock.slice": "/slice.parquet"
        }

    def test_the_returned_map_is_the_one_handed_to_the_writer(self, registered, paraviewExecution):
        outputs = {"ExtractBlock.slice": "/slice.parquet"}
        returned = registered._buildAndExecuteParaViewPipeline(
            ["ExtractBlock.slice"], outputs, [0.0], None, False, False
        )
        assert returned == paraviewExecution["writeCase"][0]["filtersDict"]

    def test_the_time_list_fields_overwrite_and_mesh_kind_reach_the_writer(
        self, registered, paraviewExecution
    ):
        registered._buildAndExecuteParaViewPipeline(
            ["ExtractBlock.slice"], {"ExtractBlock.slice": "/out.zarr"}, [0.0, 10.0], ["U"], True, True
        )
        kwargs = paraviewExecution["writeCase"][0]
        assert kwargs["timeList"] == [0.0, 10.0]
        assert kwargs["fieldnames"] == ["U"]
        assert kwargs["overwrite"] is True
        assert kwargs["regularMesh"] is True

    def test_the_pipelines_timestep_block_size_reaches_the_writer(self, registered, paraviewExecution):
        registered.tsBlockNum = 7
        registered._buildAndExecuteParaViewPipeline(
            ["ExtractBlock.slice"], {"ExtractBlock.slice": "/out.parquet"}, [0.0], None, False, False
        )
        assert paraviewExecution["writeCase"][0]["tsBlockNum"] == 7

    def test_every_paraview_source_is_deleted_afterwards(
        self, registered, paraviewExecution, paraviewLog, monkeypatch
    ):
        monkeypatch.setattr(
            pipelineModule.pvsimple, "GetSources", lambda: {"a": "proxyA", "b": "proxyB"}
        )
        registered._buildAndExecuteParaViewPipeline(
            ["ExtractBlock.slice"], {"ExtractBlock.slice": "/out.parquet"}, [0.0], None, False, False
        )
        assert [entry for entry in paraviewLog if entry[0] == "delete"] == [
            ("delete", "proxyA"),
            ("delete", "proxyB"),
        ]

    def test_a_filter_to_process_with_no_output_path_is_not_caught(self, registered, paraviewExecution):
        """Characterisation: filtersOutputFilename is indexed without a guard."""
        with pytest.raises(KeyError, match="unknown"):
            registered._buildAndExecuteParaViewPipeline(
                ["unknown"], {}, [0.0], None, False, False
            )


# ---------------------------------------------------------------------------
# getData and its two convenience wrappers
# ---------------------------------------------------------------------------

@pytest.fixture()
def executedPipeline(registered, monkeypatch, paraviewLog):
    """``registered`` with a fixed case time list and a writer that really writes.

    ``OFToolkit.getTimeList`` is patched on the class so the timesteps are
    deterministic (the real one walks the case directory), and
    ``writeCase`` is replaced by a stub that writes a genuine parquet file
    per filter so that the final ``document.getData()`` reload is real.
    """
    calls = {"writeCase": []}

    monkeypatch.setattr(
        type(registered.datalayer), "getTimeList", lambda self, *args, **kwargs: [0.0, 10.0]
    )

    def _initializeReader(self, readerName="reader"):
        return _RecordingProxy(paraviewLog, f"__{readerName}__")

    def _writeCase(self, **kwargs):
        calls["writeCase"].append(kwargs)
        for filterName, outputPath in kwargs["filtersDict"].items():
            os.makedirs(os.path.dirname(outputPath), exist_ok=True)
            pandas.DataFrame({"time": kwargs["timeList"], "filter": filterName}).to_parquet(
                outputPath
            )

    monkeypatch.setattr(paraviewOpenFOAM, "initializeReader", _initializeReader)
    monkeypatch.setattr(paraviewOpenFOAM, "writeCase", _writeCase)
    return registered, calls


@pytest.mark.unit
class TestGetData:
    def test_it_returns_one_entry_per_write_enabled_filter(self, executedPipeline):
        registered, _ = executedPipeline
        assert list(registered.getData(regularMesh=False)) == ["ExtractBlock.slice"]

    def test_the_returned_value_is_whatever_the_writer_produced(self, executedPipeline):
        registered, _ = executedPipeline
        frame = registered.getData(regularMesh=False)["ExtractBlock.slice"].compute()
        assert sorted(frame["time"]) == [0.0, 10.0]

    def test_it_records_the_computed_timesteps_in_the_cache(self, of, executedPipeline):
        registered, _ = executedPipeline
        registered.getData(regularMesh=False)
        docList = of.getCacheDocuments(type=TYPE_VTK_FILTER)
        assert len(docList) == 1
        assert docList[0].desc["filterName"] == "ExtractBlock.slice"
        assert docList[0].desc["simulation"]["timeList"] == [0.0, 10.0]

    def test_the_cached_document_points_at_the_file_that_was_written(self, of, executedPipeline):
        registered, _ = executedPipeline
        registered.getData(regularMesh=False)
        resource = of.getCacheDocuments(type=TYPE_VTK_FILTER)[0].resource
        assert os.path.exists(resource)
        assert os.path.basename(resource) == "ExtractBlock_slice_0.parquet"

    def test_the_case_time_list_is_used_when_no_times_are_requested(self, executedPipeline):
        registered, calls = executedPipeline
        registered.getData(regularMesh=False)
        assert calls["writeCase"][0]["timeList"] == [0.0, 10.0]

    def test_the_latest_time_flag_keeps_only_the_last_timestep(self, executedPipeline):
        registered, calls = executedPipeline
        registered.getData(regularMesh=False, latestTime=True)
        assert calls["writeCase"][0]["timeList"] == [10.0]

    def test_an_explicit_time_list_is_honoured(self, executedPipeline):
        registered, calls = executedPipeline
        registered.getData(regularMesh=False, timeList=[10.0])
        assert calls["writeCase"][0]["timeList"] == [10.0]

    def test_a_single_named_filter_overrides_the_write_flags(self, executedPipeline):
        registered, calls = executedPipeline
        result = registered.getData(regularMesh=False, filterName="ExtractBlock")
        assert list(result) == ["ExtractBlock"]
        assert list(calls["writeCase"][0]["filtersDict"]) == ["ExtractBlock"]

    def test_the_field_names_are_passed_to_the_writer(self, executedPipeline):
        registered, calls = executedPipeline
        registered.getData(regularMesh=False, fieldNames=["U"])
        assert calls["writeCase"][0]["fieldnames"] == ["U"]

    def test_a_second_identical_request_is_served_from_the_cache(self, of, executedPipeline):
        registered, calls = executedPipeline
        registered.getData(regularMesh=False)
        registered.getData(regularMesh=False)
        assert len(calls["writeCase"]) == 1
        assert len(of.getCacheDocuments(type=TYPE_VTK_FILTER)) == 1

    def test_overwrite_recomputes_an_already_cached_filter(self, executedPipeline):
        registered, calls = executedPipeline
        registered.getData(regularMesh=False)
        registered.getData(regularMesh=False, overwrite=True)
        assert len(calls["writeCase"]) == 2
        assert calls["writeCase"][1]["overwrite"] is True

    def test_the_paraview_pipeline_is_not_built_when_nothing_needs_computing(
        self, executedPipeline, paraviewLog
    ):
        registered, _ = executedPipeline
        registered.getData(regularMesh=False)
        del paraviewLog[:]
        registered.getData(regularMesh=False)
        assert paraviewLog == []

    def test_a_filter_starved_of_timesteps_is_still_recorded_as_cached(self, of, executedPipeline):
        """Characterisation of B319 as it surfaces through getData: the
        already-cached filter empties the timestep list, so the brand-new
        filter requested alongside it is "computed" over zero timesteps and
        still gets a cache document claiming so.
        """
        registered, calls = executedPipeline
        registered.getData(regularMesh=False)
        registered.getData(regularMesh=False, filterName=["ExtractBlock.slice", "ghost"])

        assert calls["writeCase"][1]["timeList"] == []
        ghost = [
            doc
            for doc in of.getCacheDocuments(type=TYPE_VTK_FILTER)
            if doc.desc["filterName"] == "ghost"
        ]
        assert len(ghost) == 1
        assert ghost[0].desc["simulation"]["timeList"] == []


@pytest.mark.unit
class TestGetDataWrappers:
    @pytest.fixture()
    def recordedGetData(self, monkeypatch):
        calls = []
        monkeypatch.setattr(
            registeredVTKPipeLine, "getData", lambda self, **kwargs: calls.append(kwargs)
        )
        return calls

    def test_get_regular_data_asks_for_a_regular_mesh(self, registered, recordedGetData):
        registered.getRegularData()
        assert recordedGetData[0]["regularMesh"] is True

    def test_get_non_regular_data_asks_for_an_irregular_mesh(self, registered, recordedGetData):
        registered.getNonRegularData()
        assert recordedGetData[0]["regularMesh"] is False

    def test_get_regular_data_forwards_every_argument_it_takes(self, registered, recordedGetData):
        registered.getRegularData(
            filterName="f", timeList=[1.0], fieldNames=["U"], overwrite=True
        )
        assert recordedGetData[0] == dict(
            regularMesh=True, filterName="f", timeList=[1.0], fieldNames=["U"], overwrite=True
        )

    def test_get_non_regular_data_forwards_every_argument_it_takes(self, registered, recordedGetData):
        registered.getNonRegularData(
            filterName="f", timeList=[1.0], fieldNames=["U"], overwrite=True
        )
        assert recordedGetData[0] == dict(
            regularMesh=False, filterName="f", timeList=[1.0], fieldNames=["U"], overwrite=True
        )

    def test_neither_wrapper_exposes_the_latest_time_flag(self, registered, recordedGetData):
        """Characterisation: latestTime is reachable only through getData."""
        registered.getRegularData()
        assert "latestTime" not in recordedGetData[0]
        with pytest.raises(TypeError):
            registered.getRegularData(latestTime=True)
