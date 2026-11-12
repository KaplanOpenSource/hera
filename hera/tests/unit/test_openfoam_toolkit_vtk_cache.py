"""openFoam/toolkit.py: OFToolkit construction and the VTK pipeline cache
query/table/clear methods. Now exercisable with real hermes+argos on the
path (see conftest -- the stub layer prefers a real local install).
"""
import pytest


@pytest.fixture()
def of(unit_toolkit_factory):
    from hera import toolkitHome

    return unit_toolkit_factory(toolkitHome.SIMULATIONS_OPENFOAM)


def _add_filter_doc(of, **desc_overrides):
    desc = dict(regularMesh=True, filterName="f1", simulation=dict(workflowName="wf1", groupName="g1"))
    desc.update(desc_overrides)
    of.addCacheDocument(resource="/x", dataFormat="string", type="vtk_filter", desc=desc)


@pytest.mark.unit
class TestOFToolkitConstruction:
    def test_it_builds_the_analysis_and_presentation_layers(self, of):
        assert of.analysis is not None
        assert of.presentation is not None

    def test_it_builds_the_object_home(self, of):
        assert of.OFObjectHome is not None


@pytest.mark.unit
class TestVTKPipelineCache:
    def test_an_empty_project_has_no_cached_documents(self, of):
        assert list(of.getVTKPipelineCacheDocuments()) == []

    def test_it_filters_by_regular_mesh(self, of):
        _add_filter_doc(of, regularMesh=True)
        assert len(of.getVTKPipelineCacheDocuments(regularMesh=True)) == 1
        assert len(of.getVTKPipelineCacheDocuments(regularMesh=False)) == 0

    def test_it_filters_by_workflow_and_group_name(self, of):
        _add_filter_doc(of)
        assert len(of.getVTKPipelineCacheDocuments(workflowName="wf1")) == 1
        assert len(of.getVTKPipelineCacheDocuments(workflowName="noSuchWorkflow")) == 0
        assert len(of.getVTKPipelineCacheDocuments(groupName="g1")) == 1

    def test_the_table_lists_filter_workflow_and_group_names(self, of):
        _add_filter_doc(of)
        table = of.getVTKPipelineCacheTable()
        assert list(table.iloc[0][["filterName", "workflowName", "groupName"]]) == ["f1", "wf1", "g1"]

    @pytest.mark.xfail(
        strict=True,
        reason="B101: clearVTKPipelineCache's debug log line reads "
               "doc['desc']['pipeline']['filters'] unconditionally (an "
               "f-string is built eagerly regardless of log level), but "
               "'pipeline' is not a required/documented desc key -- any "
               "cache doc without one (e.g. one built via addCacheDocument "
               "directly, matching this toolkit's own documented schema for "
               "regularMesh/filterName/simulation) crashes with KeyError. "
               "See the consolidated findings issue.",
    )
    def test_clear_removes_the_matching_cache_documents(self, of):
        _add_filter_doc(of, regularMesh=True)
        of.clearVTKPipelineCache(regularMesh=True)
        assert list(of.getVTKPipelineCacheDocuments(regularMesh=True)) == []

    def test_clear_currently_raises_without_a_pipeline_key(self, of):
        """Characterisation of B101."""
        _add_filter_doc(of, regularMesh=True)
        with pytest.raises(KeyError, match="pipeline"):
            of.clearVTKPipelineCache(regularMesh=True)
