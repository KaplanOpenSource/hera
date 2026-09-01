"""abstractToolkit's document-adding overrides.

The three overrides exist to tag every document with the toolkit's name, which
is what keeps two toolkits sharing one project from seeing each other's data.
That makes the tagging a correctness property, not a convenience.
"""
import pytest

from hera import toolkitHome
from hera.toolkit import TOOLKIT_TOOLKITNAME_FIELD

TOPOGRAPHY = toolkitHome.GIS_RASTER_TOPOGRAPHY
LANDCOVER = toolkitHome.GIS_LANDCOVER


@pytest.fixture()
def toolkit(unit_toolkit_factory):
    return unit_toolkit_factory(TOPOGRAPHY)


@pytest.mark.unit
class TestAutomaticTagging:
    def test_cache_documents_are_tagged(self, toolkit):
        toolkit.addCacheDocument(resource="/x", dataFormat="string", type="T")
        document = toolkit.getCacheDocuments(type="T")[0]
        assert document.desc[TOOLKIT_TOOLKITNAME_FIELD] == toolkit.toolkitName

    def test_measurements_documents_are_tagged(self, toolkit):
        toolkit.addMeasurementsDocument(resource="/x", dataFormat="string", type="T")
        document = toolkit.getMeasurementsDocuments(type="T")[0]
        assert document.desc[TOOLKIT_TOOLKITNAME_FIELD] == toolkit.toolkitName

    def test_simulations_documents_are_tagged(self, toolkit):
        toolkit.addSimulationsDocument(resource="/x", dataFormat="string", type="T")
        document = toolkit.getSimulationsDocuments(type="T")[0]
        assert document.desc[TOOLKIT_TOOLKITNAME_FIELD] == toolkit.toolkitName

    def test_supplied_metadata_is_preserved_alongside_the_tag(self, toolkit):
        toolkit.addCacheDocument(
            resource="/x", dataFormat="string", type="T", desc={"station": "A"}
        )
        desc = toolkit.getCacheDocuments(type="T")[0].desc
        assert desc["station"] == "A"
        assert desc[TOOLKIT_TOOLKITNAME_FIELD] == toolkit.toolkitName

    def test_an_explicit_toolkit_field_is_not_overwritten(self, toolkit):
        """setdefault, not assignment -- a caller may name the toolkit itself."""
        toolkit.addCacheDocument(
            resource="/x",
            dataFormat="string",
            type="T",
            desc={TOOLKIT_TOOLKITNAME_FIELD: "SomethingElse"},
        )
        document = toolkit.getCacheDocuments(type="T")[0]
        assert document.desc[TOOLKIT_TOOLKITNAME_FIELD] == "SomethingElse"

    def test_omitting_desc_entirely_still_works(self, toolkit):
        """desc defaults to None and must become a fresh dict, not a shared one."""
        toolkit.addCacheDocument(resource="/a", dataFormat="string", type="T")
        toolkit.addCacheDocument(resource="/b", dataFormat="string", type="T")
        assert len(toolkit.getCacheDocuments(type="T")) == 2

    def test_the_three_collections_stay_separate(self, toolkit):
        """CLAUDE.md: simulation results must never land in Measurements."""
        toolkit.addCacheDocument(resource="/c", dataFormat="string", type="T")
        toolkit.addMeasurementsDocument(resource="/m", dataFormat="string", type="T")
        toolkit.addSimulationsDocument(resource="/s", dataFormat="string", type="T")

        assert [d.resource for d in toolkit.getCacheDocuments(type="T")] == ["/c"]
        assert [d.resource for d in toolkit.getMeasurementsDocuments(type="T")] == ["/m"]
        assert [d.resource for d in toolkit.getSimulationsDocuments(type="T")] == ["/s"]


@pytest.mark.unit
class TestCallerDictIsNotMutated:
    """The overrides must not write into the dict they were handed.

    This is not style: because the tag is applied with setdefault, a mutated
    dict carries the first toolkit's name into every later call, and the
    document ends up attributed to the wrong toolkit.
    """

    @pytest.mark.xfail(
        strict=True,
        reason="B19: addCacheDocument/addMeasurementsDocument/addSimulationsDocument "
               "call desc.setdefault on the caller's dict, injecting a 'toolkit' key "
               "into it. See the consolidated findings issue.",
    )
    def test_the_supplied_dict_is_left_alone(self, toolkit):
        metadata = {"station": "A"}
        toolkit.addCacheDocument(
            resource="/x", dataFormat="string", type="T", desc=metadata
        )
        assert metadata == {"station": "A"}

    @pytest.mark.xfail(
        strict=True,
        reason="B19: a metadata dict reused across two toolkits keeps the first "
               "toolkit's tag, because setdefault will not overwrite the key the "
               "first call injected. The second toolkit's document is stored under "
               "the wrong toolkit name. See the consolidated findings issue.",
    )
    def test_reusing_a_metadata_dict_does_not_misattribute_documents(
        self, unit_toolkit_factory
    ):
        """The scenario, verified: two toolkits, one shared metadata dict.

            topography.addCacheDocument(desc=metadata)   # injects 'TopographyToolkit'
            landcover.addCacheDocument(desc=metadata)    # keeps  'TopographyToolkit'

        The landcover document is then invisible to landcover's own queries and
        shows up in topography's instead.
        """
        topography = unit_toolkit_factory(TOPOGRAPHY)
        landcover = unit_toolkit_factory(LANDCOVER)

        metadata = {"station": "A"}
        topography.addCacheDocument(
            resource="/topo", dataFormat="string", type="T", desc=metadata
        )
        landcover.addCacheDocument(
            resource="/land", dataFormat="string", type="T", desc=metadata
        )

        stored = {
            document.resource: document.desc[TOOLKIT_TOOLKITNAME_FIELD]
            for document in landcover.getCacheDocuments(type="T")
        }
        assert stored["/land"] == landcover.toolkitName


@pytest.mark.unit
class TestToolkitIdentity:
    def test_toolkit_name_is_exposed(self, toolkit):
        assert toolkit.toolkitName == "TopographyToolkit"

    def test_project_name_is_exposed(self, toolkit):
        assert toolkit.projectName == "UNIT_TEST_PROJECT"

    def test_analysis_layer_is_attached(self, toolkit):
        assert toolkit.analysis is not None

    def test_analysis_is_not_a_toolkit_subclass(self, toolkit):
        """CLAUDE.md: analysis and presentation are composition, not inheritance."""
        from hera.toolkit import abstractToolkit

        assert not isinstance(toolkit.analysis, abstractToolkit)

    def test_a_toolkit_is_a_project(self, toolkit):
        """The inheritance that makes every toolkit need a database."""
        from hera.datalayer import Project

        assert isinstance(toolkit, Project)

    def test_two_toolkits_report_different_names(self, unit_toolkit_factory):
        topography = unit_toolkit_factory(TOPOGRAPHY)
        landcover = unit_toolkit_factory(LANDCOVER)
        assert topography.toolkitName != landcover.toolkitName


@pytest.mark.unit
class TestProjectConfig:
    def test_config_round_trips(self, toolkit):
        toolkit.setConfig(alpha=1, beta="two")
        config = toolkit.getConfig()
        assert config["alpha"] == 1
        assert config["beta"] == "two"

    def test_setting_config_twice_merges_rather_than_replacing(self, toolkit):
        toolkit.setConfig(alpha=1)
        toolkit.setConfig(beta=2)
        config = toolkit.getConfig()
        assert config["alpha"] == 1 and config["beta"] == 2

    def test_config_is_shared_by_toolkits_in_one_project(self, unit_toolkit_factory):
        """Config is per project, not per toolkit -- worth knowing before relying on it."""
        topography = unit_toolkit_factory(TOPOGRAPHY)
        landcover = unit_toolkit_factory(LANDCOVER)

        topography.setConfig(shared="value")
        assert landcover.getConfig()["shared"] == "value"
