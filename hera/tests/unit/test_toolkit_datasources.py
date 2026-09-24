"""abstractToolkit: data-source management against an in-memory database.

These are not mock tests.  The mongomock seam means the real datalayer runs --
query construction, desc filtering, version resolution, document tagging --
only the storage engine is in memory.  hera/toolkit.py was the largest gap in
this area at 30% coverage.
"""
import pytest

from hera import toolkitHome
from hera.toolkit import TOOLKIT_DATASOURCE_TYPE

TOPOGRAPHY = toolkitHome.GIS_RASTER_TOPOGRAPHY
LANDCOVER = toolkitHome.GIS_LANDCOVER


@pytest.fixture()
def toolkit(unit_toolkit_factory):
    return unit_toolkit_factory(TOPOGRAPHY)


@pytest.mark.unit
class TestAddDataSource:
    def test_creates_a_retrievable_document(self, toolkit):
        toolkit.addDataSource("SRC", "/data/a.tif", "GeoTIFF")
        document = toolkit.getDataSourceDocument("SRC")
        assert document is not None
        assert document.resource == "/data/a.tif"
        assert document.dataFormat == "GeoTIFF"

    def test_default_version_is_0_0_1(self, toolkit):
        toolkit.addDataSource("SRC", "/data/a.tif", "GeoTIFF")
        assert toolkit.getDataSourceDocument("SRC").desc["version"] == [0, 0, 1]

    def test_document_type_is_the_datasource_type(self, toolkit):
        toolkit.addDataSource("SRC", "/data/a.tif", "GeoTIFF")
        assert toolkit.getDataSourceDocument("SRC")["type"] == TOOLKIT_DATASOURCE_TYPE

    def test_toolkit_name_is_tagged_automatically(self, toolkit):
        toolkit.addDataSource("SRC", "/data/a.tif", "GeoTIFF")
        assert toolkit.getDataSourceDocument("SRC").desc["toolkit"] == toolkit.toolkitName

    def test_extra_metadata_reaches_the_description(self, toolkit):
        toolkit.addDataSource("SRC", "/data/a.tif", "GeoTIFF", station="YAVNEEL", crs=2039)
        desc = toolkit.getDataSourceDocument("SRC").desc
        assert desc["station"] == "YAVNEEL"
        assert desc["crs"] == 2039

    def test_extra_metadata_is_queryable(self, toolkit):
        toolkit.addDataSource("A", "/data/a.tif", "GeoTIFF", station="YAVNEEL")
        toolkit.addDataSource("B", "/data/b.tif", "GeoTIFF", station="TELAVIV")
        assert toolkit.getDataSourceList(station="YAVNEEL") == ["A"]

    def test_duplicate_without_overwrite_raises(self, toolkit):
        toolkit.addDataSource("SRC", "/data/a.tif", "GeoTIFF")
        with pytest.raises(ValueError, match="already exists"):
            toolkit.addDataSource("SRC", "/data/other.tif", "GeoTIFF")

    def test_duplicate_error_names_the_source_and_version(self, toolkit):
        toolkit.addDataSource("SRC", "/data/a.tif", "GeoTIFF", version=(1, 2, 3))
        with pytest.raises(ValueError, match=r"SRC.*\(1, 2, 3\)"):
            toolkit.addDataSource("SRC", "/data/a.tif", "GeoTIFF", version=(1, 2, 3))

    def test_overwrite_replaces_rather_than_duplicating(self, toolkit):
        toolkit.addDataSource("SRC", "/data/a.tif", "GeoTIFF")
        toolkit.addDataSource("SRC", "/data/new.tif", "GeoTIFF", overwrite=True)

        documents = toolkit.getDataSourceDocumentsList(datasourceName="SRC")
        assert len(documents) == 1, "overwrite left a duplicate behind"
        assert documents[0].resource == "/data/new.tif"

    def test_a_second_version_does_not_replace_the_first(self, toolkit):
        toolkit.addDataSource("SRC", "/data/v1.tif", "GeoTIFF", version=(0, 0, 1))
        toolkit.addDataSource("SRC", "/data/v2.tif", "GeoTIFF", version=(0, 0, 2))
        assert len(toolkit.getDataSourceDocumentsList(datasourceName="SRC")) == 2


@pytest.mark.unit
class TestVersionResolution:
    """With no explicit version, the latest wins -- compared componentwise."""

    def test_latest_of_several_is_returned(self, toolkit):
        for version in [(0, 0, 1), (0, 0, 2), (0, 0, 3)]:
            toolkit.addDataSource("SRC", f"/data/{version}.tif", "GeoTIFF", version=version)
        assert toolkit.getDataSourceDocument("SRC").desc["version"] == [0, 0, 3]

    def test_ordering_is_numeric_not_lexicographic(self, toolkit):
        """[0,0,10] must beat [0,0,9]; string ordering would get this wrong."""
        for version in [(0, 0, 9), (0, 0, 10)]:
            toolkit.addDataSource("SRC", f"/data/{version}.tif", "GeoTIFF", version=version)
        assert toolkit.getDataSourceDocument("SRC").desc["version"] == [0, 0, 10]

    def test_minor_outranks_patch(self, toolkit):
        for version in [(0, 0, 99), (0, 1, 0)]:
            toolkit.addDataSource("SRC", f"/data/{version}.tif", "GeoTIFF", version=version)
        assert toolkit.getDataSourceDocument("SRC").desc["version"] == [0, 1, 0]

    def test_major_outranks_minor(self, toolkit):
        for version in [(0, 99, 99), (1, 0, 0)]:
            toolkit.addDataSource("SRC", f"/data/{version}.tif", "GeoTIFF", version=version)
        assert toolkit.getDataSourceDocument("SRC").desc["version"] == [1, 0, 0]

    def test_explicit_version_overrides_the_latest(self, toolkit):
        for version in [(0, 0, 1), (0, 0, 2)]:
            toolkit.addDataSource("SRC", f"/data/{version}.tif", "GeoTIFF", version=version)
        document = toolkit.getDataSourceDocument("SRC", version=(0, 0, 1))
        assert document.resource == "/data/(0, 0, 1).tif"

    def test_tuple_and_list_versions_are_interchangeable(self, toolkit):
        """They travel through different query shapes and must still agree.

        dictToMongoQuery indexes lists into version__0/__1/__2 but passes a
        tuple through as a single equality, so these are two distinct queries
        for the same thing.
        """
        toolkit.addDataSource("SRC", "/data/a.tif", "GeoTIFF", version=(0, 0, 7))
        assert toolkit.getDataSourceDocument("SRC", version=(0, 0, 7)) is not None
        assert toolkit.getDataSourceDocument("SRC", version=[0, 0, 7]) is not None

    def test_version_is_stored_as_a_list(self, toolkit):
        """BSON has no tuple, so a tuple comes back as a list."""
        toolkit.addDataSource("SRC", "/data/a.tif", "GeoTIFF", version=(0, 0, 1))
        assert toolkit.getDataSourceDocument("SRC").desc["version"] == [0, 0, 1]

    def test_unknown_source_returns_none(self, toolkit):
        assert toolkit.getDataSourceDocument("NOTHERE") is None

    def test_unknown_version_returns_none(self, toolkit):
        toolkit.addDataSource("SRC", "/data/a.tif", "GeoTIFF", version=(0, 0, 1))
        assert toolkit.getDataSourceDocument("SRC", version=(9, 9, 9)) is None

    def test_lookup_does_not_write_to_the_database(self, toolkit):
        """A getter must not persist anything.

        The docstring still promises that the latest version is "automatically
        persisted as the default in the project config"; the code removed that
        deliberately, with a comment calling it a hidden side effect.  The code
        is right and the docstring is stale -- this test pins the code.
        """
        for version in [(0, 0, 1), (0, 0, 2)]:
            toolkit.addDataSource("SRC", f"/data/{version}.tif", "GeoTIFF", version=version)

        before = dict(toolkit.getConfig())
        toolkit.getDataSourceDocument("SRC")
        assert dict(toolkit.getConfig()) == before


@pytest.mark.unit
class TestDefaultVersion:
    def test_setting_a_default_changes_what_lookup_returns(self, toolkit):
        for version in [(0, 0, 1), (0, 0, 2)]:
            toolkit.addDataSource("SRC", f"/data/{version}.tif", "GeoTIFF", version=version)

        toolkit.setDataSourceDefaultVersion("SRC", [0, 0, 1])
        assert toolkit.getDataSourceDocument("SRC").desc["version"] == [0, 0, 1]

    def test_the_default_is_recorded_in_the_project_config(self, toolkit):
        toolkit.addDataSource("SRC", "/data/a.tif", "GeoTIFF", version=(0, 0, 1))
        toolkit.setDataSourceDefaultVersion("SRC", [0, 0, 1])
        assert toolkit.getConfig()["SRC_defaultVersion"] == [0, 0, 1]

    def test_an_explicit_version_still_wins_over_the_default(self, toolkit):
        for version in [(0, 0, 1), (0, 0, 2)]:
            toolkit.addDataSource("SRC", f"/data/{version}.tif", "GeoTIFF", version=version)
        toolkit.setDataSourceDefaultVersion("SRC", [0, 0, 1])
        document = toolkit.getDataSourceDocument("SRC", version=(0, 0, 2))
        assert document.desc["version"] == [0, 0, 2]

    def test_setting_a_default_for_a_missing_version_raises(self, toolkit):
        toolkit.addDataSource("SRC", "/data/a.tif", "GeoTIFF", version=(0, 0, 1))
        with pytest.raises(ValueError, match="No DataSource with name=SRC"):
            toolkit.setDataSourceDefaultVersion("SRC", [9, 9, 9])

    def test_setting_a_default_for_a_missing_source_raises(self, toolkit):
        with pytest.raises(ValueError, match="No DataSource with name=NOTHERE"):
            toolkit.setDataSourceDefaultVersion("NOTHERE", [0, 0, 1])

    def test_defaults_are_per_source(self, toolkit):
        toolkit.addDataSource("A", "/data/a.tif", "GeoTIFF", version=(0, 0, 1))
        toolkit.addDataSource("A", "/data/a2.tif", "GeoTIFF", version=(0, 0, 2))
        toolkit.addDataSource("B", "/data/b.tif", "GeoTIFF", version=(0, 0, 1))
        toolkit.addDataSource("B", "/data/b2.tif", "GeoTIFF", version=(0, 0, 2))

        toolkit.setDataSourceDefaultVersion("A", [0, 0, 1])
        assert toolkit.getDataSourceDocument("A").desc["version"] == [0, 0, 1]
        assert toolkit.getDataSourceDocument("B").desc["version"] == [0, 0, 2]


@pytest.mark.unit
class TestListingDataSources:
    def test_map_reports_format_resource_and_description(self, toolkit):
        toolkit.addDataSource("SRC", "/data/a.tif", "GeoTIFF", station="YAVNEEL")
        entry = toolkit.getDataSourceMap()[0]
        assert entry["dataFormat"] == "GeoTIFF"
        assert entry["resource"] == "/data/a.tif"
        assert entry["station"] == "YAVNEEL"
        assert entry["datasourceName"] == "SRC"

    def test_map_has_one_entry_per_version(self, toolkit):
        """Documented as "all data sources AND their versions"."""
        for version in [(0, 0, 1), (0, 0, 2)]:
            toolkit.addDataSource("SRC", f"/data/{version}.tif", "GeoTIFF", version=version)
        assert len(toolkit.getDataSourceMap()) == 2

    def test_empty_project_gives_empty_listings(self, toolkit):
        assert toolkit.getDataSourceList() == []
        assert toolkit.getDataSourceMap() == []

    def test_table_is_empty_for_an_empty_project(self, toolkit):
        table = toolkit.getDataSourceTable()
        assert table.empty

    def test_table_has_a_row_per_document(self, toolkit):
        toolkit.addDataSource("A", "/data/a.tif", "GeoTIFF")
        toolkit.addDataSource("B", "/data/b.tif", "GeoTIFF")
        table = toolkit.getDataSourceTable()
        assert len(table) == 2
        assert set(table["datasourceName"]) == {"A", "B"}

    def test_documents_helper_returns_a_single_element_list(self, toolkit):
        toolkit.addDataSource("SRC", "/data/a.tif", "GeoTIFF")
        assert len(toolkit.getDataSourceDocuments("SRC")) == 1

    def test_documents_helper_returns_empty_for_a_miss(self, toolkit):
        assert toolkit.getDataSourceDocuments("NOTHERE") == []

    @pytest.mark.xfail(
        strict=True,
        reason="B18: getDataSourceList is documented as returning 'data source "
               "names' but yields one entry per document, so a source with four "
               "versions appears four times. getDataSourceMap is the one "
               "documented as covering versions. See the consolidated findings issue.",
    )
    def test_list_reports_each_name_once(self, toolkit):
        for version in [(0, 0, 1), (0, 0, 2), (0, 0, 3)]:
            toolkit.addDataSource("SRC", f"/data/{version}.tif", "GeoTIFF", version=version)
        assert toolkit.getDataSourceList() == ["SRC"]


@pytest.mark.unit
class TestDeleteDataSource:
    def test_deleting_removes_the_document(self, toolkit):
        toolkit.addDataSource("SRC", "/data/a.tif", "GeoTIFF")
        toolkit.deleteDataSource("SRC")
        assert toolkit.getDataSourceDocument("SRC") is None

    def test_deleting_returns_the_removed_document(self, toolkit):
        toolkit.addDataSource("SRC", "/data/a.tif", "GeoTIFF")
        removed = toolkit.deleteDataSource("SRC")
        assert removed.resource == "/data/a.tif"

    def test_deleting_something_absent_returns_none(self, toolkit):
        assert toolkit.deleteDataSource("NOTHERE") is None

    def test_deleting_one_version_leaves_the_others(self, toolkit):
        for version in [(0, 0, 1), (0, 0, 2)]:
            toolkit.addDataSource("SRC", f"/data/{version}.tif", "GeoTIFF", version=version)

        toolkit.deleteDataSource("SRC", version=(0, 0, 1))
        remaining = toolkit.getDataSourceDocumentsList(datasourceName="SRC")
        assert len(remaining) == 1
        assert remaining[0].desc["version"] == [0, 0, 2]

    def test_delete_then_add_the_same_name_is_allowed(self, toolkit):
        toolkit.addDataSource("SRC", "/data/a.tif", "GeoTIFF")
        toolkit.deleteDataSource("SRC")
        toolkit.addDataSource("SRC", "/data/b.tif", "GeoTIFF")
        assert toolkit.getDataSourceDocument("SRC").resource == "/data/b.tif"


@pytest.mark.unit
class TestToolkitIsolation:
    """Two toolkits share a project; neither may see the other's sources."""

    def test_sources_are_scoped_to_their_toolkit(self, unit_toolkit_factory):
        topography = unit_toolkit_factory(TOPOGRAPHY)
        landcover = unit_toolkit_factory(LANDCOVER)

        topography.addDataSource("SHARED_NAME", "/data/topo.tif", "GeoTIFF")
        landcover.addDataSource("SHARED_NAME", "/data/land.tif", "GeoTIFF")

        assert topography.getDataSourceDocument("SHARED_NAME").resource == "/data/topo.tif"
        assert landcover.getDataSourceDocument("SHARED_NAME").resource == "/data/land.tif"

    def test_listings_do_not_leak_across_toolkits(self, unit_toolkit_factory):
        topography = unit_toolkit_factory(TOPOGRAPHY)
        landcover = unit_toolkit_factory(LANDCOVER)

        topography.addDataSource("ONLY_TOPO", "/data/topo.tif", "GeoTIFF")
        assert topography.getDataSourceList() == ["ONLY_TOPO"]
        assert landcover.getDataSourceList() == []

    def test_deleting_in_one_toolkit_leaves_the_other(self, unit_toolkit_factory):
        topography = unit_toolkit_factory(TOPOGRAPHY)
        landcover = unit_toolkit_factory(LANDCOVER)

        topography.addDataSource("SHARED_NAME", "/data/topo.tif", "GeoTIFF")
        landcover.addDataSource("SHARED_NAME", "/data/land.tif", "GeoTIFF")

        topography.deleteDataSource("SHARED_NAME")
        assert topography.getDataSourceDocument("SHARED_NAME") is None
        assert landcover.getDataSourceDocument("SHARED_NAME") is not None
