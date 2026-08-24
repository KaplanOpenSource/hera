"""The seam must route the entire Project/toolkit stack at mongomock."""
import pytest


@pytest.mark.unit
def test_project_constructs_without_mongodb(unit_project):
    assert unit_project.projectName == "UNIT_TEST_PROJECT"


@pytest.mark.unit
def test_document_roundtrip(unit_project):
    unit_project.addMeasurementsDocument(
        resource="/synthetic/a.parquet",
        dataFormat="parquet",
        type="Probe",
        desc=dict(station="A"),
    )
    found = unit_project.getMeasurementsDocuments(type="Probe", station="A")
    assert len(found) == 1
    assert found[0].resource == "/synthetic/a.parquet"


@pytest.mark.unit
def test_config_roundtrip(unit_project):
    unit_project.setConfig(alpha=1)
    assert unit_project.getConfig()["alpha"] == 1


@pytest.mark.unit
def test_real_toolkit_is_usable(unit_toolkit_factory):
    """A real toolkit, built through toolkitHome, backed by mongomock."""
    from hera import toolkitHome

    toolkit = unit_toolkit_factory(toolkitHome.GIS_RASTER_TOPOGRAPHY)
    toolkit.addDataSource("SRC", "/rel/path.tif", "GeoTIFF", version=[0, 0, 1])
    names = [d.desc.get("datasourceName") for d in toolkit.getDataSourceDocuments("SRC")]
    assert names == ["SRC"]


@pytest.mark.unit
def test_module_level_singletons_are_rebound():
    """abstractcalculator.py uses datalayer.Cache directly; it must be live."""
    import getpass

    import hera.datalayer as dl
    from hera.datalayer.document import getDBObject

    assert dl.Cache._metadataCol is getDBObject("Cache", getpass.getuser())
