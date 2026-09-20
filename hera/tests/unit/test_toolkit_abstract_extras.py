"""hera/toolkit.py: three small abstractToolkit methods not covered
elsewhere -- presentation (the property, exercised directly rather than
through a concrete toolkit's own presentation layer), classLoggerName,
and getToolkitDocument (a second, independent implementation of the same
lookup idea as hera.utils.data.toolkit_repository.ToolkitRepository,
covered in an earlier batch)."""
import pytest


@pytest.mark.unit
class TestAbstractToolkitExtras:
    def test_presentation_defaults_to_none(self, unit_toolkit_factory):
        from hera import toolkitHome

        toolkit = unit_toolkit_factory(toolkitHome.GIS_RASTER_TOPOGRAPHY)
        assert toolkit.presentation is None

    def test_class_logger_name_includes_the_full_qualified_class_path(self, unit_toolkit_factory):
        from hera import toolkitHome

        toolkit = unit_toolkit_factory(toolkitHome.GIS_RASTER_TOPOGRAPHY)
        assert "TopographyToolkit" in toolkit.classLoggerName
        assert toolkit.classLoggerName.endswith("{the_function_name}")

    def test_get_toolkit_document_returns_none_when_unregistered(self, unit_toolkit_factory):
        from hera import toolkitHome

        toolkit = unit_toolkit_factory(toolkitHome.GIS_RASTER_TOPOGRAPHY)
        assert toolkit.getToolkitDocument("noSuchToolkit") is None

    def test_get_toolkit_document_finds_a_registered_datasource(self, unit_toolkit_factory):
        from hera import toolkitHome

        toolkit = unit_toolkit_factory(toolkitHome.GIS_RASTER_TOPOGRAPHY)
        toolkit.addMeasurementsDocument(
            resource="/x", dataFormat="string", type="ToolkitDataSource",
            desc={"datasourceName": "myTool", "toolkit": "myTool"},
        )
        found = toolkit.getToolkitDocument("myTool")
        assert found is not None
        assert found.desc["datasourceName"] == "myTool"
