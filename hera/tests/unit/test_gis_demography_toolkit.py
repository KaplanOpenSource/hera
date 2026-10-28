"""DemographyToolkit: construction, config plumbing, and the deprecated
polygon-population wrapper.

Constructing it also constructs a real ``BuildingsToolkit`` (mongomock
handles both). The heavier analysis/presentation pipelines
(``calculatePopulationInPolygon``, all the plotting methods) need real
population shapefiles and are left for an integration test.

One defect:

* B81: the ``shapes`` property returns ``self._shapes``, but nothing in
  the class -- not a class attribute, not ``__init__`` -- ever sets
  ``_shapes``. Every access raises ``AttributeError``.
"""
import pytest

from hera import toolkitHome


@pytest.fixture()
def demography(unit_toolkit_factory):
    return unit_toolkit_factory(toolkitHome.GIS_DEMOGRAPHY)


@pytest.mark.unit
class TestConstruction:
    def test_it_also_builds_a_real_buildings_toolkit(self, demography):
        from hera.measurements.GIS.vector.buildings.toolkit import BuildingsToolkit

        assert isinstance(demography.buildings, BuildingsToolkit)

    def test_the_population_type_map_matches_the_documented_groups(self, demography):
        assert demography.populationTypes == {
            "All": "total_pop", "Children": "age_0_14", "Youth": "age_15_19",
            "YoungAdults": "age_20_29", "Adults": "age_30_64", "Elderly": "age_65_up",
        }

    def test_analysis_and_presentation_layers_point_back_at_the_toolkit(self, demography):
        assert demography.analysis.datalayer is demography
        assert demography.presentation.datalayer is demography


@pytest.mark.unit
class TestShapesIsBroken:
    @pytest.mark.xfail(
        strict=True,
        reason="B81: the shapes property returns self._shapes, which "
               "nothing in the class ever assigns -- neither a class "
               "attribute nor anywhere in __init__. Every access raises "
               "AttributeError. See the consolidated findings issue.",
    )
    def test_shapes_should_return_something(self, demography):
        demography.shapes

    def test_shapes_currently_always_raises(self, demography):
        """Characterisation of B81."""
        with pytest.raises(AttributeError, match="_shapes"):
            demography.shapes


@pytest.mark.unit
class TestFilesDirectory:
    def test_set_default_directory_creates_a_missing_directory(self, demography, tmp_path):
        target = tmp_path / "new" / "nested"
        demography.setDefaultDirectory(str(target), create=True)
        assert target.is_dir()
        assert demography.FilesDirectory == str(target)

    def test_set_default_directory_without_create_raises_on_a_missing_path(self, demography, tmp_path):
        target = tmp_path / "does_not_exist"
        with pytest.raises(NotADirectoryError):
            demography.setDefaultDirectory(str(target), create=False)

    def test_an_existing_directory_is_accepted_without_create(self, demography, tmp_path):
        demography.setDefaultDirectory(str(tmp_path), create=False)
        assert demography.FilesDirectory == str(tmp_path)


@pytest.mark.unit
class TestLoadData:
    def test_it_registers_a_datasource_with_the_geopandas_format(self, demography, tmp_path):
        from hera.datalayer import datatypes

        demography.loadData("myRegion", str(tmp_path / "data.geojson"), version=(0, 0, 1))
        sources = demography.getDataSourceMap()
        assert sources[0]["datasourceName"] == "myRegion"
        assert sources[0]["dataFormat"] == datatypes.GEOPANDAS
        assert sources[0]["version"] == [0, 0, 1]

    def test_extra_metadata_is_passed_through_to_the_datasource(self, demography, tmp_path):
        demography.loadData(
            "myRegion", str(tmp_path / "data.geojson"), metadata={"population_year": 2020}
        )
        sources = demography.getDataSourceMap()
        assert sources[0]["population_year"] == 2020


@pytest.mark.unit
class TestProjectPolygonOnPopulationDelegates:
    def test_it_warns_that_it_is_deprecated(self, demography):
        with pytest.warns(DeprecationWarning, match="calculatePopulationInPolygon"):
            demography.analysis.calculatePopulationInPolygon = lambda **kwargs: None
            demography.projectPolygonOnPopulation(shapelyPolygon="POLY", dataSourceOrData="DS")

    def test_it_forwards_all_arguments_to_the_analysis_layer(self, demography):
        received = {}
        demography.analysis.calculatePopulationInPolygon = lambda **kwargs: received.update(kwargs)
        demography.projectPolygonOnPopulation(
            shapelyPolygon="POLY", dataSourceOrData="DS", populationTypes="Children", dataSourceVersion=(1, 0, 0)
        )
        assert received == {
            "shapelyPolygon": "POLY", "dataSourceOrData": "DS",
            "populationTypes": "Children", "dataSourceVersion": (1, 0, 0),
        }
