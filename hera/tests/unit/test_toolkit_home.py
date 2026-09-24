"""ToolkitHome: the registry that hands out toolkits and remembers repositories.

CLAUDE.md forbids instantiating toolkits directly, so everything goes through
here.  That makes the registry's behaviour -- which class you get, what a bad
name does, where the repository default lives -- part of the public contract.
"""
import pytest

from hera import toolkitHome
from hera.datalayer import datatypes

PROJECT = "UNIT_TEST_PROJECT"


@pytest.mark.unit
class TestGetToolkit:
    def test_returns_the_expected_class(self, unit_toolkit_factory):
        from hera.measurements.GIS.raster.topography import TopographyToolkit

        toolkit = unit_toolkit_factory(toolkitHome.GIS_RASTER_TOPOGRAPHY)
        assert isinstance(toolkit, TopographyToolkit)

    def test_the_toolkit_knows_its_own_name(self, unit_toolkit_factory):
        toolkit = unit_toolkit_factory(toolkitHome.GIS_RASTER_TOPOGRAPHY)
        assert toolkit.toolkitName == "TopographyToolkit"

    def test_the_toolkit_is_bound_to_the_requested_project(self, unit_toolkit_factory):
        toolkit = unit_toolkit_factory(toolkitHome.GIS_RASTER_TOPOGRAPHY)
        assert toolkit.projectName == PROJECT

    def test_two_calls_give_independent_instances(self, unit_toolkit_factory):
        """No caching promised, so two calls must not alias one object."""
        first = unit_toolkit_factory(toolkitHome.GIS_RASTER_TOPOGRAPHY)
        second = unit_toolkit_factory(toolkitHome.GIS_RASTER_TOPOGRAPHY)
        assert first is not second

    def test_different_constants_give_different_classes(self, unit_toolkit_factory):
        topography = unit_toolkit_factory(toolkitHome.GIS_RASTER_TOPOGRAPHY)
        landcover = unit_toolkit_factory(toolkitHome.GIS_LANDCOVER)
        assert type(topography) is not type(landcover)

    def test_an_unknown_name_is_reported_rather_than_returning_none(
        self, unit_files_directory
    ):
        """A silent None here becomes an AttributeError somewhere far away."""
        with pytest.raises(Exception) as raised:
            toolkitHome.getToolkit(
                "NoSuchToolkitAnywhere",
                projectName=PROJECT,
                filesDirectory=unit_files_directory,
            )
        assert "NoSuchToolkitAnywhere" in str(raised.value)


@pytest.mark.unit
class TestToolkitNameConstants:
    """The constants exist so callers never spell a toolkit name by hand."""

    CONSTANTS = [
        "GIS_RASTER_TOPOGRAPHY",
        "GIS_LANDCOVER",
        "GIS_BUILDINGS",
        "GIS_DEMOGRAPHY",
        "GIS_TILES",
        "GIS_VECTOR_TOPOGRAPHY",
        "METEOROLOGY_LOWFREQ",
        "METEOROLOGY_HIGHFREQ",
        "RISKASSESSMENT",
        "LSM",
        "GAUSSIANDISPERSION",
        "WINDPROFILE",
        "SIMULATIONS_OPENFOAM",
        "SIMULATIONS_WORKFLOWS",
        "EXPERIMENT",
        "DATA",
        "MACHINELEARNING_DEEPLEARNING",
    ]

    @pytest.mark.parametrize("name", CONSTANTS)
    def test_the_constant_exists_and_is_a_string(self, name):
        assert isinstance(getattr(toolkitHome, name), str)

    def test_the_constants_are_distinct(self):
        values = [getattr(toolkitHome, name) for name in self.CONSTANTS]
        assert len(values) == len(set(values))

    def test_the_default_project_constant_is_exposed(self):
        assert toolkitHome.DEFAULTPROJECT == "defaultProject"

    @pytest.mark.parametrize(
        "name, expected",
        [
            ("TOOLKIT_SAVEMODE_NOSAVE", None),
            ("TOOLKIT_SAVEMODE_ONLYFILE", "File"),
            ("TOOLKIT_SAVEMODE_ONLYFILE_REPLACE", "File_overwrite"),
            ("TOOLKIT_SAVEMODE_FILEANDDB", "DB"),
            ("TOOLKIT_SAVEMODE_FILEANDDB_REPLACE", "DB_overwrite"),
        ],
    )
    def test_save_mode_constants_keep_their_values(self, name, expected):
        """These strings are persisted, so changing one is a data migration."""
        assert getattr(toolkitHome, name) == expected


@pytest.mark.unit
class TestDefaultRepository:
    """Documented as persisting a repository name FOR A PROJECT.

    Neither half of that works today, and the tests below say which is which:
    B21 -- ToolkitHome._get_data_toolkit declares a projectName parameter and
           ignores it, returning dataToolkit(), which always operates on the
           default project.  So the setting is global, not per project.
    B22 -- because the write therefore lands on defaultProject, which
           Project guards as read-only (project.py:757), every call raises.
           registerToolkit sets _allowWritingToDefaultProject before writing
           (toolkit.py:1271); setDefaultRepository does not.
    """

    def test_unset_reads_as_an_empty_string(self):
        assert toolkitHome.getDefaultRepository(projectName="NO_SUCH_PROJECT") == ""

    def test_setting_without_a_project_name_raises(self):
        with pytest.raises(ValueError, match="'projectName' is required"):
            toolkitHome.setDefaultRepository(projectName="", repositoryName="x")

    def test_setting_without_a_repository_name_raises(self):
        with pytest.raises(ValueError, match="'repositoryName' is required"):
            toolkitHome.setDefaultRepository(projectName=PROJECT, repositoryName="")

    def test_reading_without_a_project_name_raises(self):
        with pytest.raises(ValueError, match="'projectName' is required"):
            toolkitHome.getDefaultRepository(projectName="")

    def test_the_helper_ignores_the_project_it_is_given(self):
        """Characterisation of B21's mechanism, so the cause is unambiguous."""
        fromProject = toolkitHome._get_data_toolkit(projectName=PROJECT)
        fromNothing = toolkitHome._get_data_toolkit()
        assert fromProject.projectName == fromNothing.projectName
        assert fromProject.projectName == toolkitHome.DEFAULTPROJECT

    def test_the_json_format_lookup_is_dead(self):
        """setDefaultRepository tries three names that do not exist.

        `getattr(datatypes, "JSON") or getattr(datatypes, "json") or
        getattr(datatypes, "TEXT")` -- the real constants are JSON_DICT,
        JSON_PANDAS and JSON_GEOPANDAS, so the block never sets a format.
        Reported as dead code rather than a wrong result: the payload lives in
        desc, so the format is not load-bearing here.
        """
        for name in ("JSON", "json", "TEXT"):
            assert not hasattr(datatypes, name), (
                f"datatypes.{name} now exists, so the lookup in "
                "setDefaultRepository is no longer dead -- revisit it"
            )

    @pytest.mark.xfail(
        strict=True,
        reason="B22: setDefaultRepository writes through a dataToolkit bound to "
               "defaultProject and does not set _allowWritingToDefaultProject, so "
               "it always raises RuntimeError: project defaultProject is "
               "read-only. See the consolidated findings issue.",
    )
    def test_round_trip(self):
        toolkitHome.setDefaultRepository(
            projectName=PROJECT, repositoryName="meteo_data_v1"
        )
        assert toolkitHome.getDefaultRepository(projectName=PROJECT) == "meteo_data_v1"

    @pytest.mark.xfail(
        strict=True,
        reason="B21: the setting is not scoped to a project, because "
               "_get_data_toolkit ignores its projectName argument. Blocked behind "
               "B22, which makes the write fail first. "
               "See the consolidated findings issue.",
    )
    def test_defaults_are_per_project(self):
        toolkitHome.setDefaultRepository(projectName=PROJECT, repositoryName="mine")
        assert toolkitHome.getDefaultRepository(projectName="OTHER_PROJECT") == ""


@pytest.mark.unit
class TestRegisterToolkit:
    def test_an_empty_name_is_rejected(self, tmp_path):
        with pytest.raises(ValueError, match="toolkit_name must be provided"):
            toolkitHome.registerToolkit(toolkit_name="", toolkit_path=str(tmp_path))

    def test_an_empty_path_is_rejected(self):
        with pytest.raises(ValueError, match="toolkit_path must be provided"):
            toolkitHome.registerToolkit(toolkit_name="MyToolkit", toolkit_path="")

    def test_a_directory_is_stored_as_given(self, tmp_path):
        document = toolkitHome.registerToolkit(
            toolkit_name="DirToolkit", toolkit_path=str(tmp_path)
        )
        assert document.resource == str(tmp_path)

    def test_a_file_is_stored_as_its_parent_directory(self, tmp_path):
        """getToolkit puts the resource on sys.path, so it must be a directory."""
        module = tmp_path / "MyToolkit.py"
        module.write_text("class MyToolkit: pass\n", encoding="utf-8")

        document = toolkitHome.registerToolkit(
            toolkit_name="FileToolkit", toolkit_path=str(module)
        )
        assert document.resource == str(tmp_path)

    def test_a_relative_path_is_made_absolute(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "sub").mkdir()

        document = toolkitHome.registerToolkit(
            toolkit_name="RelToolkit", toolkit_path="sub"
        )
        assert document.resource == str(tmp_path / "sub")

    def test_params_are_stored_in_the_description(self, tmp_path):
        document = toolkitHome.registerToolkit(
            toolkit_name="ParamToolkit",
            toolkit_path=str(tmp_path),
            params={"alpha": 1},
        )
        assert document.desc["params"] == {"alpha": 1}

    def test_omitted_params_become_an_empty_dict(self, tmp_path):
        document = toolkitHome.registerToolkit(
            toolkit_name="NoParamToolkit", toolkit_path=str(tmp_path)
        )
        assert document.desc["params"] == {}

    def test_registering_twice_without_overwrite_raises(self, tmp_path):
        toolkitHome.registerToolkit(toolkit_name="TwiceToolkit", toolkit_path=str(tmp_path))
        with pytest.raises(ValueError, match="already exists"):
            toolkitHome.registerToolkit(
                toolkit_name="TwiceToolkit", toolkit_path=str(tmp_path)
            )

    def test_overwrite_replaces_the_registration(self, tmp_path):
        other = tmp_path / "other"
        other.mkdir()

        toolkitHome.registerToolkit(toolkit_name="OverToolkit", toolkit_path=str(tmp_path))
        document = toolkitHome.registerToolkit(
            toolkit_name="OverToolkit", toolkit_path=str(other), overwrite=True
        )
        assert document.resource == str(other)

    def test_the_registration_is_findable_afterwards(self, tmp_path):
        toolkitHome.registerToolkit(
            toolkit_name="FindableToolkit", toolkit_path=str(tmp_path)
        )
        from hera.utils.data.toolkit import dataToolkit

        found = dataToolkit().getDataSourceDocument("FindableToolkit")
        assert found is not None
        assert found.resource == str(tmp_path)
