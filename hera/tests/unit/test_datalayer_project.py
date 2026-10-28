"""Project: configuration, counters, and document loading from dicts.

Everything here runs the real Project against the in-memory store, so
constructor side effects, config persistence and counter arithmetic are
exercised for real -- only the storage engine is substituted.
"""
import json

import pytest

from hera.datalayer import Project, createProjectDirectory


@pytest.mark.unit
class TestConfig:
    def test_config_starts_empty_apart_from_the_files_directory(self, unit_project):
        """Project.__init__ persists filesDirectory, so that key is expected."""
        assert set(unit_project.getConfig()) <= {"filesDirectory"}

    def test_set_then_get(self, unit_project):
        unit_project.setConfig(alpha=1)
        assert unit_project.getConfig()["alpha"] == 1

    def test_setting_a_second_key_keeps_the_first(self, unit_project):
        unit_project.setConfig(alpha=1)
        unit_project.setConfig(beta=2)
        config = unit_project.getConfig()
        assert config["alpha"] == 1 and config["beta"] == 2

    def test_setting_the_same_key_twice_takes_the_new_value(self, unit_project):
        unit_project.setConfig(alpha=1)
        unit_project.setConfig(alpha=2)
        assert unit_project.getConfig()["alpha"] == 2

    def test_nested_values_survive(self, unit_project):
        unit_project.setConfig(nested={"a": [1, 2]})
        assert unit_project.getConfig()["nested"] == {"a": [1, 2]}

    def test_init_config_does_not_overwrite_an_existing_key(self, unit_project):
        """Documented: initConfig sets only keys that do not exist yet."""
        unit_project.setConfig(alpha=1)
        unit_project.initConfig(alpha=99)
        assert unit_project.getConfig()["alpha"] == 1

    def test_init_config_sets_a_missing_key(self, unit_project):
        unit_project.initConfig(fresh=5)
        assert unit_project.getConfig()["fresh"] == 5

    def test_config_survives_reopening_the_project(self, unit_files_directory):
        """Config lives in the database, not in the instance."""
        first = Project(projectName="CONFIG_PERSIST", filesDirectory=unit_files_directory)
        first.setConfig(alpha=1)

        second = Project(projectName="CONFIG_PERSIST", filesDirectory=unit_files_directory)
        assert second.getConfig()["alpha"] == 1

    def test_config_is_per_project(self, unit_files_directory):
        one = Project(projectName="CONFIG_A", filesDirectory=unit_files_directory)
        two = Project(projectName="CONFIG_B", filesDirectory=unit_files_directory)

        one.setConfig(alpha=1)
        assert "alpha" not in two.getConfig()


@pytest.mark.unit
class TestDefaultProjectRefusesConfig:
    """The default project has no config document, and says so clearly."""

    @pytest.fixture()
    def default_project(self, unit_files_directory):
        return Project(projectName="defaultProject", filesDirectory=unit_files_directory)

    def test_get_config_explains_itself(self, default_project):
        with pytest.raises(ValueError, match="Default project cannot use configuration"):
            default_project.getConfig()

    def test_init_config_explains_itself(self, default_project):
        with pytest.raises(ValueError, match="Default project cannot use configuration"):
            default_project.initConfig(alpha=1)

    def test_set_config_explains_how_to_fix_it(self, default_project):
        """The long message names the two ways out; it must actually be raised."""
        with pytest.raises(ValueError, match="caseConfiguration.json"):
            default_project.setConfig(alpha=1)


@pytest.mark.unit
class TestCounters:
    """Documented contract: the first call DEFINES the counter at 0 and adds
    nothing; every later call adds and returns the new value.  So N calls with
    the default addition leave the counter at N-1, which is surprising enough
    to be worth pinning explicitly."""

    def test_an_unset_counter_reads_as_none(self, unit_project):
        assert unit_project.getCounter("runs") is None

    def test_the_first_call_returns_zero_and_adds_nothing(self, unit_project):
        assert unit_project.getCounterAndAdd("runs") == 0
        assert unit_project.getCounter("runs") == 0

    def test_the_second_call_returns_the_incremented_value(self, unit_project):
        unit_project.getCounterAndAdd("runs")
        assert unit_project.getCounterAndAdd("runs") == 1
        assert unit_project.getCounter("runs") == 1

    def test_increments_accumulate_after_the_definition(self, unit_project):
        for _ in range(4):
            unit_project.getCounterAndAdd("runs")
        assert unit_project.getCounter("runs") == 3

    def test_a_custom_addition_applies_from_the_second_call(self, unit_project):
        assert unit_project.getCounterAndAdd("runs", 5) == 0
        assert unit_project.getCounterAndAdd("runs", 5) == 5
        assert unit_project.getCounter("runs") == 5

    def test_counters_are_independent(self, unit_project):
        unit_project.getCounterAndAdd("a")
        unit_project.getCounterAndAdd("a")
        unit_project.getCounterAndAdd("b", 7)
        assert unit_project.getCounter("a") == 1
        assert unit_project.getCounter("b") == 0

    def test_dots_in_a_name_are_normalised_to_underscores(self, unit_project):
        """MongoDB keys cannot contain dots, so the name is rewritten."""
        unit_project.getCounterAndAdd("group.runs")
        assert unit_project.getCounter("group_runs") == 0

    def test_a_dotted_name_and_its_normalised_form_are_one_counter(self, unit_project):
        unit_project.getCounterAndAdd("group.runs")
        unit_project.getCounterAndAdd("group_runs")
        assert unit_project.getCounter("group.runs") == 1

    @pytest.mark.parametrize("name", ["a__b", "a..b", "a._b", "a_.b"])
    def test_names_that_would_collide_with_mongo_are_rejected(self, unit_project, name):
        with pytest.raises(RuntimeError, match="cannot contain"):
            unit_project.getCounter(name)

    def test_the_deprecated_alias_delegates(self, unit_project):
        """addCounter is documented as deprecated in favour of getCounterAndAdd.

        It must behave identically, first-call semantics included.
        """
        assert unit_project.addCounter("runs") == 0
        assert unit_project.addCounter("runs") == 1
        assert unit_project.getCounter("runs") == 1


@pytest.mark.unit
class TestAddDocumentFromDict:
    """Used by project import; the _cls field selects the collection."""

    def test_a_measurements_document_lands_in_measurements(self, unit_project):
        unit_project.addDocumentFromDict(
            {
                "_cls": "Metadata.Measurements",
                "desc": {"station": "A"},
                "type": "Imported",
                "resource": "/data/a.parquet",
                "dataFormat": "parquet",
            }
        )
        found = unit_project.getMeasurementsDocuments(type="Imported")
        assert len(found) == 1
        assert found[0].resource == "/data/a.parquet"

    def test_a_simulations_document_lands_in_simulations(self, unit_project):
        unit_project.addDocumentFromDict(
            {
                "_cls": "Metadata.Simulations",
                "desc": {},
                "type": "Imported",
                "resource": "/sim/1",
                "dataFormat": "string",
            }
        )
        assert len(unit_project.getSimulationsDocuments(type="Imported")) == 1
        assert len(unit_project.getMeasurementsDocuments(type="Imported")) == 0

    def test_a_cache_document_lands_in_cache(self, unit_project):
        unit_project.addDocumentFromDict(
            {
                "_cls": "Metadata.Cache",
                "desc": {},
                "type": "Imported",
                "resource": "/cache/1",
                "dataFormat": "string",
            }
        )
        assert len(unit_project.getCacheDocuments(type="Imported")) == 1

    def test_an_incoming_project_name_is_dropped(self, unit_project):
        """The document must join THIS project, not the one it came from."""
        unit_project.addDocumentFromDict(
            {
                "_cls": "Metadata.Measurements",
                "projectName": "SOME_OTHER_PROJECT",
                "desc": {},
                "type": "Imported",
                "resource": "/data/a.parquet",
                "dataFormat": "parquet",
            }
        )
        found = unit_project.getMeasurementsDocuments(type="Imported")
        assert len(found) == 1
        assert found[0].projectName == unit_project.projectName

    def test_an_unknown_collection_raises(self, unit_project):
        with pytest.raises(AttributeError):
            unit_project.addDocumentFromDict(
                {
                    "_cls": "Metadata.NoSuchCollection",
                    "desc": {},
                    "type": "Imported",
                    "resource": "/x",
                    "dataFormat": "string",
                }
            )


@pytest.mark.unit
class TestCreateProjectDirectory:
    def test_writes_a_case_configuration_naming_the_project(self, tmp_path):
        target = tmp_path / "case"
        createProjectDirectory(str(target), projectName="MY_PROJECT")

        written = json.loads((target / "caseConfiguration.json").read_text(encoding="utf-8"))
        assert written == {"projectName": "MY_PROJECT"}

    def test_creates_the_directory_if_it_is_missing(self, tmp_path):
        target = tmp_path / "deep" / "case"
        createProjectDirectory(str(target), projectName="MY_PROJECT")
        assert target.is_dir()

    def test_refuses_to_write_for_the_default_project(self, tmp_path):
        """A caseConfiguration naming defaultProject would be a trap."""
        target = tmp_path / "case"
        createProjectDirectory(str(target), projectName="defaultProject")
        assert not (target / "caseConfiguration.json").exists()


@pytest.mark.unit
class TestProjectNameFromDirectory:
    """With projectName=None, the name comes from caseConfiguration.json."""

    def test_name_is_read_from_the_configuration_file(self, tmp_path, unit_files_directory):
        createProjectDirectory(str(tmp_path), projectName="FROM_DISK")
        project = Project(
            configurationPath=str(tmp_path), filesDirectory=unit_files_directory
        )
        assert project.projectName == "FROM_DISK"

    def test_a_missing_file_falls_back_to_the_default_project(
        self, tmp_path, unit_files_directory
    ):
        project = Project(
            configurationPath=str(tmp_path), filesDirectory=unit_files_directory
        )
        assert project.projectName == project.DEFAULTPROJECT

    def test_a_file_without_the_key_is_rejected_with_an_example(
        self, tmp_path, unit_files_directory
    ):
        (tmp_path / "caseConfiguration.json").write_text('{"other": 1}', encoding="utf-8")
        with pytest.raises(ValueError, match="projectName"):
            Project(configurationPath=str(tmp_path), filesDirectory=unit_files_directory)
