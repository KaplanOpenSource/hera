"""LSMToolkit: construction/property validation, template loading, and
getTemplatesTable.

getSimulations/getSimulationsList/prepareSlurmLSMExecution are left
uncovered: they need real SingleSimulation documents or a working Slurm
environment, out of scope for a pure unit test.
"""
import json

import pytest

from hera import toolkit as toolkit_mod
from hera.simulations.LSM.template import LSMTemplate


@pytest.fixture()
def lsm(unit_toolkit_factory):
    from hera import toolkitHome

    return unit_toolkit_factory(toolkitHome.LSM)


@pytest.fixture()
def template_file(tmp_path):
    path = tmp_path / "tpl.json"
    path.write_text(json.dumps({"name": "myTemplate", "params": {"a": 1}}))
    return str(path)


@pytest.mark.unit
class TestConstruction:
    def test_defaults_are_stored(self, lsm):
        assert lsm.to_xarray is True
        assert lsm.to_database is False
        assert lsm.forceKeep is False

    def test_analysis_is_wired_up(self, lsm):
        assert lsm.analysis is not None
        assert lsm.analysis.datalayer is lsm

    def test_single_simulation_exposes_the_class(self, lsm):
        from hera.simulations.LSM.singleSimulation import SingleSimulation

        assert lsm.singleSimulation is SingleSimulation


@pytest.mark.unit
class TestPropertyValidation:
    def test_to_xarray_rejects_non_boolean(self, lsm):
        with pytest.raises(ValueError, match="to_xarray must be boolean"):
            lsm.to_xarray = "yes"

    def test_to_database_rejects_non_boolean(self, lsm):
        with pytest.raises(ValueError, match="to_xarray must be boolean"):
            lsm.to_database = "yes"

    def test_force_keep_rejects_non_boolean(self, lsm):
        with pytest.raises(ValueError, match="to_xarray must be boolean"):
            lsm.forceKeep = "yes"

    def test_valid_boolean_values_are_accepted(self, lsm):
        lsm.to_xarray = False
        assert lsm.to_xarray is False


@pytest.mark.unit
class TestAnalysisWiring:
    def test_coordinate_handler_is_constructed(self, lsm):
        assert lsm.analysis.coordinateHandler is not None


@pytest.mark.unit
class TestLoadDataNoSave:
    def test_it_builds_a_template_without_touching_the_db(self, lsm, template_file):
        tpl = lsm.loadData(template_file, saveMode=toolkit_mod.TOOLKIT_SAVEMODE_NOSAVE)
        assert isinstance(tpl, LSMTemplate)
        assert lsm.getTemplates() == []

    def test_a_non_file_string_raises(self, lsm):
        with pytest.raises(ValueError, match="must be a JSON template file"):
            lsm.loadData("/no/such/file.json")

    def test_a_non_string_argument_raises(self, lsm):
        with pytest.raises(ValueError, match="must be a JSON template file"):
            lsm.loadData(123)


@pytest.mark.unit
class TestLoadDataFileAndDB:
    def test_it_registers_a_findable_datasource(self, lsm, template_file):
        lsm.loadData(template_file)
        templates = lsm.getTemplates()
        assert len(templates) == 1
        assert isinstance(templates[0], LSMTemplate)

    def test_get_template_by_name_finds_it(self, lsm, template_file):
        lsm.loadData(template_file)
        found = lsm.getTemplateByName("myTemplate")
        assert found is not None

    def test_get_template_by_name_returns_none_for_an_unknown_name(self, lsm):
        assert lsm.getTemplateByName("noSuchTemplate") is None

    def test_loading_the_same_template_twice_raises(self, lsm, template_file):
        lsm.loadData(template_file)
        with pytest.raises(ValueError, match="already exists in the DB"):
            lsm.loadData(template_file)


@pytest.mark.unit
class TestGetTemplatesTable:
    def test_an_empty_project_yields_an_empty_table(self, lsm):
        table = lsm.getTemplatesTable()
        assert len(table) == 0

    def test_it_should_list_a_registered_template(self, lsm, template_file):
        lsm.loadData(template_file)
        table = lsm.getTemplatesTable()
        assert len(table) == 1

    def test_the_template_json_reaches_the_stored_document(self, lsm, template_file):
        """loadData must persist the parsed template JSON, not only the
        caller's own kwargs -- everything downstream reads it from there."""
        lsm.loadData(template_file)
        stored = lsm.getDataSourceDocument("myTemplate")
        assert stored.desc["params"] == {"a": 1}
        assert lsm.getTemplateByName("myTemplate").params == {"a": 1}
