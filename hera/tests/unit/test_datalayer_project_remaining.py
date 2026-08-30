"""Project: the query/save helpers batch 2 didn't reach.

B88: Project.getProjectList (the classmethod) is decorated @staticmethod
but declares an unused `cls` leading parameter, so the intended call
`Project.getProjectList()` raises TypeError. See TestClassmethodGetProjectListIsBroken.
"""
import pandas
import pytest

from hera.datalayer import Project
from hera.datalayer.project import getProjectList


@pytest.mark.unit
class TestFilesDirectoryAliases:
    def test_filesdirectory_and_filesDirectory_agree(self, unit_project):
        assert unit_project.FilesDirectory == unit_project.filesDirectory


@pytest.mark.unit
class TestUpdateProjectNameOnDoc:
    def test_it_stamps_the_projects_own_name(self, unit_project):
        result = unit_project.updateProjectNameOnDoc({"a": 1})
        assert result == {"a": 1, "projectName": unit_project.projectName}

    def test_it_overwrites_an_existing_projectName_key(self, unit_project):
        result = unit_project.updateProjectNameOnDoc({"projectName": "someone else"})
        assert result["projectName"] == unit_project.projectName


@pytest.mark.unit
class TestSaveDataFamily:
    def test_save_measurement_data_registers_a_measurements_document(self, unit_project):
        doc = unit_project.saveMeasurementData("myData", pandas.DataFrame({"a": [1, 2]}), desc={"k": "v"})
        assert doc.type == "myData"
        assert doc.desc["k"] == "v"

    def test_save_cache_data_registers_a_cache_document(self, unit_project):
        doc = unit_project.saveCacheData("myCache", pandas.DataFrame({"a": [1]}), desc={})
        found = unit_project.getCacheDocuments(type="myCache")
        assert len(found) == 1

    def test_save_simulation_data_registers_a_simulations_document(self, unit_project):
        doc = unit_project.saveSimulationData("mySim", pandas.DataFrame({"a": [1]}), desc={})
        found = unit_project.getSimulationsDocuments(type="mySim")
        assert len(found) == 1

    def test_an_explicit_type_overrides_the_name(self, unit_project):
        doc = unit_project.saveMeasurementData("myData", pandas.DataFrame({"a": [1]}), desc={}, type="explicitType")
        assert doc.type == "explicitType"

    def test_the_saved_document_can_be_read_back_by_id(self, unit_project):
        doc = unit_project.saveMeasurementData("myData", pandas.DataFrame({"a": [1]}), desc={})
        assert unit_project.getDocumentByID(doc.id).id == doc.id


@pytest.mark.unit
class TestQueryHelpers:
    def test_get_metadata_returns_one_row_per_document(self, unit_project):
        unit_project.saveMeasurementData("a", pandas.DataFrame({"x": [1]}), desc={"k": 1})
        unit_project.saveMeasurementData("b", pandas.DataFrame({"x": [2]}), desc={"k": 2})
        metadata = unit_project.getMetadata()
        assert len(metadata) >= 2

    def test_get_all_documents_spans_every_collection(self, unit_project):
        unit_project.saveMeasurementData("m", pandas.DataFrame({"x": [1]}), desc={})
        unit_project.saveCacheData("c", pandas.DataFrame({"x": [1]}), desc={})
        unit_project.saveSimulationData("s", pandas.DataFrame({"x": [1]}), desc={})
        all_docs = unit_project.getAllDocuments(type="m")
        assert len(all_docs) == 1

    def test_get_measurements_documents_as_dict_matches_the_query(self, unit_project):
        unit_project.saveMeasurementData("m", pandas.DataFrame({"x": [1]}), desc={})
        result = unit_project.getMeasurementsDocumentsAsDict(type="m")
        assert len(result["documents"]) == 1
        assert result["documents"][0]["type"] == "m"

    def test_get_simulations_documents_as_dict_matches_the_query(self, unit_project):
        unit_project.saveSimulationData("s", pandas.DataFrame({"x": [1]}), desc={})
        result = unit_project.getSimulationsDocumentsAsDict(type="s")
        assert len(result["documents"]) == 1

    def test_get_cache_documents_as_dict_matches_the_query(self, unit_project):
        unit_project.saveCacheData("c", pandas.DataFrame({"x": [1]}), desc={})
        result = unit_project.getCacheDocumentsAsDict(type="c")
        assert len(result["documents"]) == 1

    def test_delete_measurements_documents_removes_matching_docs(self, unit_project):
        unit_project.saveMeasurementData("m", pandas.DataFrame({"x": [1]}), desc={})
        deleted = unit_project.deleteMeasurementsDocuments(type="m")
        assert len(deleted) == 1
        # a mongoengine QuerySet never equals a plain list, even empty --
        # must be materialised first.
        assert list(unit_project.getMeasurementsDocuments(type="m")) == []

    def test_delete_simulations_documents_removes_matching_docs(self, unit_project):
        unit_project.saveSimulationData("s", pandas.DataFrame({"x": [1]}), desc={})
        deleted = unit_project.deleteSimulationsDocuments(type="s")
        assert len(deleted) == 1

    def test_delete_cache_documents_removes_matching_docs(self, unit_project):
        unit_project.saveCacheData("c", pandas.DataFrame({"x": [1]}), desc={})
        deleted = unit_project.deleteCacheDocuments(type="c")
        assert len(deleted) == 1


@pytest.mark.unit
class TestSetCounter:
    def test_get_counter_and_add_returns_the_value_after_adding(self, unit_project):
        """getCounterAndAdd's docstring reads as "return the value, then
        add" but it actually increments first and returns the new value --
        so after setCounter(defaultValue=5) the first read is 6, not 5."""
        unit_project.setCounter("myCounter", defaultValue=5)
        assert unit_project.getCounterAndAdd("myCounter") == 6

    def test_the_counter_increments_on_each_read(self, unit_project):
        unit_project.setCounter("myCounter", defaultValue=0)
        first = unit_project.getCounterAndAdd("myCounter")
        second = unit_project.getCounterAndAdd("myCounter")
        assert second == first + 1


@pytest.mark.unit
class TestGetProjectList:
    def test_the_current_project_appears_in_the_module_level_list(self, unit_project):
        assert unit_project.projectName in getProjectList()


@pytest.mark.unit
class TestClassmethodGetProjectListIsBroken:
    """B88: Project.getProjectList is decorated @staticmethod but its
    signature is (cls, user=None), with `cls` never used in the body.
    Calling it the intended way -- Project.getProjectList() -- raises
    TypeError for a missing positional argument."""

    @pytest.mark.xfail(
        strict=True,
        reason="B88: decorated @staticmethod but declares an unused `cls` "
               "leading positional parameter -- Project.getProjectList() "
               "raises TypeError instead of returning the project list. "
               "See the consolidated findings issue.",
    )
    def test_it_should_be_callable_with_no_arguments(self, unit_project):
        assert unit_project.projectName in Project.getProjectList()

    def test_it_currently_raises_for_the_intended_call(self):
        """Characterisation of B88."""
        with pytest.raises(TypeError, match="cls"):
            Project.getProjectList()

    def test_it_only_works_if_the_caller_smuggles_something_into_cls(self, unit_project):
        """Characterisation of B88: the unused `cls` slot can be filled
        with anything (even None) and the call proceeds normally."""
        assert unit_project.projectName in Project.getProjectList(None)
