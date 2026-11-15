"""Project.rename and Project.relocate (#962).

rename moves every document of the project to the new name, the config
document included (its type carries the project name).  relocate moves the
files directory and rewrites the resources that point into it.

Both run the real Project against the in-memory store.  The DB is inspected
through AbstractCollection rather than by reopening the old project, because
opening a project creates its config document.
"""
import os

import pytest

from hera.datalayer import Project
from hera.datalayer.collection import AbstractCollection


def _documentsOf(projectName):
    return list(AbstractCollection().getDocuments(projectName=projectName))


def _addOneOfEach(project):
    project.addMeasurementsDocument(resource="m", type="M")
    project.addSimulationsDocument(resource="s", type="S")
    project.addCacheDocument(resource="c", type="C")


@pytest.mark.unit
class TestRename:
    def test_every_document_moves_to_the_new_name(self, unit_project):
        _addOneOfEach(unit_project)
        before = len(_documentsOf(unit_project.projectName))
        oldName = unit_project.projectName

        unit_project.rename("RENAMED")

        assert _documentsOf(oldName) == []
        assert len(_documentsOf("RENAMED")) == before

    def test_the_config_and_counters_follow_the_project(self, unit_project, unit_files_directory):
        unit_project.setConfig(answer=42)
        unit_project.defineCounter("runs", 7)

        unit_project.rename("RENAMED")

        reopened = Project(projectName="RENAMED")
        assert reopened.getConfig()["answer"] == 42
        assert reopened.getConfig()["filesDirectory"] == unit_files_directory
        assert reopened.getCounter("runs") == 7

    def test_the_config_document_type_carries_the_new_name(self, unit_project):
        unit_project.setConfig(answer=42)

        unit_project.rename("RENAMED")

        assert [doc.type for doc in _documentsOf("RENAMED")] == ["RENAMED__config__"]

    def test_the_instance_follows_the_new_name(self, unit_project):
        unit_project.addMeasurementsDocument(resource="m", type="M")

        unit_project.rename("RENAMED")

        assert unit_project.projectName == "RENAMED"
        assert [doc.resource for doc in unit_project.getMeasurementsDocuments(type="M")] == ["m"]

    def test_other_projects_are_left_alone(self, unit_project, unit_files_directory):
        other = Project(projectName="OTHER", filesDirectory=unit_files_directory)
        other.addMeasurementsDocument(resource="o", type="M")
        before = len(_documentsOf("OTHER"))

        unit_project.rename("RENAMED")

        assert len(_documentsOf("OTHER")) == before

    def test_a_name_that_is_already_taken_is_refused_and_nothing_moves(
            self, unit_project, unit_files_directory):
        Project(projectName="TAKEN", filesDirectory=unit_files_directory)
        oldName = unit_project.projectName
        before = len(_documentsOf(oldName))

        with pytest.raises(ValueError, match="TAKEN"):
            unit_project.rename("TAKEN")

        assert len(_documentsOf(oldName)) == before
        assert unit_project.projectName == oldName

    @pytest.mark.parametrize("newName", ["", Project.DEFAULTPROJECT])
    def test_an_empty_or_default_new_name_is_refused(self, unit_project, newName):
        with pytest.raises(ValueError):
            unit_project.rename(newName)

    def test_renaming_to_the_same_name_is_refused(self, unit_project):
        with pytest.raises(ValueError):
            unit_project.rename(unit_project.projectName)

    def test_the_default_project_cannot_be_renamed(self, unit_files_directory):
        default = Project(projectName=Project.DEFAULTPROJECT, filesDirectory=unit_files_directory)
        with pytest.raises(ValueError):
            default.rename("RENAMED")


@pytest.mark.unit
class TestRelocate:
    @pytest.fixture()
    def target(self, tmp_path):
        return str(tmp_path / "moved" / "heraFiles")

    def test_the_files_move_to_the_new_directory(self, unit_project, unit_files_directory, target):
        with open(os.path.join(unit_files_directory, "data.txt"), "w") as f:
            f.write("x")

        unit_project.relocate(target)

        assert not os.path.exists(unit_files_directory)
        with open(os.path.join(target, "data.txt")) as f:
            assert f.read() == "x"

    def test_resources_inside_the_old_directory_point_into_the_new_one(
            self, unit_project, unit_files_directory, target):
        unit_project.addMeasurementsDocument(
            resource=os.path.join(unit_files_directory, "sub", "a.parquet"), type="M")
        unit_project.addCacheDocument(resource=unit_files_directory, type="C")

        unit_project.relocate(target)

        assert unit_project.getMeasurementsDocuments(type="M")[0].resource == \
            os.path.join(target, "sub", "a.parquet")
        assert unit_project.getCacheDocuments(type="C")[0].resource == target

    def test_resources_outside_the_old_directory_are_left_alone(
            self, unit_project, unit_files_directory, target):
        """Including a sibling whose name merely starts with the old directory's."""
        sibling = unit_files_directory + "2/a.parquet"
        unit_project.addMeasurementsDocument(resource="/elsewhere/a.parquet", type="OUT")
        unit_project.addMeasurementsDocument(resource=sibling, type="SIBLING")
        unit_project.addMeasurementsDocument(resource={"not": "a path"}, type="DICT", dataFormat="JSON_dict")

        unit_project.relocate(target)

        assert unit_project.getMeasurementsDocuments(type="OUT")[0].resource == "/elsewhere/a.parquet"
        assert unit_project.getMeasurementsDocuments(type="SIBLING")[0].resource == sibling
        assert unit_project.getMeasurementsDocuments(type="DICT")[0].resource == {"not": "a path"}

    def test_other_projects_resources_are_left_alone(self, unit_project, unit_files_directory, target):
        other = Project(projectName="OTHER", filesDirectory=unit_files_directory)
        inside = os.path.join(unit_files_directory, "a.parquet")
        other.addMeasurementsDocument(resource=inside, type="M")

        unit_project.relocate(target)

        assert other.getMeasurementsDocuments(type="M")[0].resource == inside

    def test_the_new_directory_is_saved_in_the_config(self, unit_project, target):
        unit_project.relocate(target)

        assert unit_project.filesDirectory == target
        assert Project(projectName=unit_project.projectName).filesDirectory == target

    def test_an_existing_target_is_used_as_is(self, unit_project, unit_files_directory, target):
        """The files were already copied there (e.g. with rsync): only the
        references change, and neither directory is touched."""
        os.makedirs(target)
        inside = os.path.join(unit_files_directory, "a.parquet")
        unit_project.addMeasurementsDocument(resource=inside, type="M")

        unit_project.relocate(target)

        assert os.path.isdir(unit_files_directory)
        assert unit_project.getMeasurementsDocuments(type="M")[0].resource == os.path.join(target, "a.parquet")

    def test_a_target_inside_the_files_directory_is_refused(self, unit_project, unit_files_directory):
        with pytest.raises(ValueError):
            unit_project.relocate(os.path.join(unit_files_directory, "inner"))

    def test_the_default_project_cannot_be_relocated(self, unit_files_directory, target):
        default = Project(projectName=Project.DEFAULTPROJECT, filesDirectory=unit_files_directory)
        with pytest.raises(ValueError):
            default.relocate(target)
