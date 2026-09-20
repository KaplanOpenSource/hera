"""ToolkitRepository: registering/looking-up dynamic toolkits as
ToolkitDataSource measurements documents. Backed by the session-wide
mongomock reroute, so a plain Project(projectName=...) here already talks
to the in-memory database -- no fixture needed beyond a unique project name
per test to avoid cross-test document leakage."""
import os

import pytest

from hera.datalayer import Project
from hera.utils.data.toolkit_repository import ToolkitRepository


class _DummyToolkit:
    pass


@pytest.fixture(autouse=True)
def _toolkit_repository_uses_a_temp_files_directory(monkeypatch, unit_files_directory):
    """ToolkitRepository always builds its own Project(projectName=...) with
    no way to pass filesDirectory -- redirect it to the per-test directory
    so nothing lands in the guarded fake home."""
    original = Project

    def _project_with_files_directory(projectName=None, **kwargs):
        kwargs.setdefault("filesDirectory", unit_files_directory)
        return original(projectName=projectName, **kwargs)

    monkeypatch.setattr(
        "hera.utils.data.toolkit_repository.Project", _project_with_files_directory
    )


@pytest.mark.unit
class TestRegisterToolkitViaClass:
    def test_it_derives_classpath_from_the_class(self):
        repo = ToolkitRepository(project_name="TK_REPO_1")
        doc = repo.registerToolkit(toolkitclass=_DummyToolkit, datasource_name="myTool")
        assert doc.desc["classpath"].endswith("_DummyToolkit")

    def test_it_derives_resource_as_the_containing_directory(self):
        repo = ToolkitRepository(project_name="TK_REPO_2")
        doc = repo.registerToolkit(toolkitclass=_DummyToolkit, datasource_name="myTool")
        assert doc.resource == os.path.dirname(__file__)

    def test_it_requires_either_a_class_or_resource_and_classpath(self):
        repo = ToolkitRepository(project_name="TK_REPO_3")
        with pytest.raises(TypeError, match="registerToolkit requires"):
            repo.registerToolkit(datasource_name="myTool")


@pytest.mark.unit
class TestRegisterToolkitViaExplicitPath:
    def test_it_accepts_resource_and_classpath_directly(self):
        repo = ToolkitRepository(project_name="TK_REPO_4")
        doc = repo.registerToolkit(
            datasource_name="myTool", resource="/some/dir", classpath="pkg.mod.Cls"
        )
        assert doc.desc["classpath"] == "pkg.mod.Cls"
        assert doc.resource == "/some/dir"

    def test_the_default_version_is_0_0_1(self):
        repo = ToolkitRepository(project_name="TK_REPO_5")
        doc = repo.registerToolkit(
            datasource_name="myTool", resource="/some/dir", classpath="pkg.mod.Cls"
        )
        assert doc.desc["version"] == [0, 0, 1]

    def test_an_explicit_version_tuple_is_stored_as_a_list(self):
        repo = ToolkitRepository(project_name="TK_REPO_6")
        doc = repo.registerToolkit(
            datasource_name="myTool", resource="/some/dir", classpath="pkg.mod.Cls",
            version=(1, 2, 3),
        )
        assert doc.desc["version"] == [1, 2, 3]


@pytest.mark.unit
class TestRegisterToolkitOverwrite:
    def test_overwrite_removes_the_previous_registration(self):
        repo = ToolkitRepository(project_name="TK_REPO_7")
        repo.registerToolkit(
            datasource_name="myTool", resource="/v1", classpath="pkg.mod.Cls"
        )
        repo.registerToolkit(
            datasource_name="myTool", resource="/v2", classpath="pkg.mod.Cls", overwrite=True
        )
        table = repo.getToolkitTable()
        assert len(table) == 1
        assert table.iloc[0]["cls"] == "pkg.mod.Cls"

    def test_without_overwrite_both_registrations_are_kept(self):
        repo = ToolkitRepository(project_name="TK_REPO_8")
        repo.registerToolkit(
            datasource_name="myTool", resource="/v1", classpath="pkg.mod.Cls"
        )
        repo.registerToolkit(
            datasource_name="myTool", resource="/v2", classpath="pkg.mod.Cls"
        )
        assert len(repo.getToolkitTable()) == 2


@pytest.mark.unit
class TestGetToolkitDocument:
    def test_it_finds_a_registered_toolkit_by_datasource_name(self):
        repo = ToolkitRepository(project_name="TK_REPO_9")
        repo.registerToolkit(datasource_name="myTool", resource="/v1", classpath="pkg.mod.Cls")
        found = repo.getToolkitDocument("myTool")
        assert found is not None
        assert found.desc["classpath"] == "pkg.mod.Cls"

    def test_it_returns_none_for_an_unregistered_name(self):
        repo = ToolkitRepository(project_name="TK_REPO_10")
        assert repo.getToolkitDocument("noSuchTool") is None


@pytest.mark.unit
class TestGetToolkitTable:
    def test_an_empty_project_yields_an_empty_table(self):
        repo = ToolkitRepository(project_name="TK_REPO_11")
        table = repo.getToolkitTable()
        assert len(table) == 0

    def test_each_registration_becomes_one_row(self):
        repo = ToolkitRepository(project_name="TK_REPO_12")
        repo.registerToolkit(datasource_name="toolA", resource="/a", classpath="pkg.A")
        repo.registerToolkit(datasource_name="toolB", resource="/b", classpath="pkg.B")
        table = repo.getToolkitTable()
        assert sorted(table["toolkit"]) == ["toolA", "toolB"]
