"""The project entrypoints: the list, one project, toolkits, create, delete, notebooks."""
import json
import sys
import types

import pytest
from fastapi import HTTPException

import dependencies
from project_api.project_api import DEFAULT_PROJECT, ProjectApi
from project_api.project_files import ProjectFiles
from project_api.toolkit_documents import toolkit_documents
from routers import projects as routes


class FakeDatabase:
    """The little bit of hera the project entrypoints touch, in memory."""

    def __init__(self, projects, documents=None, filesDirectories=None):
        self.projects = list(projects)
        self.documents = {name: list(docs) for name, docs in (documents or {}).items()}
        self.filesDirectories = dict(filesDirectories or {})
        self.created = []
        self.deletedIds = []

    # --- hera.datalayer.All ---

    def getDocumentsAsDict(self, projectName, with_id=False):
        return {"documents": self.documents.get(projectName, [])}

    def deleteDocumentByID(self, docId):
        self.deletedIds.append(docId)
        for name, docs in self.documents.items():
            gone = [d for d in docs if self._idOf(d) == docId]
            self.documents[name] = [d for d in docs if self._idOf(d) != docId]
            # Losing the config document is what makes a project stop existing.
            if any(d.get("type") == f"{name}__config__" for d in gone) and name in self.projects:
                self.projects.remove(name)

    @staticmethod
    def _idOf(doc):
        return doc.get("_id", {}).get("$oid")

    # --- hera.datalayer.project ---

    def getProjectList(self):
        return list(self.projects)

    def Project(self, projectName, filesDirectory=None):
        if filesDirectory:
            self.filesDirectories[projectName] = filesDirectory
        return types.SimpleNamespace(
            filesDirectory=self.filesDirectories.get(projectName, ""))

    # --- hera.datalayer.Cache_Collection ---

    def Cache_Collection(self):
        return types.SimpleNamespace(
            getDocuments=self._notebookDocuments,
            addDocument=self._addDocument)

    def _notebookDocuments(self, projectName, type):
        docs = [d for d in self.documents.get(projectName, []) if d.get("type") == type]
        return [types.SimpleNamespace(resource=d.get("resource")) for d in docs]

    def _addDocument(self, projectName, resource, dataFormat, type, desc):
        doc = {
            "_id": {"$oid": f"id{len(self.documents.get(projectName, []))}"},
            "resource": resource,
            "dataFormat": dataFormat,
            "type": type,
            "desc": desc,
        }
        self.documents.setdefault(projectName, []).append(doc)

    # --- hera.utils.data.CLI ---

    def project_create(self, args):
        self.created.append(args)
        self.projects.append(args.projectName)
        self.filesDirectories[args.projectName] = args.directory


@pytest.fixture
def fake_hera(monkeypatch):
    """Install the fake hera and return the database behind it."""
    def _install(projects, documents=None, filesDirectories=None):
        db = FakeDatabase(projects, documents, filesDirectories)

        datalayer = types.ModuleType("hera.datalayer")
        datalayer.All = types.SimpleNamespace(
            getDocumentsAsDict=db.getDocumentsAsDict,
            deleteDocumentByID=db.deleteDocumentByID)
        datalayer.Cache_Collection = db.Cache_Collection

        project_module = types.ModuleType("hera.datalayer.project")
        project_module.getProjectList = db.getProjectList
        project_module.Project = db.Project

        cli_module = types.ModuleType("hera.utils.data.CLI")
        cli_module.project_create = db.project_create

        hera = types.ModuleType("hera")
        hera.__version__ = "9.9.9"
        hera.toolkitHome = types.SimpleNamespace(
            getToolkitDocuments=lambda: [{"toolkit": "LSM", "desc": {"classpath": "lsm.cls"}}])
        hera.datalayer = datalayer

        monkeypatch.setitem(sys.modules, "hera", hera)
        monkeypatch.setitem(sys.modules, "hera.datalayer", datalayer)
        monkeypatch.setitem(sys.modules, "hera.datalayer.project", project_module)
        monkeypatch.setitem(sys.modules, "hera.utils", types.ModuleType("hera.utils"))
        monkeypatch.setitem(sys.modules, "hera.utils.data", types.ModuleType("hera.utils.data"))
        monkeypatch.setitem(sys.modules, "hera.utils.data.CLI", cli_module)
        return db

    return _install


@pytest.fixture
def warmed(monkeypatch):
    monkeypatch.setattr(dependencies.warmup, "_ready", True)


def _notebook(path):
    return {"_id": {"$oid": "nb"}, "type": "notebook", "resource": path}


class TestProjectList:
    def test_default_project_comes_first(self, fake_hera):
        fake_hera(["zeta", DEFAULT_PROJECT, "Alpha"])

        assert [p["name"] for p in ProjectApi.names()] == [DEFAULT_PROJECT, "Alpha", "zeta"]

    def test_no_projects_is_an_empty_list(self, fake_hera):
        fake_hera([])

        assert ProjectApi.names() == []

    def test_route_answers_the_list(self, fake_hera, warmed):
        fake_hera(["Alpha", DEFAULT_PROJECT])

        assert routes.names() == {"projects": [{"name": DEFAULT_PROJECT}, {"name": "Alpha"}]}

    def test_exists_only_for_a_project_in_the_list(self, fake_hera):
        fake_hera(["Alpha"])

        assert ProjectApi("Alpha").exists()
        assert not ProjectApi("Beta").exists()


class TestProjectDocuments:
    def test_answers_the_name_and_the_documents(self, fake_hera):
        fake_hera(["Alpha"], {"Alpha": [{"x": 1}]})

        assert ProjectApi("Alpha").documents() == {"name": "Alpha", "documents": [{"x": 1}]}

    def test_route_rejects_an_unknown_project(self, fake_hera, warmed):
        fake_hera(["Alpha"], {"Alpha": []})

        with pytest.raises(HTTPException) as raised:
            routes.documents("Nope")

        assert raised.value.status_code == 404

    def test_route_rejects_an_empty_name(self, fake_hera, warmed):
        fake_hera(["Alpha"], {"Alpha": []})

        with pytest.raises(HTTPException) as raised:
            routes.documents("")

        assert raised.value.status_code == 400

    def test_route_answers_while_warm(self, fake_hera, warmed):
        fake_hera(["Alpha"], {"Alpha": [{"x": 1}]})

        assert routes.documents("Alpha") == {"project": {"name": "Alpha", "documents": [{"x": 1}]}}

    def test_route_says_warming_up(self, fake_hera, monkeypatch):
        fake_hera(["Alpha"], {"Alpha": []})
        monkeypatch.setattr(dependencies.warmup, "_ready", False)

        with pytest.raises(HTTPException) as raised:
            routes.documents("Alpha")

        assert raised.value.status_code == 503
        assert raised.value.detail == "WARMING_UP"


class TestToolkits:
    def test_answers_the_toolkit_documents(self, fake_hera, warmed):
        fake_hera(["Alpha"])

        assert toolkit_documents()[0]["toolkit"] == "LSM"
        assert routes.toolkits()["toolkitDocs"][0]["toolkit"] == "LSM"


class TestDetectNotebooks:
    def test_registers_untracked_notebooks(self, fake_hera, tmp_path):
        db = fake_hera(["Alpha"], {"Alpha": []})
        (tmp_path / "a.ipynb").write_text("{}")
        (tmp_path / "sub").mkdir()
        (tmp_path / "sub" / "b.ipynb").write_text("{}")

        added = ProjectApi("Alpha").detect_notebooks(str(tmp_path))

        assert added == 2
        names = sorted(d["desc"]["datasourceName"] for d in db.documents["Alpha"])
        assert names == ["a", "b"]

    def test_skips_notebooks_already_tracked(self, fake_hera, tmp_path):
        notebook = tmp_path / "a.ipynb"
        notebook.write_text("{}")
        fake_hera(["Alpha"], {"Alpha": [_notebook(str(notebook))]})

        assert ProjectApi("Alpha").detect_notebooks(str(tmp_path)) == 0

    def test_skips_checkpoint_copies(self, fake_hera, tmp_path):
        checkpoints = tmp_path / ".ipynb_checkpoints"
        checkpoints.mkdir()
        (checkpoints / "a-checkpoint.ipynb").write_text("{}")
        fake_hera(["Alpha"], {"Alpha": []})

        assert ProjectApi("Alpha").detect_notebooks(str(tmp_path)) == 0

    def test_route_answers_the_project_and_the_count(self, fake_hera, warmed, tmp_path):
        fake_hera(["Alpha"], {"Alpha": []})
        (tmp_path / "a.ipynb").write_text("{}")

        answer = routes.detect(types.SimpleNamespace(
            projectName="Alpha", filesDirectory=str(tmp_path)))

        assert answer["addedNotebooks"] == 1
        assert answer["project"]["name"] == "Alpha"
        assert len(answer["project"]["documents"]) == 1

    def test_route_rejects_a_project_with_no_files_directory(self, fake_hera, warmed):
        fake_hera(["Alpha"], {"Alpha": []})

        with pytest.raises(HTTPException) as raised:
            routes.detect(types.SimpleNamespace(projectName="Alpha", filesDirectory=""))

        assert raised.value.status_code == 400


class TestCreateProject:
    def test_creates_the_project_in_the_default_directory(self, fake_hera):
        db = fake_hera([])

        project = ProjectApi("NewOne")

        directory = project.create()

        assert directory == project.default_files_directory()
        assert project.files_directory() == directory
        assert db.created[0].projectName == "NewOne"
        assert db.created[0].overwrite is False
        assert db.created[0].loadRepositories is True
        assert db.projects == ["NewOne"]

    def test_uses_the_directory_the_user_gave(self, fake_hera):
        db = fake_hera([])

        directory = ProjectApi("NewOne").create("/data/mine", loadRepositories=False)

        assert directory == "/data/mine"
        assert db.created[0].directory == "/data/mine"
        assert db.created[0].loadRepositories is False

    def test_route_also_detects_the_notebooks_already_there(self, fake_hera, warmed, tmp_path):
        fake_hera([])
        (tmp_path / "old.ipynb").write_text("{}")

        answer = routes.create(types.SimpleNamespace(
            projectName=" NewOne ", filesDirectory=str(tmp_path), loadRepositories=True))

        assert answer["addedNotebooks"] == 1
        # Both, so creating a project costs one call: the list to show it in, and
        # the project itself so the client need not load it again.
        assert {p["name"] for p in answer["projects"]} == {"NewOne"}
        assert answer["project"]["name"] == "NewOne"

    def test_route_rejects_a_blank_name(self, fake_hera, warmed):
        fake_hera([])

        with pytest.raises(HTTPException) as raised:
            routes.create(types.SimpleNamespace(
                projectName="  ", filesDirectory="", loadRepositories=True))

        assert raised.value.status_code == 400


class TestDeleteProject:
    def test_deletes_every_document_config_last(self, fake_hera):
        db = fake_hera(["Alpha"], {"Alpha": [
            {"_id": {"$oid": "cfg"}, "type": "Alpha__config__", "resource": ""},
            {"_id": {"$oid": "one"}, "type": "notebook", "resource": ""},
        ]})

        ProjectApi("Alpha").delete()

        assert db.deletedIds == ["one", "cfg"]

    def test_keeps_the_files_by_default(self, fake_hera, tmp_path):
        data = tmp_path / "a.parquet"
        data.write_text("x")
        fake_hera(["Alpha"], {"Alpha": [
            {"_id": {"$oid": "one"}, "type": "measure", "resource": str(data)},
        ]}, {"Alpha": str(tmp_path)})

        ProjectApi("Alpha").delete(deleteFiles=False)

        assert data.exists()

    def test_deletes_the_files_inside_the_project_directory(self, fake_hera, tmp_path):
        inside = tmp_path / "a.parquet"
        inside.write_text("x")
        folder = tmp_path / "run"
        folder.mkdir()
        (folder / "f").write_text("y")
        fake_hera(["Alpha"], {"Alpha": [
            {"_id": {"$oid": "one"}, "type": "measure", "resource": str(inside)},
            {"_id": {"$oid": "two"}, "type": "simulation", "resource": str(folder)},
        ]}, {"Alpha": str(tmp_path)})

        ProjectApi("Alpha").delete(deleteFiles=True)

        assert not inside.exists()
        assert not folder.exists()
        assert not tmp_path.exists()

    def test_leaves_files_outside_the_project_directory(self, fake_hera, tmp_path):
        outside = tmp_path / "shared.parquet"
        outside.write_text("x")
        projectDir = tmp_path / "proj"
        projectDir.mkdir()
        fake_hera(["Alpha"], {"Alpha": [
            {"_id": {"$oid": "one"}, "type": "measure", "resource": str(outside)},
        ]}, {"Alpha": str(projectDir)})

        ProjectApi("Alpha").delete(deleteFiles=True)

        assert outside.exists()

    def test_removes_the_case_configuration_it_wrote(self, fake_hera, tmp_path):
        (tmp_path / "caseConfiguration.json").write_text(json.dumps({}))
        fake_hera(["Alpha"], {"Alpha": []}, {"Alpha": str(tmp_path)})

        ProjectApi("Alpha").delete(deleteFiles=True)

        assert not tmp_path.exists()

    def test_keeps_a_directory_that_still_holds_other_files(self, fake_hera, tmp_path):
        (tmp_path / "someones.txt").write_text("x")
        fake_hera(["Alpha"], {"Alpha": []}, {"Alpha": str(tmp_path)})

        ProjectApi("Alpha").delete(deleteFiles=True)

        assert tmp_path.is_dir()

    def test_a_project_without_a_files_directory_touches_no_files(self, fake_hera, tmp_path):
        data = tmp_path / "a.parquet"
        data.write_text("x")
        fake_hera(["Alpha"], {"Alpha": [
            {"_id": {"$oid": "one"}, "type": "measure", "resource": str(data)},
        ]})

        ProjectApi("Alpha").delete(deleteFiles=True)

        assert data.exists()

    def test_route_answers_the_remaining_list_and_the_next_project(self, fake_hera, warmed):
        fake_hera(["Alpha", "Beta"],
                  {"Alpha": [{"_id": {"$oid": "cfg"}, "type": "Alpha__config__", "resource": ""}],
                   "Beta": [{"x": 1}]})

        answer = routes.delete(types.SimpleNamespace(projectName="Alpha", deleteFiles=False))

        assert answer["project"]["name"] == "Beta"
        assert {p["name"] for p in answer["projects"]} == {"Beta"}

    def test_route_answers_no_project_when_the_last_one_goes(self, fake_hera, warmed):
        fake_hera(["Alpha"], {"Alpha": [{"_id": {"$oid": "cfg"}, "type": "Alpha__config__", "resource": ""}]})

        answer = routes.delete(types.SimpleNamespace(projectName="Alpha", deleteFiles=False))

        assert answer["projects"] == []
        assert answer["project"] is None

    def test_route_rejects_an_unknown_project(self, fake_hera, warmed):
        fake_hera(["Alpha"], {"Alpha": []})

        with pytest.raises(HTTPException) as raised:
            routes.delete(types.SimpleNamespace(projectName="Nope", deleteFiles=False))

        assert raised.value.status_code == 404

    def test_route_says_so_when_the_project_survives_the_delete(self, fake_hera, warmed, monkeypatch):
        fake_hera(["Alpha"], {"Alpha": []})
        # A delete that quietly does nothing must not look like a success.
        monkeypatch.setattr(ProjectApi, "delete", lambda self, deleteFiles=False: None)

        with pytest.raises(HTTPException) as raised:
            routes.delete(types.SimpleNamespace(projectName="Alpha", deleteFiles=False))

        assert raised.value.status_code == 500
        assert "still in the project list" in raised.value.detail


class TestProjectFiles:
    def test_holds_only_what_is_inside(self, tmp_path):
        files = ProjectFiles(str(tmp_path))

        assert files.holds(str(tmp_path / "a.parquet"))
        assert files.holds(str(tmp_path / "sub" / "a.parquet"))
        assert not files.holds(str(tmp_path.parent / "other.parquet"))
        assert not files.holds("")

    def test_a_project_with_no_directory_holds_nothing(self, tmp_path):
        files = ProjectFiles()

        assert not files.holds(str(tmp_path / "a.parquet"))
        assert not files.remove(str(tmp_path))
        assert not files.remove_if_empty()

    def test_removes_a_file_and_a_folder_inside(self, tmp_path):
        oneFile = tmp_path / "a.parquet"
        oneFile.write_text("x")
        folder = tmp_path / "run"
        folder.mkdir()
        (folder / "f").write_text("y")
        files = ProjectFiles(str(tmp_path))

        assert files.remove(str(oneFile))
        assert files.remove(str(folder))
        assert not oneFile.exists()
        assert not folder.exists()

    def test_leaves_a_file_outside_alone(self, tmp_path):
        outside = tmp_path / "shared.parquet"
        outside.write_text("x")
        inside = tmp_path / "proj"
        inside.mkdir()
        files = ProjectFiles(str(inside))

        assert not files.remove(str(outside))
        assert outside.exists()

    def test_removing_a_path_that_is_already_gone_says_so(self, tmp_path):
        files = ProjectFiles(str(tmp_path))

        assert not files.remove(str(tmp_path / "never-was"))

    def test_drops_the_directory_once_only_leftovers_remain(self, tmp_path):
        (tmp_path / ProjectFiles.CASE_CONFIG).write_text("{}")
        files = ProjectFiles(str(tmp_path))

        assert files.remove_if_empty()
        assert not tmp_path.exists()

    def test_keeps_a_directory_holding_someone_elses_file(self, tmp_path):
        (tmp_path / "someones.txt").write_text("x")
        files = ProjectFiles(str(tmp_path))

        assert not files.remove_if_empty()
        assert tmp_path.is_dir()


# Every route reaches hera, so each one refuses while hera is still warming.
WARMING_CALLS = {
    "names": lambda: routes.names(),
    "documents": lambda: routes.documents("Alpha"),
    "toolkits": lambda: routes.toolkits(),
    "create": lambda: routes.create(types.SimpleNamespace(
        projectName="Alpha", filesDirectory="", loadRepositories=True)),
    "delete": lambda: routes.delete(types.SimpleNamespace(
        projectName="Alpha", deleteFiles=False)),
    "detectNotebooks": lambda: routes.detect(types.SimpleNamespace(
        projectName="Alpha", filesDirectory="/tmp")),
}


@pytest.mark.parametrize("route", sorted(WARMING_CALLS))
def test_every_route_says_warming_up(route, fake_hera, monkeypatch):
    fake_hera(["Alpha"], {"Alpha": []})
    monkeypatch.setattr(dependencies.warmup, "_ready", False)

    with pytest.raises(HTTPException) as raised:
        WARMING_CALLS[route]()

    assert raised.value.status_code == 503
    assert raised.value.detail == "WARMING_UP"
