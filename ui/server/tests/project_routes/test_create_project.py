"""Creating a project."""
import types

import pytest
from fastapi import HTTPException

from project_api.project_api import ProjectApi
from routers import projects as routes


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
