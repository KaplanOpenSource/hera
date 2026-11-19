"""Deleting a project and its files."""
import json
import types

import pytest
from fastapi import HTTPException

from project_api.project_api import ProjectApi
from routers import projects as routes


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
