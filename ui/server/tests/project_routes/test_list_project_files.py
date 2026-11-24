"""Listing the files in a project's folder."""
import pytest
from fastapi import HTTPException

from project_api.project_files import ProjectFiles
from routers import projects as routes


class TestListProjectFiles:
    def test_lists_folders_first_then_files(self, tmp_path):
        (tmp_path / "b.txt").write_text("x")
        (tmp_path / "A.txt").write_text("x")
        (tmp_path / "run").mkdir()

        assert ProjectFiles(str(tmp_path)).list() == [
            {"name": "run", "isFolder": True},
            {"name": "A.txt", "isFolder": False},
            {"name": "b.txt", "isFolder": False},
        ]

    def test_lists_a_subfolder(self, tmp_path):
        (tmp_path / "run").mkdir()
        (tmp_path / "run" / "out.nc").write_text("x")

        assert ProjectFiles(str(tmp_path)).list("run") == [{"name": "out.nc", "isFolder": False}]

    def test_nothing_above_the_project_folder(self, tmp_path):
        inside = tmp_path / "proj"
        inside.mkdir()
        (tmp_path / "secret.txt").write_text("x")

        assert ProjectFiles(str(inside)).list("..") == []
        assert ProjectFiles(str(inside)).list(str(tmp_path)) == []

    def test_nothing_through_a_link_out(self, tmp_path):
        inside = tmp_path / "proj"
        inside.mkdir()
        (tmp_path / "secret.txt").write_text("x")
        (inside / "out").symlink_to(tmp_path)

        assert ProjectFiles(str(inside)).list("out") == []

    def test_nothing_for_a_path_that_does_not_exist(self, tmp_path):
        assert ProjectFiles(str(tmp_path)).list("never-was") == []

    def test_nothing_for_a_file(self, tmp_path):
        (tmp_path / "a.txt").write_text("x")

        assert ProjectFiles(str(tmp_path)).list("a.txt") == []

    def test_nothing_with_no_project_folder(self, tmp_path):
        assert ProjectFiles().list() == []
        assert ProjectFiles(str(tmp_path / "gone")).list() == []

    def test_route_reads_the_current_folder(self, fake_hera, warmed, tmp_path):
        (tmp_path / "a.txt").write_text("x")
        db = fake_hera(["Alpha"], filesDirectories={"Alpha": str(tmp_path / "old")})
        db.filesDirectories["Alpha"] = str(tmp_path)

        assert routes.files("Alpha", "") == {"files": [{"name": "a.txt", "isFolder": False}]}

    def test_route_rejects_an_unknown_project(self, fake_hera, warmed):
        fake_hera(["Alpha"])

        with pytest.raises(HTTPException) as raised:
            routes.files("Nope", "")

        assert raised.value.status_code == 404
