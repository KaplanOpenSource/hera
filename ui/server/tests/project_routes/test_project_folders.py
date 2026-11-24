"""Each project with its files folder."""
import sys

from project_api.project_api import DEFAULT_PROJECT, ProjectApi
from routers import projects as routes


class TestProjectFolders:
    def test_each_project_with_its_folder_from_its_settings_record(self, fake_hera, warmed):
        fake_hera(["Beta", "Alpha"], documents={
            "Alpha": [{"type": "Alpha__config__", "desc": {"filesDirectory": "/data/a"}}],
            "Beta": [{"type": "Beta__config__", "desc": {"filesDirectory": "/data/b"}},
                     {"type": "Station", "desc": {"filesDirectory": "/not/this"}}],
        })

        assert routes.folders() == {"projects": [
            {"name": "Alpha", "filesDirectory": "/data/a", "folderExists": False, "documentCount": 0},
            {"name": "Beta", "filesDirectory": "/data/b", "folderExists": False, "documentCount": 1},
        ]}

    def test_a_project_without_a_settings_record_has_no_folder(self, fake_hera):
        fake_hera([DEFAULT_PROJECT])

        assert ProjectApi.folders() == [{"name": DEFAULT_PROJECT, "filesDirectory": "", "folderExists": False,
                                          "documentCount": 0}]

    def test_says_whether_the_folder_is_there(self, fake_hera, tmp_path):
        fake_hera(["Alpha", "Beta"], documents={
            "Alpha": [{"type": "Alpha__config__", "desc": {"filesDirectory": str(tmp_path)}}],
            "Beta": [{"type": "Beta__config__", "desc": {"filesDirectory": str(tmp_path / "gone")}}],
        })

        assert ProjectApi.folders() == [
            {"name": "Alpha", "filesDirectory": str(tmp_path), "folderExists": True, "documentCount": 0},
            {"name": "Beta", "filesDirectory": str(tmp_path / "gone"), "folderExists": False, "documentCount": 0},
        ]

    def test_counts_the_documents_but_not_the_settings_record(self, fake_hera):
        fake_hera(["Alpha"], documents={
            "Alpha": [{"type": "Alpha__config__", "desc": {"filesDirectory": "/data/a"}},
                      {"type": "Station", "desc": {}},
                      {"type": "Station", "desc": {}}],
        })

        assert ProjectApi.folders()[0]["documentCount"] == 2

    def test_does_not_open_or_change_any_project(self, fake_hera, warmed, monkeypatch):
        db = fake_hera(["Alpha"], documents={
            "Alpha": [{"type": "Alpha__config__", "desc": {"filesDirectory": "/data/a"}}],
        })
        def no_open(*args, **kwargs):
            raise AssertionError("folders opened a project")
        monkeypatch.setattr(sys.modules["hera.datalayer.project"], "Project", no_open)

        assert routes.folders() == {"projects": [{"name": "Alpha", "filesDirectory": "/data/a", "folderExists": False,
                                                       "documentCount": 0}]}
        assert db.documents == {"Alpha": [{"type": "Alpha__config__", "desc": {"filesDirectory": "/data/a"}}]}
        assert db.created == []
        assert db.filesDirectories == {}
