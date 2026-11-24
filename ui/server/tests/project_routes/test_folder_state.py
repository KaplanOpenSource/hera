"""The state of one files folder."""
import os

import pytest

from project_api.folder_state import folder_state
from project_api.folder_state_enum import FolderStateEnum


class TestFolderState:
    def test_no_folder_set(self):
        assert folder_state("") == FolderStateEnum.NO_FOLDER

    def test_nothing_at_the_path(self, tmp_path):
        assert folder_state(str(tmp_path / "gone")) == FolderStateEnum.MISSING

    def test_a_file_at_the_path(self, tmp_path):
        (tmp_path / "file.txt").write_text("x")

        assert folder_state(str(tmp_path / "file.txt")) == FolderStateEnum.NOT_A_FOLDER

    def test_an_empty_folder(self, tmp_path):
        assert folder_state(str(tmp_path)) == FolderStateEnum.EMPTY

    def test_a_folder_with_only_the_case_config(self, tmp_path):
        (tmp_path / "caseConfiguration.json").write_text("{}")

        assert folder_state(str(tmp_path)) == FolderStateEnum.ONLY_CASE_CONFIG

    def test_a_folder_with_files(self, tmp_path):
        (tmp_path / "caseConfiguration.json").write_text("{}")
        (tmp_path / "data.parquet").write_text("x")

        assert folder_state(str(tmp_path)) == FolderStateEnum.HAS_FILES

    def test_a_folder_with_only_a_subfolder(self, tmp_path):
        (tmp_path / "run").mkdir()

        assert folder_state(str(tmp_path)) == FolderStateEnum.HAS_FILES

    @pytest.mark.skipif(os.geteuid() == 0, reason="root can list any folder")
    def test_a_folder_that_cannot_be_listed(self, tmp_path):
        locked = tmp_path / "locked"
        locked.mkdir()
        locked.chmod(0o000)
        try:
            assert folder_state(str(locked)) == FolderStateEnum.UNREADABLE
        finally:
            locked.chmod(0o700)

    def test_the_folders_answer_carries_the_state_string(self, fake_hera, tmp_path):
        from project_api.project_api import ProjectApi
        (tmp_path / "caseConfiguration.json").write_text("{}")
        fake_hera(["Alpha"], documents={
            "Alpha": [{"type": "Alpha__config__", "desc": {"filesDirectory": str(tmp_path)}}],
        })

        assert ProjectApi.folders()[0]["folderState"] == "onlyCaseConfig"
