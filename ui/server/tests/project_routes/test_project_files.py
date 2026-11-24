"""Which files belong to a project, and removing them."""
from project_api.project_files import ProjectFiles


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
