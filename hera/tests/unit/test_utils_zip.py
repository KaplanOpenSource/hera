"""zipUtils: building zip archives from files, directories and dicts.

These touch the filesystem, but only through tmp_path -- no test data set and
no fixed paths, so they stay in the hermetic layer.
"""
import json
import zipfile

import pytest

from hera.utils.zipUtils import add_directory_to_zip, list_json_files_in_zip, zip_items


@pytest.fixture()
def tree(tmp_path):
    """A small directory tree with a nested level, to check arcname handling."""
    root = tmp_path / "payload"
    (root / "inner").mkdir(parents=True)
    (root / "top.txt").write_text("top", encoding="utf-8")
    (root / "inner" / "deep.txt").write_text("deep", encoding="utf-8")
    return root


@pytest.mark.unit
class TestZipItems:
    def test_zip_extension_is_added_when_missing(self, tmp_path):
        target = tmp_path / "archive"
        zip_items(str(target), [{"a.txt": "content"}])
        assert (tmp_path / "archive.zip").exists()
        assert not target.exists()

    def test_existing_zip_extension_is_not_doubled(self, tmp_path):
        target = tmp_path / "archive.zip"
        zip_items(str(target), [{"a.txt": "content"}])
        assert target.exists()
        assert not (tmp_path / "archive.zip.zip").exists()

    def test_a_dict_becomes_files_named_by_its_keys(self, tmp_path):
        target = tmp_path / "d.zip"
        zip_items(str(target), [{"a.txt": "alpha", "b.txt": "beta"}])
        with zipfile.ZipFile(target) as archive:
            assert sorted(archive.namelist()) == ["a.txt", "b.txt"]
            assert archive.read("a.txt").decode() == "alpha"

    def test_a_file_is_stored_under_its_basename(self, tmp_path):
        source = tmp_path / "some" / "nested" / "file.txt"
        source.parent.mkdir(parents=True)
        source.write_text("payload", encoding="utf-8")

        target = tmp_path / "f.zip"
        zip_items(str(target), [str(source)])
        with zipfile.ZipFile(target) as archive:
            assert archive.namelist() == ["file.txt"]

    def test_a_directory_is_stored_with_its_own_name_as_root(self, tree, tmp_path):
        """arcname is relative to the directory's PARENT, so the name is kept."""
        target = tmp_path / "t.zip"
        zip_items(str(target), [str(tree)])
        with zipfile.ZipFile(target) as archive:
            names = sorted(archive.namelist())
        assert names == ["payload/inner/deep.txt", "payload/top.txt"]

    def test_mixed_items_in_one_archive(self, tree, tmp_path):
        loose = tmp_path / "loose.txt"
        loose.write_text("loose", encoding="utf-8")

        target = tmp_path / "m.zip"
        zip_items(str(target), [str(loose), str(tree), {"gen.txt": "generated"}])
        with zipfile.ZipFile(target) as archive:
            names = set(archive.namelist())
        assert {"loose.txt", "gen.txt", "payload/top.txt"} <= names

    def test_empty_item_list_makes_an_empty_archive(self, tmp_path):
        target = tmp_path / "e.zip"
        zip_items(str(target), [])
        with zipfile.ZipFile(target) as archive:
            assert archive.namelist() == []


@pytest.mark.unit
class TestZipItemsErrors:
    def test_non_string_dict_key_raises_typeerror(self, tmp_path):
        with pytest.raises(TypeError, match="Dict key must be a string"):
            zip_items(str(tmp_path / "x.zip"), [{1: "content"}])

    def test_non_string_dict_value_raises_typeerror(self, tmp_path):
        with pytest.raises(TypeError, match="Dict value must be a string"):
            zip_items(str(tmp_path / "x.zip"), [{"a.txt": 42}])

    def test_unsupported_item_type_raises_typeerror(self, tmp_path):
        with pytest.raises(TypeError, match="Unsupported item type"):
            zip_items(str(tmp_path / "x.zip"), [42])

    @pytest.mark.xfail(
        strict=True,
        reason="B5: the branch does `raise(f\"Type Error: ...\")` -- raising a "
               "string. Python reports 'exceptions must derive from BaseException' "
               "and the message naming the bad path never reaches the caller. "
               "See the consolidated findings issue.",
    )
    def test_a_path_that_is_neither_file_nor_directory_names_itself(self, tmp_path):
        """The caller needs to know WHICH item was wrong."""
        with pytest.raises(Exception, match="not a file or a directory"):
            zip_items(str(tmp_path / "x.zip"), [str(tmp_path / "absent")])


@pytest.mark.unit
class TestAddDirectoryToZip:
    def test_paths_are_relative_to_the_folder(self, tree, tmp_path):
        target = tmp_path / "a.zip"
        with zipfile.ZipFile(target, "w") as archive:
            add_directory_to_zip(archive, str(tree))
        with zipfile.ZipFile(target) as archive:
            assert sorted(archive.namelist()) == ["inner/deep.txt", "top.txt"]

    def test_zip_path_prefixes_every_entry(self, tree, tmp_path):
        target = tmp_path / "b.zip"
        with zipfile.ZipFile(target, "w") as archive:
            add_directory_to_zip(archive, str(tree), zip_path="under")
        with zipfile.ZipFile(target) as archive:
            assert sorted(archive.namelist()) == [
                "under/inner/deep.txt",
                "under/top.txt",
            ]


@pytest.mark.unit
class TestListJsonFilesInZip:
    def test_returns_name_and_parsed_content(self, tmp_path):
        target = tmp_path / "j.zip"
        zip_items(str(target), [{"conf.json": json.dumps({"a": 1})}])
        found = list_json_files_in_zip(str(target))
        assert found == [{"name": "conf.json", "content": {"a": 1}}]

    def test_non_json_members_are_skipped(self, tmp_path):
        target = tmp_path / "j.zip"
        zip_items(str(target), [{"conf.json": json.dumps({"a": 1}), "notes.txt": "hi"}])
        assert [entry["name"] for entry in list_json_files_in_zip(str(target))] == [
            "conf.json"
        ]

    def test_an_archive_without_json_yields_an_empty_list(self, tmp_path):
        target = tmp_path / "j.zip"
        zip_items(str(target), [{"notes.txt": "hi"}])
        assert list_json_files_in_zip(str(target)) == []

    def test_every_json_member_is_returned(self, tmp_path):
        target = tmp_path / "j.zip"
        zip_items(
            str(target),
            [{"a.json": json.dumps({"n": 1}), "b.json": json.dumps({"n": 2})}],
        )
        found = {entry["name"]: entry["content"]["n"] for entry in list_json_files_in_zip(str(target))}
        assert found == {"a.json": 1, "b.json": 2}

    def test_substring_match_picks_up_non_json_names(self, tmp_path):
        """Selection is `".json" in file_name`, not an extension check.

        Pinned because it means a backup file like conf.json.bak is parsed as
        JSON. Asserting it keeps the looseness visible rather than incidental.
        """
        target = tmp_path / "j.zip"
        zip_items(str(target), [{"conf.json.bak": json.dumps({"a": 1})}])
        assert [entry["name"] for entry in list_json_files_in_zip(str(target))] == [
            "conf.json.bak"
        ]

    @pytest.mark.xfail(
        strict=True,
        reason="B6: on UnicodeDecodeError the handler prints and falls through, "
               "leaving `content` unbound, so the next line raises "
               "UnboundLocalError instead of skipping the member. "
               "See the consolidated findings issue.",
    )
    def test_an_undecodable_json_member_is_skipped(self, tmp_path):
        """The code says 'Binary file, skipped reading content' -- so skip it."""
        target = tmp_path / "bad.zip"
        with zipfile.ZipFile(target, "w") as archive:
            archive.writestr("bad.json", b"\xff\xfe\x00binary")
            archive.writestr("good.json", json.dumps({"a": 1}))

        found = list_json_files_in_zip(str(target))
        assert [entry["name"] for entry in found] == ["good.json"]
