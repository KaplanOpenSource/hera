"""datalayer/document/__init__.py: addOrUpdateDatabase/removeConnection,
which read and write $HOME/.pyhera/config.json. HOME is monkeypatched to
a tmp_path for every test here so the real user's config is never touched,
and the config file is pre-created since addOrUpdateDatabase and
removeConnection both raise FileNotFoundError if it doesn't already exist."""
import json

import pytest

from hera.datalayer.document import addOrUpdateDatabase, removeConnection


@pytest.fixture()
def config_path(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    cfg_dir = tmp_path / ".pyhera"
    cfg_dir.mkdir()
    path = cfg_dir / "config.json"
    path.write_text("{}")
    return path


@pytest.mark.unit
class TestAddOrUpdateDatabase:
    def test_it_adds_a_new_connection_entry(self, config_path):
        addOrUpdateDatabase("myconn", "u", "p", "1.2.3.4", "mydb")
        data = json.loads(config_path.read_text())
        assert data["myconn"] == {
            "username": "u", "password": "p", "dbIP": "1.2.3.4", "dbName": "mydb",
        }

    def test_it_overwrites_an_existing_entry_with_the_same_name(self, config_path):
        addOrUpdateDatabase("myconn", "u", "p", "1.2.3.4", "mydb")
        addOrUpdateDatabase("myconn", "u2", "p2", "5.6.7.8", "otherdb")
        data = json.loads(config_path.read_text())
        assert data["myconn"]["username"] == "u2"
        assert data["myconn"]["dbName"] == "otherdb"

    def test_it_leaves_other_connections_untouched(self, config_path):
        addOrUpdateDatabase("first", "u1", "p1", "1.1.1.1", "db1")
        addOrUpdateDatabase("second", "u2", "p2", "2.2.2.2", "db2")
        data = json.loads(config_path.read_text())
        assert set(data) == {"first", "second"}

    def test_a_missing_config_file_raises(self, tmp_path, monkeypatch):
        monkeypatch.setenv("HOME", str(tmp_path))
        with pytest.raises(FileNotFoundError):
            addOrUpdateDatabase("myconn", "u", "p", "1.2.3.4", "mydb")


@pytest.mark.unit
class TestRemoveConnection:
    def test_it_removes_the_named_connection(self, config_path):
        addOrUpdateDatabase("myconn", "u", "p", "1.2.3.4", "mydb")
        removeConnection("myconn")
        assert json.loads(config_path.read_text()) == {}

    def test_it_leaves_other_connections_in_place(self, config_path):
        addOrUpdateDatabase("first", "u1", "p1", "1.1.1.1", "db1")
        addOrUpdateDatabase("second", "u2", "p2", "2.2.2.2", "db2")
        removeConnection("first")
        assert set(json.loads(config_path.read_text())) == {"second"}

    def test_a_missing_config_file_raises(self, tmp_path, monkeypatch):
        monkeypatch.setenv("HOME", str(tmp_path))
        with pytest.raises(FileNotFoundError):
            removeConnection("myconn")
