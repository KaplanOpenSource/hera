"""utils/logging/helpers.py: add_FileHandler/add_formatter (pure dict builders)."""
import pytest

from hera.utils.logging import helpers


@pytest.mark.unit
class TestAddFileHandler:
    def test_it_builds_a_file_handler_config_tuple(self):
        name, cfg, section = helpers.add_FileHandler("h1", "out.log", mode="a")
        assert (name, section) == ("h1", "handlers")
        assert cfg == {"class": "logging.FileHandler", "filename": "out.log", "mode": "a", "formatter": "default"}

    def test_default_mode_is_overwrite(self):
        _, cfg, _ = helpers.add_FileHandler("h1", "out.log")
        assert cfg["mode"] == "w"


@pytest.mark.unit
class TestAddFormatter:
    def test_it_builds_a_formatter_config_tuple(self):
        name, cfg, section = helpers.add_formatter("f1", "%(message)s")
        assert (name, section) == ("f1", "formatters")
        assert cfg["format"] == "%(message)s"
        assert cfg["datefmt"] == "%Y-%m-%d %H:%M:%S"
