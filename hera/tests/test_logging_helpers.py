import warnings

from hera.utils.logging import helpers


class TestGetDefaultLoggingConfig:
    def test_no_deprecation_warning_when_creating_default_config(self, tmp_path, monkeypatch):
        """Loading the packaged default config must not use deprecated importlib.resources APIs."""
        log_dir = tmp_path / "log"
        log_dir.mkdir()
        monkeypatch.setattr(helpers, "HERA_DEFAULT_LOG_DIR", log_dir)

        with warnings.catch_warnings():
            warnings.simplefilter("error", DeprecationWarning)
            config = helpers.get_default_logging_config()

        assert isinstance(config, dict)
        assert (log_dir / "heraLogging.config").is_file()
