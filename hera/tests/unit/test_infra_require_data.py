"""The missing-data gate must skip locally and fail only when armed."""
import pytest


class _FakeConfig:
    def __init__(self, require_data):
        self._require_data = require_data

    def getoption(self, name):
        if name == "--require-data":
            return self._require_data
        raise KeyError(name)


class _FakeRequest:
    def __init__(self, require_data):
        self.config = _FakeConfig(require_data)


@pytest.mark.unit
def test_skips_by_default(monkeypatch):
    from hera.tests.conftest import _missing_test_data

    monkeypatch.delenv("HERA_REQUIRE_TEST_DATA", raising=False)
    with pytest.raises(pytest.skip.Exception):
        _missing_test_data(_FakeRequest(False), "no data")


@pytest.mark.unit
def test_fails_when_cli_flag_is_set(monkeypatch):
    from hera.tests.conftest import _missing_test_data

    monkeypatch.delenv("HERA_REQUIRE_TEST_DATA", raising=False)
    with pytest.raises(pytest.fail.Exception):
        _missing_test_data(_FakeRequest(True), "no data")


@pytest.mark.unit
@pytest.mark.parametrize("value", ["1", "true", "True", "yes"])
def test_fails_when_env_var_is_set(monkeypatch, value):
    from hera.tests.conftest import _missing_test_data

    monkeypatch.setenv("HERA_REQUIRE_TEST_DATA", value)
    with pytest.raises(pytest.fail.Exception):
        _missing_test_data(_FakeRequest(False), "no data")


@pytest.mark.unit
@pytest.mark.parametrize("value", ["", "0", "false", "no"])
def test_still_skips_for_falsey_env_values(monkeypatch, value):
    from hera.tests.conftest import _missing_test_data

    monkeypatch.setenv("HERA_REQUIRE_TEST_DATA", value)
    with pytest.raises(pytest.skip.Exception):
        _missing_test_data(_FakeRequest(False), "no data")
