"""The notebook server route: start one, reuse it, or restart it for a new folder."""
import pytest

import args_settings
from api_models import JupyterStartPayload
from routers import jupyter


class FakeJupyter:
    """Stands in for JupyterServerThread, so no real notebook server starts."""

    started = []

    def __init__(self, root_dir, port, ip="127.0.0.1", dark=False):
        self.root_dir = root_dir
        self.port = port
        self.ip = ip
        self.dark = dark
        self.running = True
        self.ready_waited = False
        FakeJupyter.started.append(self)

    def wait_until_ready(self):
        self.ready_waited = True
        return True

    def is_running(self):
        return self.running

    def set_theme(self, dark):
        self.dark = dark

    def stop(self):
        self.running = False


@pytest.fixture
def fake_jupyter(monkeypatch):
    FakeJupyter.started = []
    monkeypatch.setattr(jupyter, "JupyterServerThread", FakeJupyter)
    monkeypatch.setattr(jupyter, "jupyter", None)
    monkeypatch.setattr(args_settings, "_settings", args_settings.parse_args_settings(["--jupyter-port", "9999"]))
    return FakeJupyter


def _use_settings(monkeypatch, argv):
    settings = args_settings.parse_args_settings(argv)
    monkeypatch.setattr(args_settings, "_settings", settings)


def test_starts_a_server_and_waits_for_it(fake_jupyter):
    answer = jupyter.jupyter_ensure(JupyterStartPayload(root_dir="/data", dark=True))

    assert answer == {"port": 9999, "root_dir": "/data"}
    assert len(fake_jupyter.started) == 1
    started = fake_jupyter.started[0]
    assert started.ready_waited
    assert started.dark is True


def test_same_folder_reuses_the_server_and_sets_the_theme(fake_jupyter):
    jupyter.jupyter_ensure(JupyterStartPayload(root_dir="/data", dark=False))

    answer = jupyter.jupyter_ensure(JupyterStartPayload(root_dir="/data", dark=True))

    assert answer == {"port": 9999, "root_dir": "/data"}
    assert len(fake_jupyter.started) == 1
    assert fake_jupyter.started[0].dark is True


def test_new_folder_stops_the_old_server_and_starts_another(fake_jupyter):
    jupyter.jupyter_ensure(JupyterStartPayload(root_dir="/a"))

    answer = jupyter.jupyter_ensure(JupyterStartPayload(root_dir="/b"))

    assert answer == {"port": 9999, "root_dir": "/b"}
    first, second = fake_jupyter.started
    assert first.running is False
    assert second.root_dir == "/b"


def test_a_dead_server_is_started_again(fake_jupyter):
    jupyter.jupyter_ensure(JupyterStartPayload(root_dir="/data"))
    fake_jupyter.started[0].running = False

    jupyter.jupyter_ensure(JupyterStartPayload(root_dir="/data"))

    assert len(fake_jupyter.started) == 2


def test_stays_on_localhost_without_cors(fake_jupyter, monkeypatch):
    _use_settings(monkeypatch, ["--host", "0.0.0.0"])

    jupyter.jupyter_ensure(JupyterStartPayload(root_dir="/data"))

    assert fake_jupyter.started[0].ip == "127.0.0.1"


def test_binds_to_the_host_with_cors(fake_jupyter, monkeypatch):
    _use_settings(monkeypatch, ["--host", "0.0.0.0", "--cors", "10.0.0.1"])

    jupyter.jupyter_ensure(JupyterStartPayload(root_dir="/data"))

    assert fake_jupyter.started[0].ip == "0.0.0.0"
