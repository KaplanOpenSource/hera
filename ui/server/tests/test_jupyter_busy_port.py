"""Issue #1123: someone else's Jupyter holds the port, so the UI must be sent to ours."""
import socket
import subprocess
import sys
import time
import urllib.request

import pytest

import args_settings
from api_models import JupyterStartPayload
from routers import jupyter

pytest.importorskip("jupyterlab")


def _free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def _frame_rule(port):
    answer = urllib.request.urlopen(f"http://localhost:{port}/api/status", timeout=5)
    return answer.headers.get("Content-Security-Policy")


@pytest.fixture
def other_jupyter(tmp_path):
    """A plain Jupyter with default settings, like one a user left running."""
    port = _free_port()
    process = subprocess.Popen(
        [sys.executable, "-m", "jupyterlab", f"--port={port}", "--ServerApp.ip=127.0.0.1", "--no-browser",
         f"--notebook-dir={tmp_path}", "--ServerApp.token=", "--ServerApp.password="],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    deadline = time.time() + 30
    while time.time() < deadline:
        try:
            _frame_rule(port)
            break
        except OSError:
            time.sleep(0.5)
    yield port
    process.terminate()
    process.wait()


def test_the_ui_is_sent_to_our_jupyter_not_the_other_one(other_jupyter, tmp_path, monkeypatch):
    monkeypatch.setattr(jupyter, "jupyter", None)
    monkeypatch.setattr(args_settings, "_settings", args_settings.parse_args_settings(["--jupyter-port", str(other_jupyter)]))
    assert _frame_rule(other_jupyter).startswith("frame-ancestors 'self'")

    answer = jupyter.jupyter_ensure(JupyterStartPayload(root_dir=str(tmp_path)))

    assert answer["port"] != other_jupyter
    assert _frame_rule(answer["port"]).startswith("frame-ancestors *")
    jupyter.jupyter.stop()
