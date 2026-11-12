"""A real notebook server: it must report the port it really listens on."""
import socket
import urllib.request

import pytest

from jupyter_server_thread import JupyterServerThread

pytest.importorskip("jupyter_server")


def _free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def test_a_busy_port_moves_jupyter_and_reports_the_new_port(tmp_path):
    busy_port = _free_port()
    blocker = socket.socket()
    blocker.bind(("127.0.0.1", busy_port))
    blocker.listen()
    jupyter = JupyterServerThread(str(tmp_path), busy_port)

    assert jupyter.wait_until_ready()

    assert jupyter.port != busy_port
    answer = urllib.request.urlopen(f"http://localhost:{jupyter.port}/api/status", timeout=5)
    assert answer.headers.get("Content-Security-Policy").startswith("frame-ancestors *")
    jupyter.stop()
    blocker.close()
