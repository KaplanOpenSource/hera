"""The unit layer must be unable to reach the network or the real home."""
import os
import socket
import tempfile

import pytest


@pytest.mark.unit
def test_network_is_blocked():
    with pytest.raises(RuntimeError, match="network access"):
        socket.socket(socket.AF_INET, socket.SOCK_STREAM)


@pytest.mark.unit
def test_matplotlib_backend_is_agg():
    import matplotlib

    assert matplotlib.get_backend().lower() == "agg"


@pytest.mark.unit
def test_home_is_the_isolated_one():
    assert os.environ["HOME"].startswith(
        tempfile.gettempdir()
    ), "HOME was not isolated; a test could pollute the real ~/.hera"
