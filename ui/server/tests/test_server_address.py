"""The startup banner lists the URLs a user can open to reach the UI."""
import socket

import server_address


def test_a_fixed_host_gives_one_url():
    assert server_address.ui_urls("10.0.0.5", 8001) == ["http://10.0.0.5:8001"]


def test_all_interfaces_gives_localhost_name_and_address(monkeypatch):
    monkeypatch.setattr(socket, "gethostname", lambda: "int-math551")
    monkeypatch.setattr(
        socket, "gethostbyname_ex",
        lambda name: (name, [], ["127.0.1.1", "192.168.14.107"]),
    )

    assert server_address.ui_urls("0.0.0.0", 8001) == [
        "http://localhost:8001",
        "http://int-math551:8001",
        "http://192.168.14.107:8001",
    ]


def test_an_unresolvable_name_still_gives_localhost_and_name(monkeypatch):
    def fail(name):
        raise socket.gaierror("no such name")

    monkeypatch.setattr(socket, "gethostname", lambda: "box")
    monkeypatch.setattr(socket, "gethostbyname_ex", fail)

    assert server_address.ui_urls("0.0.0.0", 8000) == ["http://localhost:8000", "http://box:8000"]


def test_the_banner_prints_every_url(monkeypatch, capsys):
    monkeypatch.setattr(server_address, "ui_urls", lambda host, port: ["http://a:1", "http://b:1"])

    server_address.print_ui_urls("0.0.0.0", 1)

    out = capsys.readouterr().out
    assert "http://a:1" in out
    assert "http://b:1" in out
