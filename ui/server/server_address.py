"""Where the UI server listens: pick a free port and tell the user the URLs."""
from __future__ import annotations

import socket


def find_available_port(host: str, start_port: int, attempts: int = 79) -> int:
    """Return the first free port at or after start_port, trying up to `attempts` ports."""
    for port in range(start_port, start_port + attempts):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            if sock.connect_ex((host, port)) != 0:
                return port
            print(f"Port {port} is in use, trying the next one...")
    raise SystemExit(
        f"No free port found in range {start_port}-{start_port + attempts - 1}."
    )


def ui_urls(host: str, port: int) -> list[str]:
    """Return the URLs a user can open in the browser to reach the UI."""
    if host != "0.0.0.0":
        return [f"http://{host}:{port}"]
    hostname = socket.gethostname()
    try:
        addresses = socket.gethostbyname_ex(hostname)[2]
    except socket.gaierror:
        # The machine name does not resolve; the name and localhost still work.
        addresses = []
    hosts = ["localhost", hostname] + [a for a in addresses if not a.startswith("127.")]
    return [f"http://{h}:{port}" for h in dict.fromkeys(hosts)]


def print_ui_urls(host: str, port: int) -> None:
    line = "=" * 60
    print(line)
    print("  Hera UI is starting. Open it in your browser at:")
    for url in ui_urls(host, port):
        print(f"    {url}")
    print(line, flush=True)
