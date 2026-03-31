from __future__ import annotations

import asyncio
import threading
import urllib.request

DEFAULT_JUPYTER_PORT = 8888


class JupyterServerThread:
    _port: int = DEFAULT_JUPYTER_PORT
    _root_dir: str | None = None
    _thread: threading.Thread | None = None

    @staticmethod
    def port() -> int:
        return JupyterServerThread._port

    @staticmethod
    def root_dir() -> str | None:
        return JupyterServerThread._root_dir

    @staticmethod
    def is_running() -> bool:
        return JupyterServerThread._thread is not None and JupyterServerThread._thread.is_alive()

    @staticmethod
    def start(root_dir: str, port: int = DEFAULT_JUPYTER_PORT):
        if JupyterServerThread.is_running():
            if JupyterServerThread._root_dir == root_dir:
                return
            JupyterServerThread._stop()

        JupyterServerThread._port = port
        JupyterServerThread._root_dir = root_dir
        JupyterServerThread._thread = threading.Thread(
            target=JupyterServerThread._run,
            args=(port, root_dir),
            daemon=True,
        )
        JupyterServerThread._thread.start()
        print(f"Jupyter server starting on port {port}, root_dir={root_dir}")

    @staticmethod
    def _stop():
        if not JupyterServerThread.is_running():
            return
        try:
            port = JupyterServerThread._port
            req = urllib.request.Request(
                f'http://localhost:{port}/api/shutdown',
                method='POST',
                headers={'Content-Type': 'application/json'},
                data=b'{}',
            )
            urllib.request.urlopen(req, timeout=5)
        except Exception as e:
            print(f"Jupyter shutdown request: {e}")
        if JupyterServerThread._thread:
            JupyterServerThread._thread.join(timeout=5)
        JupyterServerThread._thread = None
        JupyterServerThread._root_dir = None
        print("Jupyter server stopped")

    @staticmethod
    def _run(port: int, root_dir: str):
        try:
            asyncio.set_event_loop(asyncio.new_event_loop())
            from jupyter_server.serverapp import ServerApp
            ServerApp.clear_instance()
            ServerApp.init_signal = lambda self: None  # signals only work in main thread
            app = ServerApp.instance()
            app.initialize([
                f'--port={port}',
                f'--notebook-dir={root_dir}',
                '--no-browser',
                '--allow-root',
                '--ServerApp.token=',
                '--ServerApp.password=',
                '--ServerApp.disable_check_xsrf=True',
                '--ServerApp.allow_origin=*',
                '--ServerApp.tornado_settings={"headers": {"Content-Security-Policy": "frame-ancestors *"}}',
            ])
            app.start()
        except ImportError:
            print("WARNING: jupyterlab is not installed — notebook server disabled. Install with: pip install jupyterlab")
        except Exception as e:
            print(f"ERROR: Jupyter server failed to start: {e}")
