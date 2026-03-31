from __future__ import annotations

import asyncio
import threading
import urllib.request

DEFAULT_JUPYTER_PORT = 8888


class JupyterServerThread:
    def __init__(self, root_dir: str, port: int = DEFAULT_JUPYTER_PORT):
        self._port = port
        self._root_dir = root_dir
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()
        print(f"Jupyter server starting on port {port}, root_dir={root_dir}")

    @property
    def port(self) -> int:
        return self._port

    @property
    def root_dir(self) -> str:
        return self._root_dir

    def is_running(self) -> bool:
        return self._thread.is_alive()

    def stop(self):
        if not self.is_running():
            return
        try:
            req = urllib.request.Request(
                f'http://localhost:{self._port}/api/shutdown',
                method='POST',
                headers={'Content-Type': 'application/json'},
                data=b'{}',
            )
            urllib.request.urlopen(req, timeout=5)
        except Exception as e:
            print(f"Jupyter shutdown request: {e}")
        self._thread.join(timeout=5)
        print("Jupyter server stopped")

    def _run(self):
        try:
            asyncio.set_event_loop(asyncio.new_event_loop())
            from jupyter_server.serverapp import ServerApp
            ServerApp.clear_instance()
            ServerApp.init_signal = lambda self: None  # signals only work in main thread
            app = ServerApp.instance()
            app.initialize([
                f'--port={self._port}',
                f'--notebook-dir={self._root_dir}',
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
