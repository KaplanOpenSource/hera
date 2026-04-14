from __future__ import annotations

import asyncio
import json
import threading
import urllib.request
from pathlib import Path

DEFAULT_JUPYTER_PORT = 8888


class JupyterServerThread:
    def __init__(self, root_dir: str, port: int = DEFAULT_JUPYTER_PORT):
        self._port = port
        self._root_dir = root_dir
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()
        print(f"Jupyter server starting on port {port}, root_dir={root_dir}")

    def wait_until_ready(self, timeout: float = 30) -> bool:
        import time
        url = f"http://localhost:{self._port}/api/status"
        deadline = time.time() + timeout
        while time.time() < deadline:
            try:
                resp = urllib.request.urlopen(url, timeout=2)
                if resp.status == 200:
                    return True
            except Exception:
                pass
            time.sleep(0.5)
        return False

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

    @staticmethod
    def _enable_ai_magics():
        try:
            config_dir = Path.home() / '.ipython' / 'profile_default'
            config_dir.mkdir(parents=True, exist_ok=True)
            config_path = config_dir / 'ipython_config.py'
            line = "c.InteractiveShellApp.extensions = ['jupyter_ai_magics']"
            if config_path.exists() and line in config_path.read_text():
                return
            with open(config_path, 'a') as f:
                f.write(f"\n{line}\n")
        except OSError:
            print("WARNING: Could not write IPython config — %%ai magic requires manual %load_ext")

    @staticmethod
    def _disable_announcements():
        try:
            import sys
            settings_dir = Path(sys.prefix) / 'share' / 'jupyter' / 'lab' / 'settings'
            settings_dir.mkdir(parents=True, exist_ok=True)
            overrides_path = settings_dir / 'overrides.json'
            overrides = {}
            if overrides_path.exists():
                try:
                    overrides = json.loads(overrides_path.read_text())
                except Exception:
                    pass
            overrides['@jupyterlab/apputils-extension:notification'] = {
                'fetchNews': 'false',
                'checkForUpdates': False,
            }
            overrides['@jupyterlab/docmanager-extension:plugin'] = {
                'autosave': True,
                'autosaveInterval': 2,
            }
            overrides_path.write_text(json.dumps(overrides, indent=2))
        except OSError:
            print("WARNING: Could not write Jupyter overrides (permission denied) — news popup may appear")

    def _run(self):
        try:
            asyncio.set_event_loop(asyncio.new_event_loop())
            self._disable_announcements()
            self._enable_ai_magics()
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
