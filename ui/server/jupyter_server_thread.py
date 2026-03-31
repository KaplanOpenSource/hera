import asyncio
import threading

DEFAULT_JUPYTER_PORT = 8888


class JupyterServerThread:
    _port: int = DEFAULT_JUPYTER_PORT

    def __init__(self, port: int = DEFAULT_JUPYTER_PORT):
        JupyterServerThread._port = port
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()
        print(f"Jupyter server starting on port {port}")

    @staticmethod
    def port() -> int:
        return JupyterServerThread._port

    def _run(self):
        try:
            asyncio.set_event_loop(asyncio.new_event_loop())
            from jupyter_server.serverapp import ServerApp
            ServerApp.init_signal = lambda self: None  # signals only work in main thread
            app = ServerApp.instance()
            app.initialize([
                f'--port={self._port}',
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
