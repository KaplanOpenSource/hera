"""The notebook server the UI opens notebooks in."""
from fastapi import APIRouter

from api_models import JupyterStartPayload
from jupyter_server_thread import JupyterServerThread, DEFAULT_JUPYTER_PORT
from args_settings import get_args_settings

router = APIRouter()

jupyter: JupyterServerThread | None = None


@router.post("/jupyter/ensure")
def jupyter_ensure(payload: JupyterStartPayload) -> dict:
    global jupyter
    settings = get_args_settings()
    jupyter_port = settings.jupyter_port or DEFAULT_JUPYTER_PORT
    if jupyter and jupyter.is_running():
        if jupyter.root_dir == payload.root_dir:
            jupyter.set_theme(payload.dark)
            return {"port": jupyter.port, "root_dir": jupyter.root_dir}
        jupyter.stop()
    # Only expose the notebook server beyond localhost when CORS is enabled (i.e. the user
    # opted into remote access). With no --cors, keep it local-only regardless of --host.
    jupyter_ip = settings.host if settings.remote_allowed else '127.0.0.1'
    jupyter = JupyterServerThread(payload.root_dir, jupyter_port, ip=jupyter_ip, dark=payload.dark)
    jupyter.wait_until_ready()
    return {"port": jupyter.port, "root_dir": jupyter.root_dir}
