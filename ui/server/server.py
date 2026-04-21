from __future__ import annotations

import argparse
import mimetypes
import traceback
from pathlib import Path

from typing import Any, Optional

from fastapi import FastAPI, HTTPException
from fastapi.encoders import jsonable_encoder
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from cors_handler import CorsHandler
from jupyter_server_thread import JupyterServerThread, DEFAULT_JUPYTER_PORT
from mock_data import MOCK_PROJECTS

cors_handler = CorsHandler()
parser = argparse.ArgumentParser(description="Hera UI API server")
cors_handler.add_argument(parser)
parser.add_argument('--debug', action='store_true', help='Enable debugpy remote debugging on port 5678')
parser.add_argument('-y', '--yes', action='store_true', help='Skip confirmation prompts')
parser.add_argument('--jupyter-port', type=int, default=8888, help='Port for Jupyter server (0 to disable)')
args = parser.parse_args()

app = FastAPI(title="Hera UI API")

origins = cors_handler.get_origins(args)

# Allow local Vite dev server
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# Mock data endpoints removed per simplified API


@app.get("/healthz")
def healthz() -> dict:
    return {"status": "ok"}


jupyter: JupyterServerThread | None = None


class JupyterStartPayload(BaseModel):
    root_dir: str


@app.post("/jupyter/ensure")
def jupyter_ensure(payload: JupyterStartPayload) -> dict:
    global jupyter
    jupyter_port = args.jupyter_port or DEFAULT_JUPYTER_PORT
    if jupyter and jupyter.is_running():
        if jupyter.root_dir == payload.root_dir:
            return {"port": jupyter.port, "root_dir": jupyter.root_dir}
        jupyter.stop()
    jupyter = JupyterServerThread(payload.root_dir, jupyter_port)
    jupyter.wait_until_ready()
    return {"port": jupyter.port, "root_dir": jupyter.root_dir}



@app.get("/cors")
def cors_info() -> dict:
    print ('cors', cors_handler.custom_origins)
    return {"origins": cors_handler.custom_origins}


class ExecPayload(BaseModel):
    code: str


class Problem(BaseModel):
    error: str
    traceback: str


class ExecResponse(BaseModel):
    data: Any = None
    problem: Optional[Problem] = None


# Code execution endpoint (simple: eval expression and return its value)
@app.post("/exec", response_model=ExecResponse)
def exec_code(payload: ExecPayload) -> ExecResponse:
    # DANGER: This is a security risk. It allows arbitrary code execution.
    # Only use this in a trusted environment.
    # The `_locals` dict will be updated with any variables created in the code.
    _locals = {} # "MOCK_PROJECTS": MOCK_PROJECTS}
    print("executing: " + payload.code)
    try:
        exec(payload.code, _locals)
    except Exception as e:
        error = f"{type(e).__name__}: {e}"
        tb = traceback.format_exc()
        print("exec error:", tb)
        return ExecResponse(problem=Problem(error=error, traceback=tb))
    result = _locals.get("result", None)
    print("got:", result)
    return ExecResponse(data=jsonable_encoder(result))


@app.get("/file/{file_path:path}")
def serve_file(file_path: str):
    print('serving: ', file_path)
    full_path = Path("/") / file_path
    if not full_path.is_file():
        raise HTTPException(status_code=404, detail="File not found")
    media_type = mimetypes.guess_type(str(full_path))[0]
    print('mime: ', media_type)
    return FileResponse(str(full_path), media_type=media_type)


# Serve built frontend (Vite) in production
FRONTEND_DIST = Path(__file__).resolve().parent.parent / "client" / "bundle"
# FRONTEND_DIST = Path("/client/dist").resolve()

if FRONTEND_DIST.exists():
    assets_dir = FRONTEND_DIST / "assets"
    if assets_dir.exists():
        app.mount(
            "/assets",
            StaticFiles(directory=str(assets_dir)),
            name="assets",
        )


@app.get("/")
def serve_index_root():
    index_file = FRONTEND_DIST / "index.html"
    if index_file.exists():
        return FileResponse(str(index_file))
    return {"message": "Frontend build not found. Run 'npm run build'."}


# SPA fallback: serve index.html for other GET routes (after API routes)
@app.get("/{full_path:path}")
def spa_fallback(full_path: str):  # noqa: ARG001 (unused)
    index_file = FRONTEND_DIST / "index.html"
    if index_file.exists():
        return FileResponse(str(index_file))
    return {"message": "Not found"}


if __name__ == "__main__":
    # Use a single process: no reload watcher
    import uvicorn

    # os.system("sh hera/scripts/run_mongo.sh")
    # os.system("sh hera/scripts/add_repo.sh")

    if args.debug:
        import debugpy

        debugpy.listen(("0.0.0.0", 5678))
        print("debugpy listening on 0.0.0.0:5678 - attach your debugger")

    uvicorn.run(
        app,
        host="0.0.0.0",
        port=8000,
        reload=False,
        workers=1,
    )
