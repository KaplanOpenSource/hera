from __future__ import annotations

import argparse
import mimetypes
from pathlib import Path

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
    return {"port": jupyter.port, "root_dir": jupyter.root_dir}


class NotebookPayload(BaseModel):
    root_dir: str


@app.get("/notebooks")
def list_notebooks(root_dir: str) -> dict:
    notebooks_dir = Path(root_dir) / "notebooks"
    if not notebooks_dir.exists():
        return {"notebooks": []}
    files = sorted(
        f.name for f in notebooks_dir.iterdir()
        if f.suffix == ".ipynb" and f.is_file()
    )
    return {"notebooks": files}


@app.post("/notebooks/create")
def create_notebook(payload: NotebookPayload) -> dict:
    import json as _json
    notebooks_dir = Path(payload.root_dir) / "notebooks"
    notebooks_dir.mkdir(parents=True, exist_ok=True)
    existing = [
        int(f.stem.split("_")[1])
        for f in notebooks_dir.iterdir()
        if f.suffix == ".ipynb" and f.stem.startswith("notebook_") and f.stem.split("_")[1].isdigit()
    ]
    next_id = max(existing, default=0) + 1
    name = f"notebook_{next_id}.ipynb"
    empty_notebook = {
        "nbformat": 4,
        "nbformat_minor": 5,
        "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}},
        "cells": [],
    }
    (notebooks_dir / name).write_text(_json.dumps(empty_notebook, indent=2))
    return {"name": name}


@app.get("/cors")
def cors_info() -> dict:
    print ('cors', cors_handler.custom_origins)
    return {"origins": cors_handler.custom_origins}


class ExecPayload(BaseModel):
    code: str


# Code execution endpoint (simple: eval expression and return its value)
@app.post("/exec")
def exec_code(payload: ExecPayload):
    # DANGER: This is a security risk. It allows arbitrary code execution.
    # Only use this in a trusted environment.
    # The `_locals` dict will be updated with any variables created in the code.
    _locals = {} # "MOCK_PROJECTS": MOCK_PROJECTS}
    print("executing: " + payload.code)
    exec(payload.code, {}, _locals)
    result = _locals.get("result", None)
    print("got:", result)
    return jsonable_encoder(result)


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
