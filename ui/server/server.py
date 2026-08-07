from __future__ import annotations

import argparse, argcomplete
import mimetypes
import traceback
from pathlib import Path

from typing import Any, Optional

from fastapi import FastAPI, HTTPException
from fastapi.encoders import jsonable_encoder
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from cors_handler import CorsHandler
from jupyter_server_thread import JupyterServerThread, DEFAULT_JUPYTER_PORT
from node_catalog import get_node_catalog
from workflow_runner import WorkflowRunner

LOG_MAX_LEN = 350

cors_handler = CorsHandler()
parser = argparse.ArgumentParser(description="Hera UI API server")
cors_handler.add_argument(parser)
parser.add_argument('--debug', action='store_true', help='Enable debugpy remote debugging on port 5678')
parser.add_argument('-y', '--yes', action='store_true', help='Skip confirmation prompts')
parser.add_argument('--host', default='0.0.0.0', help='Address to bind to (default: 0.0.0.0, all interfaces; use 127.0.0.1 for local-only, e.g. behind a reverse proxy).')
parser.add_argument('--port', type=int, default=8000, help='Port for the API server')
parser.add_argument('--jupyter-port', type=int, default=8888, help='Port for Jupyter server (0 to disable)')
argcomplete.autocomplete(parser)
args = parser.parse_args()

app = FastAPI(title="Hera UI API")

origins = cors_handler.get_origins(args)


# Catch-all error reporter. Any unhandled exception (e.g. a broken import or a
# change in hera/hermes) is printed to the server console AND returned to the
# client as {error, traceback} with status 500 — the same content /exec already
# sends. Without this, FastAPI's default 500 is generated outside the CORS
# middleware, so it reaches the browser with no Access-Control-Allow-Origin
# header and shows up as a misleading "CORS error" that hides the real cause.
#
# Registered BEFORE the CORS middleware below on purpose: the last middleware
# added is the outermost, so adding CORS afterwards puts this handler *inside*
# CORS. Its JSON response then flows back out through CORS and gets the header.
@app.middleware("http")
async def report_unhandled_errors(request, call_next):
    try:
        return await call_next(request)
    except Exception as exc:
        tb = traceback.format_exc()
        print("server error:", tb)
        return JSONResponse(
            status_code=500,
            content={"error": f"{type(exc).__name__}: {exc}", "traceback": tb},
        )


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
    dark: bool = False


@app.post("/jupyter/ensure")
def jupyter_ensure(payload: JupyterStartPayload) -> dict:
    global jupyter
    jupyter_port = args.jupyter_port or DEFAULT_JUPYTER_PORT
    if jupyter and jupyter.is_running():
        if jupyter.root_dir == payload.root_dir:
            jupyter.set_theme(payload.dark)
            return {"port": jupyter.port, "root_dir": jupyter.root_dir}
        jupyter.stop()
    # Only expose the notebook server beyond localhost when CORS is enabled (i.e. the user
    # opted into remote access). With no --cors, keep it local-only regardless of --host.
    jupyter_ip = args.host if args.cors is not None else '127.0.0.1'
    jupyter = JupyterServerThread(payload.root_dir, jupyter_port, ip=jupyter_ip, dark=payload.dark)
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


def _truncate_for_log(value: Any) -> str:
    """Stringify a value for logging, cutting off long output so it doesn't flood the log."""
    text = repr(value)
    if len(text) > LOG_MAX_LEN:
        return f"{text[:LOG_MAX_LEN]}... [truncated, orig len {len(text)} chars]"
    return text


# Code execution endpoint (simple: eval expression and return its value)
@app.post("/exec", response_model=ExecResponse)
def exec_code(payload: ExecPayload) -> ExecResponse:
    # DANGER: This is a security risk. It allows arbitrary code execution.
    # Only use this in a trusted environment.
    # The `_locals` dict will be updated with any variables created in the code.
    _locals = {}
    print("executing: " + payload.code)
    try:
        exec(payload.code, _locals)
    except Exception as e:
        error = f"{type(e).__name__}: {e}"
        tb = traceback.format_exc()
        print("exec error:", tb)
        return ExecResponse(problem=Problem(error=error, traceback=tb))
    result = _locals.get("result", None)
    print("got:", _truncate_for_log(result))
    return ExecResponse(data=jsonable_encoder(result))


class RunWorkflowPayload(BaseModel):
    projectName: str
    workflowName: str


class RunWorkflowResponse(BaseModel):
    dispatch_id: str
    output: str


# Instantiated once and reused across requests (will hold run state/config later).
workflow_runner = WorkflowRunner()


@app.post("/run_workflow", response_model=RunWorkflowResponse)
def run_workflow(payload: RunWorkflowPayload) -> RunWorkflowResponse:
    """Build and execute a saved Hermes workflow from the DB (local Luigi scheduler).

    Returns the dispatch id and the run's captured console output. See WorkflowRunner.
    """
    result = workflow_runner.run(payload.projectName, payload.workflowName)
    return RunWorkflowResponse(**result)


@app.get("/node-catalog")
def node_catalog() -> list:
    """Available Hermes node types and the parameters each one accepts.

    See ``node_catalog.get_node_catalog``.
    """
    return get_node_catalog()


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


def find_available_port(host: str, start_port: int, attempts: int = 79) -> int:
    """Return the first free port at or after start_port, trying up to `attempts` ports."""
    import socket

    for port in range(start_port, start_port + attempts):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            if sock.connect_ex((host, port)) != 0:
                return port
            print(f"Port {port} is in use, trying the next one...")
    raise SystemExit(
        f"No free port found in range {start_port}-{start_port + attempts - 1}."
    )


if __name__ == "__main__":
    # Use a single process: no reload watcher
    import uvicorn

    # os.system("sh hera/scripts/run_mongo.sh")
    # os.system("sh hera/scripts/add_repo.sh")

    if args.debug:
        import debugpy

        debugpy.listen(("0.0.0.0", 5678))
        print("debugpy listening on 0.0.0.0:5678 - attach your debugger")

    # The bind address may be 0.0.0.0 (all interfaces); probe against localhost.
    probe_host = "127.0.0.1" if args.host == "0.0.0.0" else args.host
    port = find_available_port(probe_host, args.port)
    if port != args.port:
        print(f"Requested port {args.port} unavailable; starting on port {port} instead.")

    uvicorn.run(
        app,
        host=args.host,
        port=port,
        reload=False,
        workers=1,
    )
