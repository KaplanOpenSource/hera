import argparse
import os
import sys
from pathlib import Path

from fastapi import FastAPI
from fastapi.encoders import jsonable_encoder
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from mock_data import MOCK_PROJECTS

parser = argparse.ArgumentParser(description="Hera UI API server")
parser.add_argument(
    '--cors',
    nargs='?',
    const='*',
    default=None,
    metavar='ORIGINS',
    help=(
        'Enable CORS for external origins. '
        'Without a value, allows all origins (*). '
        'Pass a comma-separated list of IPs to allow specific ones '
        '(e.g. --cors 192.168.1.10,10.0.0.5). '
        'Each IP is prefixed with http:// and port 8000 automatically.'
    ),
)
parser.add_argument('--debug', action='store_true', help='Enable debugpy remote debugging on port 5678')
args = parser.parse_args()

app = FastAPI(title="Hera UI API")

LOCAL_ORIGINS = [f'http://{h}:{p}' for h in ['localhost', '127.0.0.1', '0.0.0.0'] for p in [5173, 8000]]

if args.cors is not None:
    if args.cors == '*':
        origins = ['*']
        print("WARNING: CORS is enabled for ALL origins (*). This is insecure in production.")
    else:
        cors_origins = [f'http://{ip}:8000' for ip in args.cors.split(',')]
        origins = LOCAL_ORIGINS + cors_origins
        print(f"WARNING: CORS is enabled for custom origins: {', '.join(cors_origins)}")
else:
    origins = LOCAL_ORIGINS

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
