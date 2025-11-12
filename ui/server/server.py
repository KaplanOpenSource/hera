
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

app = FastAPI(title="Hera UI API")

# Allow local Vite dev server
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:5173",
        "http://127.0.0.1:5173",
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# Mock data endpoints removed per simplified API


@app.get("/healthz")
def healthz() -> dict:
    return {"status": "ok"}


# Code execution endpoint (simple: eval expression and return its value)


class ExecPayload(BaseModel):
    code: str


@app.post("/exec")
def exec_code(payload: ExecPayload):
    # DANGER: This is a security risk. It allows arbitrary code execution.
    # Only use this in a trusted environment.
    # The `_locals` dict will be updated with any variables created in the code.
    _locals = {"MOCK_PROJECTS": MOCK_PROJECTS}
    print('executing: ' + payload.code)
    exec(payload.code, {}, _locals)
    result = _locals.get("result", None)
    print('got:', result)
    return jsonable_encoder(result)


# Serve built frontend (Vite) in production
FRONTEND_DIST = (
    Path(__file__).resolve().parent.parent / "client" / "dist"
)

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

    os.system("sh hera/scripts/run_mongo.sh")
    os.system("sh hera/scripts/add_repo.sh")

    if '--debug' in sys.argv:
        import debugpy
        debugpy.listen(("0.0.0.0", 5678))
        print("debugpy listening on 0.0.0.0:5678 - attach your debugger")

    uvicorn.run(
        "server:app",
        host="0.0.0.0",
        port=8000,
        reload=False,
        workers=1,
    )
