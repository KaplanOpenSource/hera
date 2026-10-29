"""The Hera UI API server: arguments, the app and its routes, uvicorn.

The routes themselves live in ``routers/``; this file only assembles them.
"""
from __future__ import annotations

import socket
import traceback
from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles

from args_settings import ArgsSettings, parse_args_settings, set_args_settings
from dependencies import warmup
from routers import (cors_info, documents, files, geojson, health, jupyter, nodes,
                     projects, repositories, startup, workflow)


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Server is up and /init answers (ready=false) before warming begins.
    warmup.start()
    yield


# Catch-all error reporter. Any unhandled exception (e.g. a broken import or a
# change in hera/hermes) is printed to the server console AND returned to the
# client as {error, traceback} with status 500. Without this, FastAPI's default 500 is generated outside the CORS
# middleware, so it reaches the browser with no Access-Control-Allow-Origin
# header and shows up as a misleading "CORS error" that hides the real cause.
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


def create_app(settings: ArgsSettings) -> FastAPI:
    app = FastAPI(title="Hera UI API", lifespan=lifespan)

    # Registered BEFORE the CORS middleware below on purpose: the last middleware
    # added is the outermost, so adding CORS afterwards puts this handler *inside*
    # CORS. Its JSON response then flows back out through CORS and gets the header.
    app.middleware("http")(report_unhandled_errors)

    # Allow local Vite dev server
    app.add_middleware(
        CORSMiddleware,
        allow_origins=settings.origins(),
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    app.include_router(health.router)
    app.include_router(cors_info.router)
    app.include_router(startup.router)
    app.include_router(workflow.router)
    app.include_router(jupyter.router)
    app.include_router(nodes.router)
    app.include_router(geojson.router)
    app.include_router(projects.router)
    app.include_router(documents.router)
    app.include_router(repositories.router)

    assets_dir = files.FRONTEND_DIST / "assets"
    if assets_dir.exists():
        app.mount("/assets", StaticFiles(directory=str(assets_dir)), name="assets")

    # Last: its /{full_path:path} route would otherwise swallow every route after it.
    app.include_router(files.router)

    return app


def find_available_port(host: str, start_port: int, attempts: int = 79) -> int:
    """Return the first free port at or after start_port, trying up to `attempts` ports."""
    for port in range(start_port, start_port + attempts):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            if sock.connect_ex((host, port)) != 0:
                return port
            print(f"Port {port} is in use, trying the next one...")
    raise SystemExit(
        f"No free port found in range {start_port}-{start_port + attempts - 1}."
    )


def main() -> None:
    # Use a single process: no reload watcher
    import uvicorn

    settings = parse_args_settings()
    set_args_settings(settings)
    app = create_app(settings)

    if settings.debug:
        import debugpy

        debugpy.listen(("0.0.0.0", 5678))
        print("debugpy listening on 0.0.0.0:5678 - attach your debugger")

    # The bind address may be 0.0.0.0 (all interfaces); probe against localhost.
    probe_host = "127.0.0.1" if settings.host == "0.0.0.0" else settings.host
    port = find_available_port(probe_host, settings.port)
    if port != settings.port:
        print(f"Requested port {settings.port} unavailable; starting on port {port} instead.")

    uvicorn.run(
        app,
        host=settings.host,
        port=port,
        reload=False,
        workers=1,
    )


if __name__ == "__main__":
    main()
