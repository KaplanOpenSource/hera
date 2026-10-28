"""Serving files from disk and the built client. Its catch-all must be included last."""
import mimetypes
from pathlib import Path

from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse

# The built client (npm run build), served in production.
FRONTEND_DIST = Path(__file__).resolve().parents[2] / "client" / "bundle"

router = APIRouter()


@router.get("/file/{file_path:path}")
def serve_file(file_path: str):
    print('serving: ', file_path)
    full_path = Path("/") / file_path
    if not full_path.is_file():
        raise HTTPException(status_code=404, detail="File not found")
    media_type = mimetypes.guess_type(str(full_path))[0]
    print('mime: ', media_type)
    return FileResponse(str(full_path), media_type=media_type)


@router.get("/")
def serve_index_root():
    index_file = FRONTEND_DIST / "index.html"
    if index_file.exists():
        return FileResponse(str(index_file))
    return {"message": "Frontend build not found. Run 'npm run build'."}


# SPA fallback: serve index.html for other GET routes (after API routes)
@router.get("/{full_path:path}")
def spa_fallback(full_path: str):  # noqa: ARG001 (unused)
    index_file = FRONTEND_DIST / "index.html"
    if index_file.exists():
        return FileResponse(str(index_file))
    return {"message": "Not found"}
