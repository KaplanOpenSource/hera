"""The CORS origins the server was started with."""
from fastapi import APIRouter

from args_settings import get_args_settings

router = APIRouter()


@router.get("/cors")
def cors_info() -> dict:
    origins = get_args_settings().custom_origins
    print('cors', origins)
    return {"origins": origins}
