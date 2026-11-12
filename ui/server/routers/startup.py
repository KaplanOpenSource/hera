"""The one call the UI makes to start: everything it needs, or "not ready yet"."""
from fastapi import APIRouter
from fastapi.encoders import jsonable_encoder

from dependencies import warmup
from init_data import collect_init

router = APIRouter()


@router.get("/init")
def init(project: str = "") -> dict:
    """Everything the UI needs to start, in one call. See ``init_data.collect_init``.

    Answers {"ready": false} while hera is still warming, so the client polls
    this one endpoint and the answer that says ready is also the data.
    """
    if not warmup.ready:
        return {"ready": False}
    return jsonable_encoder(collect_init(project))
