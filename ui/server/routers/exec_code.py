"""Arbitrary python execution, the UI's general purpose backdoor."""
from typing import Any
import traceback

from fastapi import APIRouter
from fastapi.encoders import jsonable_encoder

from api_models import ExecPayload, ExecResponse, Problem
from dependencies import warmup

LOG_MAX_LEN = 350

router = APIRouter()


def _truncate_for_log(value: Any) -> str:
    """Stringify a value for logging, cutting off long output so it doesn't flood the log."""
    text = repr(value)
    if len(text) > LOG_MAX_LEN:
        return f"{text[:LOG_MAX_LEN]}... [truncated, orig len {len(text)} chars]"
    return text


@router.post("/exec", response_model=ExecResponse)
def exec_code(payload: ExecPayload) -> ExecResponse:
    # DANGER: This is a security risk. It allows arbitrary code execution.
    # Only use this in a trusted environment.
    # The `_locals` dict will be updated with any variables created in the code.
    if not warmup.ready:
        # Marker so the client waits/retries instead of showing an error.
        return ExecResponse(problem=Problem(error="WARMING_UP", traceback=""))
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
