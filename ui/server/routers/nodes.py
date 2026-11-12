"""The Hermes node types the workflow editor offers, and their parameter checks."""
from fastapi import APIRouter

from api_models import ValidateNodeParamsPayload, ValidateNodeParamsResponse
from dependencies import warmup
from node_catalog import get_node_catalog
from node_params import validate_node_params

router = APIRouter()


@router.get("/node-catalog")
def node_catalog() -> list:
    """Available Hermes node types and the parameters each one accepts.

    See ``node_catalog.get_node_catalog``.
    """
    return get_node_catalog()


@router.post("/node-params/validate", response_model=ValidateNodeParamsResponse)
def node_params_validate(payload: ValidateNodeParamsPayload) -> ValidateNodeParamsResponse:
    """Test one node's parameter values with the node's own Hermes checks.

    See ``node_params.validate_node_params``. The checks reach hera and MongoDB,
    so before warmup there is nothing to report yet.
    """
    if not warmup.ready:
        return ValidateNodeParamsResponse()
    result = validate_node_params(payload.type, payload.params)
    if not result["ok"]:
        print(f'node params invalid: {payload.type}: {result["message"]}', flush=True)
    return ValidateNodeParamsResponse(**result)
