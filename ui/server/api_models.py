"""Request/response models for the Hera UI API endpoints."""

from typing import Any, Dict, List, Optional

from pydantic import BaseModel


class JupyterStartPayload(BaseModel):
    root_dir: str
    dark: bool = False


class ExecPayload(BaseModel):
    code: str


class Problem(BaseModel):
    error: str
    traceback: str


class ExecResponse(BaseModel):
    data: Any = None
    problem: Optional[Problem] = None


class RunWorkflowPayload(BaseModel):
    projectName: str
    # The whole workflow document (desc.workflow, desc.workflowName, resource, ...).
    # The client always sends it, so the run builds straight from it, no DB lookup.
    doc: Dict[str, Any]


class WorkflowChunk(BaseModel):
    # One output segment: a task's name (or "__preamble__" / "__between__") and the
    # console output captured while that segment was current.
    name: str
    text: str


class RunWorkflowResponse(BaseModel):
    # start returns token (or status "busy"); poll returns status + chunks/error.
    token: Optional[str] = None
    status: Optional[str] = None
    error: str = ""
    # Per-task output segments, in run order. Grows live while the run is going;
    # None only before any output (or when the token is unknown).
    chunks: Optional[List[WorkflowChunk]] = None


class ValidateNodeParamsPayload(BaseModel):
    # The node's dotted Hermes type, e.g. "RiskAssessment.calculateThresholds".
    type: str
    # The node's input_parameters, exactly as the editor holds them.
    params: Dict[str, Any] = {}


class ValidateNodeParamsResponse(BaseModel):
    # ok with an empty message also means "nothing to say" (unknown type, or a
    # node with no checks of its own).
    ok: bool = True
    message: str = ""
