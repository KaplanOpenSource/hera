"""Request/response models for the Hera UI API endpoints."""

from typing import Any, List, Optional

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
    workflowName: str


class WorkflowChunk(BaseModel):
    # One output segment: a task's name (or "__preamble__" / "__between__") and the
    # console output captured while that segment was current.
    name: str
    text: str


class RunWorkflowResponse(BaseModel):
    # start returns token (or status "busy"); poll returns status + output/error.
    token: Optional[str] = None
    status: Optional[str] = None
    output: str = ""
    error: str = ""
    # Per-task output segments, in run order. Filled in only once the run is done
    # (in-process path only); None while running or on the subprocess path.
    chunks: Optional[List[WorkflowChunk]] = None
