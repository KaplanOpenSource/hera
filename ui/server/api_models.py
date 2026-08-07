"""Request/response models for the Hera UI API endpoints."""

from typing import Any, Optional

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


class RunWorkflowResponse(BaseModel):
    dispatch_id: str
    output: str
