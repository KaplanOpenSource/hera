"""Starting a Hermes workflow run and polling it."""
from fastapi import APIRouter

from api_models import RunWorkflowPayload, RunWorkflowResponse
from dependencies import workflow_runner

router = APIRouter()


# Renamed from the old synchronous POST /run_workflow so a stale caller of that
# path gets a loud 404 instead of silently misreading the new token response.
@router.post("/start_workflow", response_model=RunWorkflowResponse)
def start_workflow(payload: RunWorkflowPayload) -> RunWorkflowResponse:
    """Start a saved Hermes workflow run in the background (local Luigi scheduler).

    Returns a token to poll with, or status "busy" if a run is already in progress.
    Poll GET /workflow_status/{token} for status and chunks. See WorkflowRunner.
    """
    result = workflow_runner.start(payload.projectName, payload.doc)
    return RunWorkflowResponse(**result)


@router.get("/workflow_status/{token}", response_model=RunWorkflowResponse)
def workflow_status(token: str) -> RunWorkflowResponse:
    """Report a run's status and its captured output as per-task chunks. See WorkflowRunner."""
    result = workflow_runner.poll(token)
    return RunWorkflowResponse(**result)
