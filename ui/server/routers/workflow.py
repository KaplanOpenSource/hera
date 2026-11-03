"""Starting a Hermes workflow run, polling it, stopping it, and reading its saved runs."""
from fastapi import APIRouter, HTTPException

from api_models import (
    RunWorkflowPayload,
    RunWorkflowResponse,
    WorkflowRunRecord,
    WorkflowRunsResponse,
    WorkflowRunSummary,
)
from dependencies import warmup, workflow_runner
from workflow_runner.workflow_run_history import WorkflowRunHistory

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


@router.post("/stop_workflow/{token}", response_model=RunWorkflowResponse)
def stop_workflow(token: str) -> RunWorkflowResponse:
    """Stop a running workflow by terminating its child process. See WorkflowRunner."""
    result = workflow_runner.stop(token)
    return RunWorkflowResponse(**result)


def _history(projectName: str, workflowName: str) -> WorkflowRunHistory:
    """The run history the request names. It reaches hera, so it waits for the warmup."""
    if not warmup.ready:
        raise HTTPException(status_code=503, detail="WARMING_UP")
    if not projectName:
        raise HTTPException(status_code=400, detail="no project name given")
    return WorkflowRunHistory(projectName, workflowName)


@router.get("/workflow_runs", response_model=WorkflowRunsResponse)
def workflow_runs(projectName: str, workflowName: str) -> WorkflowRunsResponse:
    """The saved runs of a workflow, newest first, without their output."""
    summaries = _history(projectName, workflowName).summaries()
    return WorkflowRunsResponse(runs=[WorkflowRunSummary(**summary) for summary in summaries])


@router.get("/workflow_run", response_model=WorkflowRunRecord)
def workflow_run(projectName: str, workflowName: str, runId: str) -> WorkflowRunRecord:
    """One saved run of a workflow, with its output."""
    record = _history(projectName, workflowName).read(runId)
    if record is None:
        raise HTTPException(status_code=404, detail=f"no such run: {runId}")
    return WorkflowRunRecord(**record)
