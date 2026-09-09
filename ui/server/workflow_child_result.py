from __future__ import annotations

from typing import NamedTuple, Union


class WorkflowChildSuccess(NamedTuple):
    """The message the workflow child puts on the queue when the run succeeds."""

    dispatch_id: str
    exec_seconds: float
    chunks: list  # per-task output segments from the in-process router


class WorkflowChildError(NamedTuple):
    """The message the workflow child puts on the queue when the run fails."""

    error: str  # a traceback string


# One message is put on the result queue per run: either a success or an error.
WorkflowChildResult = Union[WorkflowChildSuccess, WorkflowChildError]


class WorkflowRunResult(NamedTuple):
    """What ``WorkflowRunner.run`` returns to its caller: the child's output with
    timing lines appended, plus the dispatch id and the per-task chunks."""

    dispatch_id: str
    output: str
    chunks: list
