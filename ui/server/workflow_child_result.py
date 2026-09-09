from __future__ import annotations

from typing import NamedTuple, Optional, Union


class WorkflowChildSuccess(NamedTuple):
    """The message the workflow child puts on the queue when the run succeeds."""

    dispatch_id: str
    exec_seconds: float
    chunks: Optional[list]  # per-task output segments; None on the subprocess path


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
    chunks: Optional[list]
