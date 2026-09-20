from __future__ import annotations

from typing import NamedTuple, Optional, Union


class WorkflowOutput(NamedTuple):
    """A piece of captured output, tagged with the task running when it was captured.

    The child streams many of these on the result queue as the run produces output.
    ``task`` is the task name, or ``PREAMBLE`` / ``BETWEEN`` (see ``task_pointer``).
    """

    task: str
    text: str


class WorkflowDone(NamedTuple):
    """Sent once on the result queue when the run finishes successfully.

    ``dispatch_id`` is None when the run used no dispatch id (the legacy flat layout);
    see ``run_workflow_child_inprocess``.
    """

    dispatch_id: Optional[str]
    exec_seconds: float


class WorkflowError(NamedTuple):
    """Sent once on the result queue if the run raises. Carries the traceback."""

    error: str  # a traceback string


# The child streams these on the result queue: many WorkflowOutput, then exactly
# one WorkflowDone (success) or WorkflowError (failure).
WorkflowMessage = Union[WorkflowOutput, WorkflowDone, WorkflowError]


class WorkflowRunResult(NamedTuple):
    """What ``WorkflowRunner.run`` returns to its caller: the per-task chunks (with
    the timing line as a final chunk) and the dispatch id (None if unused)."""

    dispatch_id: Optional[str]
    chunks: list
