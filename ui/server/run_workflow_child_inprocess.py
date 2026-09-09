from __future__ import annotations

import os
import sys
import time
import traceback
from multiprocessing.queues import Queue

from output_router import OutputRouter
from run_chunk_state import state as chunk_state
from workflow_child_result import WorkflowChildError, WorkflowChildResult, WorkflowChildSuccess

# Number of Luigi workers when running in-process. 1 = sequential (today's behaviour).
LUIGI_WORKERS = 1


def run_workflow_child_inprocess(
    project_name: str,
    workflow_name: str,
    write_fd: int,
    result_queue: Queue[WorkflowChildResult],
) -> None:
    """Run a saved workflow inside a forked child, executing Luigi in this process.

    fd 1/2 are routed through OutputRouter so output is both forwarded to the server
    pipe (live, unchanged) and bucketed per task by the chunk pointer. The dispatch
    id, the execution time and the per-task chunks are sent to the parent through
    ``result_queue``. On failure the traceback is sent back instead so the parent can
    surface it.
    """
    router = OutputRouter(forward_fd=write_fd, state=chunk_state)
    router.start()
    os.close(write_fd)  # the router duped it; drop our extra copy

    try:
        from hera import toolkitHome

        workflow_toolkit = toolkitHome.getToolkit(
            toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS,
            projectName=project_name,
        )
        # The generated workflow module lives in the toolkit's files directory.
        os.environ["PYTHONPATH"] = workflow_toolkit.FilesDirectory + os.pathsep + os.environ.get("PYTHONPATH", "")

        started = time.perf_counter()
        from execute_workflow_inprocess import executeWorkflowFromDB_inprocess
        dispatch_id = executeWorkflowFromDB_inprocess(
            workflow_toolkit, workflow_name, workers=LUIGI_WORKERS
        )
        exec_seconds = time.perf_counter() - started

        # Stop the router first so all output is drained and bucketed before we read it.
        router.stop()
        router = None
        chunks = chunk_state.as_list()

        sys.stdout.flush()
        sys.stderr.flush()
        result_queue.put(WorkflowChildSuccess(dispatch_id=dispatch_id, exec_seconds=exec_seconds, chunks=chunks))
    except Exception:
        tb = traceback.format_exc()
        sys.stdout.flush()
        sys.stderr.flush()
        result_queue.put(WorkflowChildError(error=tb))
    finally:
        if router is not None:
            router.stop()
