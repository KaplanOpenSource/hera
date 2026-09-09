from __future__ import annotations

import os
import time
import traceback
from multiprocessing.queues import Queue

from output_router import OutputRouter
from task_pointer import task_pointer
from workflow_child_result import WorkflowDone, WorkflowError, WorkflowMessage

# Number of Luigi workers when running in-process. 1 = sequential (today's behaviour).
LUIGI_WORKERS = 1


def run_workflow_child_inprocess(
    project_name: str,
    workflow_name: str,
    result_queue: Queue[WorkflowMessage],
) -> None:
    """Run a saved workflow inside a forked child, executing Luigi in this process.

    The router captures fd 1/2 and streams the output to ``result_queue`` as
    ``WorkflowOutput`` messages, tagged per task. When the run finishes the child
    puts one ``WorkflowDone`` (dispatch id + timing); on failure it puts a
    ``WorkflowError`` with the traceback instead. The parent reads these off the queue.
    """
    router = OutputRouter(result_queue=result_queue, task_pointer=task_pointer)
    router.start()

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

        # Stop the router first so every output message is on the queue before done.
        router.stop()
        router = None
        result_queue.put(WorkflowDone(dispatch_id=dispatch_id, exec_seconds=exec_seconds))
    except Exception:
        tb = traceback.format_exc()
        if router is not None:
            router.stop()
            router = None
        result_queue.put(WorkflowError(error=tb))
    finally:
        if router is not None:
            router.stop()
