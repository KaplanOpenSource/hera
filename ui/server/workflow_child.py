import os
import sys
import time
import traceback

# Run Luigi in this process via luigi.build instead of shelling out to
# `python -m luigi`. Lets us set workers and hook Luigi events. Flip to try it.
INPROCESS_LUIGI = True
# Number of Luigi workers when running in-process. 1 = sequential (today's behaviour).
LUIGI_WORKERS = 1


def run_workflow_in_child(project_name: str, workflow_name: str, write_fd: int, result_queue) -> None:
    """Run a saved workflow inside a forked child process.

    stdout/stderr are redirected to ``write_fd`` (a pipe the parent reads), so all
    console output from Luigi/hera flows back to the server. The dispatch id and the
    execution time are sent to the parent through ``result_queue``. On failure the
    traceback is sent back instead so the parent can surface it (mirrors the old
    in-process behaviour where the error reached the HTTP 500 handler).
    """
    os.dup2(write_fd, 1)
    os.dup2(write_fd, 2)
    os.close(write_fd)

    try:
        from hera import toolkitHome

        workflow_toolkit = toolkitHome.getToolkit(
            toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS,
            projectName=project_name,
        )
        # The generated workflow module lives in the toolkit's files directory.
        os.environ["PYTHONPATH"] = workflow_toolkit.FilesDirectory + os.pathsep + os.environ.get("PYTHONPATH", "")

        started = time.perf_counter()
        if INPROCESS_LUIGI:
            from execute_workflow_inprocess import executeWorkflowFromDB_inprocess
            dispatch_id = executeWorkflowFromDB_inprocess(
                workflow_toolkit, workflow_name, workers=LUIGI_WORKERS
            )
        else:
            dispatch_id = workflow_toolkit.executeWorkflowFromDB(workflow_name, scheduler="local")
        exec_seconds = time.perf_counter() - started

        sys.stdout.flush()
        sys.stderr.flush()
        result_queue.put({"dispatch_id": dispatch_id, "exec_seconds": exec_seconds})
    except Exception:
        tb = traceback.format_exc()
        sys.stdout.flush()
        sys.stderr.flush()
        result_queue.put({"error": tb})
