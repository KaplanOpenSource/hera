import os
import sys
import time
import traceback


def run_workflow_child_subprocess(
    project_name: str,
    workflow_name: str,
    write_fd: int,
    result_queue,
) -> None:
    """Run a saved workflow inside a forked child, shelling Luigi out to a subprocess.

    stdout/stderr are redirected straight to ``write_fd`` (the pipe the parent reads),
    so all console output from Luigi/hera flows back to the server unchanged. The
    dispatch id and the execution time are sent to the parent through ``result_queue``.
    On failure the traceback is sent back instead so the parent can surface it. This
    path produces no per-task chunks (``chunks`` is None).
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
        dispatch_id = workflow_toolkit.executeWorkflowFromDB(workflow_name, scheduler="local")
        exec_seconds = time.perf_counter() - started

        sys.stdout.flush()
        sys.stderr.flush()
        result_queue.put({"dispatch_id": dispatch_id, "exec_seconds": exec_seconds, "chunks": None})
    except Exception:
        tb = traceback.format_exc()
        sys.stdout.flush()
        sys.stderr.flush()
        result_queue.put({"error": tb})
