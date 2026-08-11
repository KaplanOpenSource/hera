import time
import threading
import multiprocessing

from pipe_tee import PipeTee
from workflow_child import run_workflow_in_child


class WorkflowRunner:
    """Builds and executes saved Hermes workflows in a separate process.

    Instantiated once and reused (it will hold run state/config later). The HTTP
    entrypoint stays synchronous: ``run`` blocks until the child process finishes.
    """

    def __init__(self):
        # Serialize runs: the local Luigi scheduler / DB access is not meant to run concurrently.
        self._lock = threading.Lock()

    def run(self, project_name: str, workflow_name: str) -> dict:
        """Build and execute a saved workflow in a forked child process.

        Returns ``{"dispatch_id", "output"}``. Output is the child's captured
        console output (echoed live to the server console) with timing lines
        appended: how long the workflow itself ran and the total wall time
        including process spawn.
        """
        with self._lock:
            # Fork so the child inherits the already-warmed hera import.
            ctx = multiprocessing.get_context("fork")
            result_queue = ctx.Queue()
            tee = PipeTee()

            total_started = time.perf_counter()
            process = ctx.Process(
                target=run_workflow_in_child,
                args=(project_name, workflow_name, tee.write_fd, result_queue),
            )
            process.start()
            # The write end belongs to the child now; drop ours so the reader sees EOF.
            tee.close_write()

            result = result_queue.get()
            process.join()
            total_seconds = time.perf_counter() - total_started

            output = tee.result()

            if "error" in result:
                raise RuntimeError(result["error"])

            timing = (
                f"\n[workflow ran in {result['exec_seconds']:.2f}s; "
                f"total {total_seconds:.2f}s including process spawn]\n"
            )
            return {"dispatch_id": result["dispatch_id"], "output": output + timing}
