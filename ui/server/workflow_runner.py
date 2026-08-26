import time
import uuid
import threading
import multiprocessing

from pipe_tee import PipeTee
from workflow_child import run_workflow_in_child

# Job status values reported by poll(); "busy" is returned by start().
STATUS_RUNNING = "running"
STATUS_DONE = "done"
STATUS_ERROR = "error"
STATUS_BUSY = "busy"
STATUS_NOT_FOUND = "not_found"


def _now_readable() -> str:
    """Local time as ``YYYY-MM-DD HH:MM:SS.mmm`` (millis so quick polls differ)."""
    now = time.time()
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(now)) + f".{int((now % 1) * 1000):03d}"


class WorkflowRunner:
    """Builds and executes saved Hermes workflows in a separate process.

    ``start`` runs the workflow in the background and returns a token right away;
    the caller polls ``poll(token)`` until the run is done and then reads the whole
    output. One run at a time: ``start`` reports busy while a run is in progress.
    """

    def __init__(self):
        # Serialize runs: the local Luigi scheduler / DB access is not meant to run concurrently.
        self._lock = threading.Lock()
        # The single job slot; overwritten by the next start. None until the first run.
        self._job = None

    def start(self, project_name: str, workflow_name: str) -> dict:
        """Start a run in the background. Returns ``{"token"}`` or ``{"status": "busy"}``."""
        if self._job is not None and self._job["status"] == STATUS_RUNNING:
            return {"status": STATUS_BUSY}
        token = uuid.uuid4().hex
        self._job = {"token": token, "status": STATUS_RUNNING, "output": "", "error": ""}
        thread = threading.Thread(
            target=self._background,
            args=(token, project_name, workflow_name),
            daemon=True,
        )
        thread.start()
        return {"token": token}

    def poll(self, token: str) -> dict:
        """Return ``{"status", "output", "error"}`` for a token, or not_found if unknown."""
        job = self._job
        if job is None or job["token"] != token:
            print(f"[{_now_readable()}] poll {token}: not_found")
            return {"status": STATUS_NOT_FOUND, "output": "", "error": ""}
        print(f"[{_now_readable()}] poll {token}: {job['status']}")
        return {"status": job["status"], "output": job["output"], "error": job["error"]}

    def _background(self, token: str, project_name: str, workflow_name: str) -> None:
        # Runs in a background thread; record the outcome into the job slot for poll().
        try:
            result = self.run(project_name, workflow_name)
            status, output, error = STATUS_DONE, result["output"], ""
        except Exception as exc:
            # Surface the failure to the client via poll (this reports it, not hides it).
            status, output, error = STATUS_ERROR, "", str(exc)
            print("workflow run failed:", error)
        job = self._job
        if job is not None and job["token"] == token:
            job["output"] = output
            job["error"] = error
            job["status"] = status

    def run(self, project_name: str, workflow_name: str) -> dict:
        """Build and execute a saved workflow in a forked child process.

        Returns ``{"dispatch_id", "output"}``. Output is the child's captured
        console output with timing lines appended: how long the workflow itself
        ran and the total wall time including process spawn.
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
