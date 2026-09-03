import time
import uuid
import threading
import multiprocessing
from enum import Enum
from typing import Optional

from pipe_tee import PipeTee
from workflow_child import run_workflow_in_child


class RunStatus(str, Enum):
    """Run status reported by poll(); BUSY is returned by start().

    Subclasses ``str`` so members serialize to their plain string value over the
    wire (e.g. ``"running"``) and compare equal to it.
    """

    IDLE = "idle"  # no run started yet on this runner
    RUNNING = "running"
    DONE = "done"
    ERROR = "error"
    BUSY = "busy"
    NOT_FOUND = "not_found"


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
        # This runner holds one run. Its fields are overwritten by the next start;
        # more concurrent runs would mean more runners, not more fields here.
        # _token is None until the first run; _tee is the live output buffer, set
        # once the background thread starts the run so poll() can snapshot() it.
        self._token: Optional[str] = None
        self._status: RunStatus = RunStatus.IDLE
        self._output: str = ""
        self._error: str = ""
        self._chunks = None  # per-task output segments, filled in once the run is done
        self._tee: Optional[PipeTee] = None

    def start(self, project_name: str, workflow_name: str) -> dict:
        """Start a run in the background. Returns ``{"token"}`` or ``{"status": "busy"}``."""
        if self._status == RunStatus.RUNNING:
            return {"status": RunStatus.BUSY}
        self._token = uuid.uuid4().hex
        self._status = RunStatus.RUNNING
        self._output = ""
        self._error = ""
        self._chunks = None
        self._tee = None
        thread = threading.Thread(
            target=self._background,
            args=(self._token, project_name, workflow_name),
            daemon=True,
        )
        thread.start()
        return {"token": self._token}

    def poll(self, token: str) -> dict:
        """Return ``{"status", "output", "error"}`` for a token, or not_found if unknown."""
        # Read the token first: a mismatch means this isn't the run we hold.
        if self._token != token:
            print(f"[{_now_readable()}] poll {token}: not_found")
            return {"status": RunStatus.NOT_FOUND, "output": "", "error": ""}
        print(f"[{_now_readable()}] poll {token}: {self._status}")
        # While running, read the live buffer so the client sees output as it grows.
        # Once done, _output holds the final text (with the timing line).
        if self._status == RunStatus.RUNNING and self._tee is not None:
            output = self._tee.snapshot()
        else:
            output = self._output
        # chunks are only available once the run finished (the child returns them at
        # the end); while running they stay None and the client shows the flat log.
        return {"status": self._status, "output": output, "error": self._error, "chunks": self._chunks}

    def _background(self, token: str, project_name: str, workflow_name: str) -> None:
        # Runs in a background thread; record the outcome for poll(). The token guard
        # keeps a finished run from clobbering a newer one that took the runner over.
        # Create the tee here and store it so poll() can read partial output mid-run.
        tee = PipeTee()
        if self._token == token:
            self._tee = tee
        chunks = None
        try:
            result = self.run(project_name, workflow_name, tee)
            status, output, error = RunStatus.DONE, result["output"], ""
            chunks = result.get("chunks")
        except Exception as exc:
            # Surface the failure to the client via poll (this reports it, not hides it).
            status, output, error = RunStatus.ERROR, "", str(exc)
            print("workflow run failed:", error)
        if self._token == token:
            self._output = output
            self._error = error
            self._chunks = chunks
            self._status = status

    def run(self, project_name: str, workflow_name: str, tee: Optional[PipeTee] = None) -> dict:
        """Build and execute a saved workflow in a forked child process.

        Returns ``{"dispatch_id", "output"}``. Output is the child's captured
        console output with timing lines appended: how long the workflow itself
        ran and the total wall time including process spawn.

        ``tee`` captures the child's output; when omitted a fresh one is made.
        ``start`` passes one it also stored on the job so ``poll`` can read the
        output while the run is still going.
        """
        if tee is None:
            tee = PipeTee()
        with self._lock:
            # Fork so the child inherits the already-warmed hera import.
            ctx = multiprocessing.get_context("fork")
            result_queue = ctx.Queue()

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
            # chunks: per-task output buckets from the in-process router (None on the
            # subprocess path). Passed through so callers can show output per task.
            return {
                "dispatch_id": result["dispatch_id"],
                "output": output + timing,
                "chunks": result.get("chunks"),
            }
