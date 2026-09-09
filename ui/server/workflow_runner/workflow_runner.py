import time
import uuid
import threading
import multiprocessing
from multiprocessing.queues import Queue
from typing import Optional

from .run_status import RunStatus
from .task_pointer import BETWEEN
from time_utils import now_readable
from .run_workflow_child_inprocess import WorkflowChildInProcess
from .workflow_child_result import WorkflowError, WorkflowMessage, WorkflowOutput, WorkflowRunResult
from .workflow_log_builder import WorkflowLogBuilder


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
        # _token is None until the first run; _log is the live log builder, set once
        # the background thread starts the run so poll() can read partial output.
        self._token: Optional[str] = None
        self._status: RunStatus = RunStatus.IDLE
        self._output: str = ""
        self._error: str = ""
        self._chunks = None  # per-task output segments, filled in once the run is done
        self._log: Optional[WorkflowLogBuilder] = None

    def start(self, project_name: str, workflow_name: str) -> dict:
        """Start a run in the background. Returns ``{"token"}`` or ``{"status": "busy"}``."""
        if self._status == RunStatus.RUNNING:
            return {"status": RunStatus.BUSY}
        self._token = uuid.uuid4().hex
        self._status = RunStatus.RUNNING
        self._output = ""
        self._error = ""
        self._chunks = None
        self._log = None
        thread = threading.Thread(
            target=self._background,
            args=(self._token, project_name, workflow_name),
            daemon=True,
        )
        thread.start()
        return {"token": self._token}

    def poll(self, token: str) -> dict:
        """Return ``{"status", "output", "error", "chunks"}`` for a token, or not_found."""
        # Read the token first: a mismatch means this isn't the run we hold.
        if self._token != token:
            print(f"[{now_readable()}] poll {token}: not_found")
            return {"status": RunStatus.NOT_FOUND, "output": "", "error": ""}
        print(f"[{now_readable()}] poll {token}: {self._status}")
        # While running, read the live builder so the client sees output as it grows.
        # Once done, _output / _chunks hold the final values (with the timing line).
        if self._status == RunStatus.RUNNING and self._log is not None:
            output = self._log.output()
            chunks = self._log.chunks()
        else:
            output = self._output
            chunks = self._chunks
        return {"status": self._status, "output": output, "error": self._error, "chunks": chunks}

    def _background(self, token: str, project_name: str, workflow_name: str) -> None:
        # Runs in a background thread; record the outcome for poll(). The token guard
        # keeps a finished run from clobbering a newer one that took the runner over.
        # Create the log builder here and store it so poll() can read partial output.
        log = WorkflowLogBuilder()
        if self._token == token:
            self._log = log
        chunks = None
        try:
            result = self.run(project_name, workflow_name, log)
            status, output, error = RunStatus.DONE, result.output, ""
            chunks = result.chunks
        except Exception as exc:
            # Surface the failure to the client via poll (this reports it, not hides it).
            status, output, error = RunStatus.ERROR, "", str(exc)
            print("workflow run failed:", error)
        if self._token == token:
            self._output = output
            self._error = error
            self._chunks = chunks
            self._status = status

    def run(self, project_name: str, workflow_name: str, log: Optional[WorkflowLogBuilder] = None) -> WorkflowRunResult:
        """Build and execute a saved workflow in a forked child process.

        Returns the flat log (with timing lines appended), the per-task chunks, and
        the dispatch id. Output is the child's captured console output; timing lines
        say how long the workflow ran and the total wall time including process spawn.

        ``log`` accumulates the child's output; when omitted a fresh one is made.
        ``start`` passes one it also stored on the job so ``poll`` can read the
        output while the run is still going.
        """
        if log is None:
            log = WorkflowLogBuilder()
        with self._lock:
            # Fork so the child inherits the already-warmed hera import.
            ctx = multiprocessing.get_context("fork")
            result_queue: Queue[WorkflowMessage] = ctx.Queue()

            total_started = time.perf_counter()
            process = ctx.Process(
                target=WorkflowChildInProcess.start_child,
                args=(project_name, workflow_name, result_queue),
            )
            process.start()

            # Drain the queue live: output messages feed the log builder, until the
            # child sends its one terminal message (done or error).
            while True:
                message: WorkflowMessage = result_queue.get()
                if isinstance(message, WorkflowOutput):
                    log.add(message)
                    continue

                # Terminal message: the run is over either way.
                process.join()
                total_seconds = time.perf_counter() - total_started

                if isinstance(message, WorkflowError):
                    raise RuntimeError(message.error)

                # message is a WorkflowDone here.
                timing = (
                    f"\n[workflow ran in {message.exec_seconds:.2f}s; "
                    f"total {total_seconds:.2f}s including process spawn]\n"
                )
                # Append the timing line as a final output piece so it shows in the log.
                log.add(WorkflowOutput(task=BETWEEN, text=timing))
                return WorkflowRunResult(
                    dispatch_id=message.dispatch_id,
                    output=log.output(),
                    chunks=log.chunks(),
                )
