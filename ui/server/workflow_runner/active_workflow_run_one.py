import os
import queue
import signal
import time
import uuid
import threading
import multiprocessing
from multiprocessing.process import BaseProcess
from multiprocessing.queues import Queue
from typing import Any, Dict, Optional, cast

from .run_status import RunStatus
from .task_pointer import BETWEEN
from .run_workflow_child_inprocess import WorkflowChildInProcess
from .workflow_child_result import WorkflowError, WorkflowMessage, WorkflowOutput, WorkflowRunResult
from .workflow_log_builder import WorkflowLogBuilder
from .workflow_run_history import WorkflowRunHistory

# How long the read loop waits on the queue before it checks that the child is alive.
CHILD_CHECK_SECONDS = 0.2


class ActiveWorkflowRunOne:
    """One workflow run, in a forked child process. Saved to the run history when it ends.

    ``start`` runs it in a background thread; ``poll`` reads its state and output;
    ``stop`` ends the child and every command it started.
    """

    def __init__(self, project_name: str, doc: Dict[str, Any]):
        self.project_name = project_name
        self.doc = doc
        self.doc_id: str = doc["_id"]["$oid"]
        self.workflow_name: str = doc["desc"]["workflowName"]
        self.token = uuid.uuid4().hex
        started = time.localtime()
        # The saved run's file name; one run at a time, so the start second is unique.
        self.run_id = time.strftime("%Y-%m-%d_%H-%M-%S", started)
        self.start_time = time.strftime("%Y-%m-%d %H:%M:%S", started)
        self.end_time = ""
        self.status = RunStatus.RUNNING
        self.error = ""
        self.chunks = None  # per-task output segments, filled in once the run is done
        # The live log builder, so poll() can read partial output while the run goes.
        self.log = WorkflowLogBuilder()
        # Set once the child process has started.
        self.process: Optional[BaseProcess] = None
        self.pid: Optional[int] = None
        self.stop_requested = False

    def start(self) -> None:
        """Run in a background thread; poll() reports the outcome."""
        thread = threading.Thread(target=self._background, daemon=True)
        thread.start()

    def stop(self) -> None:
        """Terminate the child and the commands it started. The status becomes STOPPED."""
        self.stop_requested = True
        if self.pid is not None:
            os.killpg(self.pid, signal.SIGTERM)

    def files_directory(self) -> str:
        return WorkflowRunHistory.files_directory(self.project_name)

    def target_files_folder(self) -> str:
        """The Luigi target files folder. It is deleted at the start of every run."""
        return os.path.join(self.files_directory(), f"{self.workflow_name}_targetFiles")

    def run_history_folder(self) -> str:
        return WorkflowRunHistory(self.project_name, self.workflow_name).folder()

    def save_history(self, status: RunStatus) -> None:
        """Write this finished run, which ended with ``status``, to the run history folder."""
        # A project with no files folder cannot run, so there is nothing to keep.
        if not self.files_directory():
            return
        WorkflowRunHistory(self.project_name, self.workflow_name).save(self.run_id, {
            "token": self.token,
            "startTime": self.start_time,
            "endTime": self.end_time,
            "status": status.value,
            "error": self.error,
            "chunks": self.chunks,
        })

    def poll(self) -> dict:
        """Return ``{"status", "error", "chunks"}``."""
        # While running, read the live builder so the client sees chunks as they grow.
        # Once done, chunks holds the final value (with the timing line).
        if self.status == RunStatus.RUNNING:
            chunks = self.log.chunks()
        else:
            chunks = self.chunks
        return {"status": self.status, "error": self.error, "chunks": chunks}

    def _background(self) -> None:
        try:
            result = self.run()
            chunks = result.chunks
            if self.stop_requested:
                status, error = RunStatus.STOPPED, ""
            else:
                status, error = RunStatus.DONE, ""
        except Exception as exc:
            # Surface the failure to the client via poll (this reports it, not hides it).
            status, error = RunStatus.ERROR, str(exc)
            # Keep the output captured up to the failure, not just the error.
            chunks = self.log.chunks()
            print("workflow run failed:", error)
        self.error = error
        self.chunks = chunks
        self.end_time = time.strftime("%Y-%m-%d %H:%M:%S")
        # Save before the status changes, so a client that sees the end finds the file.
        try:
            self.save_history(status)
        finally:
            self.status = status

    def run(self) -> WorkflowRunResult:
        """Build and execute the workflow in a forked child process.

        Returns the per-task chunks (with the timing line as a final chunk) and the
        dispatch id. Chunks hold the child's captured console output; the timing line
        says how long the workflow ran and the total wall time including process spawn.
        """
        # Fork so the child inherits the already-warmed hera import.
        ctx = multiprocessing.get_context("fork")
        result_queue: Queue[WorkflowMessage] = ctx.Queue()

        total_started = time.perf_counter()
        process = ctx.Process(
            target=WorkflowChildInProcess.start_child,
            args=(self.project_name, self.doc, result_queue),
        )
        process.start()
        # The pid is only None before start(), so it is an int here.
        pid = cast(int, process.pid)
        # Also set the group here, in case the child has not set it yet.
        os.setpgid(pid, pid)
        self.process = process
        self.pid = pid
        # A stop that came before the process existed could not terminate it.
        if self.stop_requested:
            os.killpg(pid, signal.SIGTERM)

        # Drain the queue live: output messages feed the log builder, until the
        # child sends its one terminal message (done or error), or dies.
        while True:
            try:
                message: WorkflowMessage = result_queue.get(timeout=CHILD_CHECK_SECONDS)
            except queue.Empty:
                # Wake now and then, so a child that died (stopped) does not hang the loop.
                if process.is_alive():
                    continue
                process.join()
                if self.stop_requested:
                    self.log.add(WorkflowOutput(task=BETWEEN, text="\n[workflow stopped]\n"))
                    return WorkflowRunResult(dispatch_id=None, chunks=self.log.chunks())
                raise RuntimeError(f"workflow process ended with no result (exit code {process.exitcode})")

            if isinstance(message, WorkflowOutput):
                self.log.add(message)
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
            self.log.add(WorkflowOutput(task=BETWEEN, text=timing))
            return WorkflowRunResult(
                dispatch_id=message.dispatch_id,
                chunks=self.log.chunks(),
            )
