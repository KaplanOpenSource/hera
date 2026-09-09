from __future__ import annotations

import importlib
import os
import sys
import time
import traceback
from multiprocessing.queues import Queue

import luigi

from .output_router import OutputRouter
from .task_pointer import task_pointer
from .workflow_child_result import WorkflowDone, WorkflowError, WorkflowMessage

# Imported for its side effect: defining the class registers the Luigi event handlers.
from . import luigi_task_events  # noqa: F401

# Number of Luigi workers when running in-process. 1 = sequential (today's behaviour).
LUIGI_WORKERS = 1

# The terminal task that triggers the whole DAG.
TARGET_TASK = "finalnode_xx_0"


class WorkflowChildInProcess:
    """Build a saved workflow and run it with ``luigi.build`` inside a forked child.

    The router captures fd 1/2 and streams output to ``result_queue`` per task; the
    run ends with one ``WorkflowDone`` or a ``WorkflowError``. Prep is shared with
    the subprocess path via ``prepareWorkflowRunFromDoc``; only the execution differs
    (``luigi.build`` in this process). Running in a fresh forked process is what lets
    Luigi register its signal handlers and import the generated module fresh.
    """

    def __init__(
        self,
        project_name: str,
        workflow_name: str,
        result_queue: Queue[WorkflowMessage],
    ) -> None:
        self.project_name = project_name
        self.workflow_name = workflow_name
        self.result_queue = result_queue
        self.router = OutputRouter(result_queue=result_queue, task_pointer=task_pointer)

    @staticmethod
    def add_python_path(files_directory: str) -> None:
        # Make the generated module importable: sys.path for the in-process import,
        # PYTHONPATH for any subprocess Luigi may spawn.
        if files_directory not in sys.path:
            sys.path.insert(0, files_directory)
        os.environ["PYTHONPATH"] = files_directory + os.pathsep + os.environ.get("PYTHONPATH", "")

    def run(self) -> None:
        self.router.start()

        try:
            from hera import toolkitHome

            workflow_toolkit = toolkitHome.getToolkit(
                toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS,
                projectName=self.project_name,
            )
            # The generated workflow module lives in the toolkit's files directory.
            self.add_python_path(workflow_toolkit.FilesDirectory)

            started = time.perf_counter()

            # The UI runs a single named workflow, so exactly one document is expected.
            docList = workflow_toolkit.getWorkflowListDocumentFromDB(self.workflow_name)
            if len(docList) != 1:
                raise RuntimeError(f"expected exactly one workflow for '{self.workflow_name}', got {len(docList)}")
            doc = docList[0]

            # Steps 1-4: rebuild, build, write, clean target files (shared with the subprocess path).
            pythonFileName, workflowName, targetFilesDir = workflow_toolkit.prepareWorkflowRunFromDoc(doc)
            print(f"[in-process luigi.build] target files for {workflowName} at {targetFilesDir}", flush=True)

            # Step 5: run in-process.
            if workflowName in sys.modules:
                # Re-run in the same process: reload the just-regenerated .py.
                mod = importlib.reload(sys.modules[workflowName])
            else:
                mod = importlib.import_module(workflowName)

            print(f"[in-process luigi.build] Running {workflowName} in-process with luigi.build (pid={os.getpid()} workers={LUIGI_WORKERS})")

            taskClass = getattr(mod, TARGET_TASK)
            # No dispatch id: local scheduler, one serialized run, so nodes use the flat layout.
            finalTask = taskClass()
            ok = luigi.build([finalTask], workers=LUIGI_WORKERS, local_scheduler=True)
            if not ok:
                raise RuntimeError(f"workflow {workflowName} failed (luigi.build returned False)")

            # Step 6: Clean up the generated Python module (the workflow JSON stays).
            print(f"Cleaning the executer python for {workflowName}")
            os.remove(pythonFileName)

            exec_seconds = time.perf_counter() - started

            # Stop the router first so every output message is on the queue before done.
            self.router.stop()
            self.result_queue.put(WorkflowDone(dispatch_id=None, exec_seconds=exec_seconds))
        except Exception:
            tb = traceback.format_exc()
            self.router.stop()
            self.result_queue.put(WorkflowError(error=tb))
        finally:
            self.router.stop()

    @staticmethod
    def start_child(
        project_name: str,
        workflow_name: str,
        result_queue: Queue[WorkflowMessage],
    ) -> None:
        """Process entry point: run one saved workflow in this (forked) process."""
        WorkflowChildInProcess(project_name, workflow_name, result_queue).run()
