from __future__ import annotations

import importlib
import os
import sys
import time
import traceback
from multiprocessing.queues import Queue

import luigi

from hera.utils.logging import get_classMethod_logger

from .output_router import OutputRouter
from .task_pointer import task_pointer
from .workflow_child_result import WorkflowDone, WorkflowError, WorkflowMessage

# Imported for its side effect: defining the class registers the Luigi event handlers.
from . import luigi_task_events  # noqa: F401

# Number of Luigi workers when running in-process. 1 = sequential (today's behaviour).
LUIGI_WORKERS = 1

# The terminal task that triggers the whole DAG.
TARGET_TASK = "finalnode_xx_0"


def run_workflow_child_inprocess(
    project_name: str,
    workflow_name: str,
    result_queue: Queue[WorkflowMessage],
) -> None:
    """Build a saved workflow and run it with ``luigi.build`` inside a forked child.

    The router captures fd 1/2 and streams output to ``result_queue`` per task; the
    child ends with one ``WorkflowDone`` or a ``WorkflowError``. Prep is shared with
    the subprocess path via ``prepareWorkflowRunFromDoc``; only the execution differs
    (``luigi.build`` in this process). Running in a fresh forked process is what lets
    Luigi register its signal handlers and import the generated module fresh.
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
        # logger = get_classMethod_logger(workflow_toolkit, "run_workflow_child_inprocess")

        # The UI runs a single named workflow, so exactly one document is expected.
        docList = workflow_toolkit.getWorkflowListDocumentFromDB(workflow_name)
        if len(docList) != 1:
            raise RuntimeError(f"expected exactly one workflow for '{workflow_name}', got {len(docList)}")
        doc = docList[0]

        # Steps 1-4: rebuild, build, write, clean target files (shared with the subprocess path).
        pythonFileName, workflowName, targetFilesDir = workflow_toolkit.prepareWorkflowRunFromDoc(doc)
        print(f"[in-process luigi.build] target files for {workflowName} at {targetFilesDir}", flush=True)

        # Step 5: run in-process. The generated module needs FilesDirectory on
        # sys.path (PYTHONPATH only helps subprocesses).
        if workflow_toolkit.FilesDirectory not in sys.path:
            sys.path.insert(0, workflow_toolkit.FilesDirectory)
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
        router.stop()
        router = None
        result_queue.put(WorkflowDone(dispatch_id=None, exec_seconds=exec_seconds))
    except Exception:
        tb = traceback.format_exc()
        if router is not None:
            router.stop()
            router = None
        result_queue.put(WorkflowError(error=tb))
    finally:
        if router is not None:
            router.stop()
