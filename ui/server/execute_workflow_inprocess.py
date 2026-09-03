"""In-process twin of ``hermesWorkflowToolkit.executeWorkflowFromDB``.

Same build-and-run steps as the toolkit's method, but step 5 calls
``luigi.build(...)`` in the current process instead of shelling out to
``python3 -m luigi ...``. This keeps the toolkit code untouched while letting the
server run Luigi in-process: it can set ``workers`` for parallelism and register
Luigi event handlers (they fire in this same process).

Run this inside the forked workflow child (see ``workflow_child.py``): that makes
it the main thread of a fresh process, so Luigi's signal handler registration
works and each run imports the generated module fresh.
"""

import importlib
import os
import shutil
import sys

import luigi
from hermes import workflow

from hera.utils.logging import get_classMethod_logger

from run_chunk_state import BETWEEN, state as chunk_state


# Every Luigi event line starts with this prefix so the UI log parser can spot
# them and hide them by default, without confusing them for a task's own output.
# Kept in sync with EVENT_PREFIX in the client's classifyLog.ts.
EVENT_PREFIX = "[luigi-event]"


def _event(message):
    """Print one Luigi-event line with the shared prefix (flushed for live output)."""
    print(f"{EVENT_PREFIX} {message}", flush=True)


# Register once at import: these fire for every luigi.Task run in this process.
# START/SUCCESS/FAILURE also move the chunk pointer so the output router buckets
# each task's output under the task's name (see run_chunk_state / output_router).
@luigi.Task.event_handler(luigi.Event.START)
def _on_task_start(task):
    # Point at this task first, so the START line and the run's output land in its bucket.
    chunk_state.current = task.task_family
    _event(f"START {task.task_family}")


@luigi.Task.event_handler(luigi.Event.SUCCESS)
def _on_task_success(task):
    # Print while still on the task's bucket, then stop pointing at it.
    _event(f"SUCCESS {task.task_family}")
    chunk_state.current = BETWEEN


@luigi.Task.event_handler(luigi.Event.FAILURE)
def _on_task_failure(task, exception):
    _event(f"FAILURE {task.task_family}: {exception}")
    chunk_state.current = BETWEEN


@luigi.Task.event_handler(luigi.Event.PROCESSING_TIME)
def _on_task_time(task, seconds):
    _event(f"TIME {task.task_family}: {seconds:.2f}s")


@luigi.Task.event_handler(luigi.Event.BROKEN_TASK)
def _on_task_broken(task, exception):
    _event(f"BROKEN {task.task_family}: {exception}")
    chunk_state.current = BETWEEN


@luigi.Task.event_handler(luigi.Event.PROGRESS)
def _on_task_progress(task, progress):
    _event(f"PROGRESS {task.task_family}: {progress}")


@luigi.Task.event_handler(luigi.Event.DEPENDENCY_DISCOVERED)
def _on_dependency_discovered(task, dependency):
    _event(f"DEP DISCOVERED {task.task_family} -> {dependency.task_family}")


@luigi.Task.event_handler(luigi.Event.DEPENDENCY_MISSING)
def _on_dependency_missing(task):
    _event(f"DEP MISSING {task.task_family}")


@luigi.Task.event_handler(luigi.Event.DEPENDENCY_PRESENT)
def _on_dependency_present(task):
    _event(f"DEP PRESENT {task.task_family}")


def executeWorkflowFromDB_inprocess(
    workflow_toolkit,
    nameOrWorkflowFileOrJSONOrResource,
    dispatch_id=None,
    workers=1,
    targetTask="finalnode_xx_0",
):
    """Build a saved workflow and run it in-process with ``luigi.build``.

    Mirrors ``hermesWorkflowToolkit.executeWorkflowFromDB`` step for step; only the
    execution (step 5) differs. Returns the ``dispatch_id`` used for the run.

    Parameters
    ----------
    workflow_toolkit : hermesWorkflowToolkit
        A toolkit instance (from ``toolkitHome.getToolkit``). Used for the DB
        lookup, the files directory, and the workflow rebuild helpers.
    nameOrWorkflowFileOrJSONOrResource : str or dict
        Same accepted forms as the toolkit method (name, resource, workflow, dict).
    dispatch_id : str, optional
        Per-run id propagated to every Luigi node. A fresh uuid-like id is NOT
        generated here on purpose: pass one in if you need dispatch isolation.
        When None the nodes fall back to their empty-default (legacy flat layout).
    workers : int
        Number of Luigi workers. 1 keeps today's sequential behaviour; >1 runs
        independent tasks in parallel (each in its own forked process).
    targetTask : str
        The terminal task that triggers the whole DAG (default ``finalnode_xx_0``).
    """
    logger = get_classMethod_logger(workflow_toolkit, "executeWorkflowFromDB_inprocess")
    docList = workflow_toolkit.getWorkflowListDocumentFromDB(nameOrWorkflowFileOrJSONOrResource)

    logger.info(f"In-process execution with dispatch_id='{dispatch_id}' workers={workers}")

    for doc in docList:
        workflowJSON = doc.desc["workflow"]
        workflowName = doc.desc["workflowName"]
        logger.info(f"Processing {workflowName}")

        # Step 1: Reconstruct the hermes workflow object from the stored JSON.
        hermesWF = workflow_toolkit.getHermesWorkflowFromJSON(
            workflowJSON, name=workflowName, resource=doc["resource"]
        )

        # Step 2: Build the workflow into a Luigi task module (Python source).
        logger.info(f"Building the workflow {workflowName}")
        build = hermesWF.build(buildername=workflow.BUILDER_LUIGI)

        # Step 3: Write the workflow JSON and the generated Python module to disk.
        logger.info(f"Writing the workflow and the executer python {workflowName}")
        wfFileName = hermesWF.Resource_path
        hermesWF.write(wfFileName)

        pythonFileName = os.path.join(workflow_toolkit.FilesDirectory, f"{workflowName}.py")
        with open(pythonFileName, "w") as outFile:
            outFile.write(build)

        # Step 4: Clean previous execution artifacts (Luigi target files) so every
        # task re-runs from scratch.
        logger.debug("Removing the targetfiles")
        executionfileDir = os.path.join(
            workflow_toolkit.FilesDirectory, f"{workflowName}_targetFiles"
        )
        shutil.rmtree(executionfileDir, ignore_errors=True)

        # Step 5: Execute in-process via luigi.build instead of a subprocess.
        # The generated module lives in FilesDirectory; put it on sys.path so the
        # in-process import resolves (the PYTHONPATH env var only helps subprocesses).
        if workflow_toolkit.FilesDirectory not in sys.path:
            sys.path.insert(0, workflow_toolkit.FilesDirectory)
        logger.debug(f"Running {workflowName} in-process with luigi.build (workers={workers})")
        if workflowName in sys.modules:
            # Same-process re-run: the .py was just regenerated, so reload it rather
            # than reuse the cached (possibly stale) module object.
            mod = importlib.reload(sys.modules[workflowName])
        else:
            mod = importlib.import_module(workflowName)

        taskClass = getattr(mod, targetTask)
        finalTask = taskClass(dispatch_id=dispatch_id) if dispatch_id is not None else taskClass()
        print(f"[in-process luigi.build] pid={os.getpid()} module={workflowName} workers={workers}", flush=True)
        ok = luigi.build([finalTask], workers=workers, local_scheduler=True)
        if not ok:
            raise RuntimeError(f"workflow {workflowName} failed (luigi.build returned False)")

        # Step 6: Clean up the generated Python module (the workflow JSON stays).
        logger.info(f"Cleaning the executer python for {workflowName}")
        os.remove(pythonFileName)

    return dispatch_id
