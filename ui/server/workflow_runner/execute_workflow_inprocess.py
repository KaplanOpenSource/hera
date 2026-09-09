"""In-process twin of ``hermesWorkflowToolkit.executeWorkflowFromDB``.

Shares the build steps with the toolkit's method via
``prepareWorkflowRunFromDoc`` (rebuild, build, write, clean target files), but
runs the result by calling ``luigi.build(...)`` in the current process instead
of shelling out to ``python3 -m luigi ...``. Running Luigi in-process lets the
server set ``workers`` for parallelism and register Luigi event handlers (they
fire in this same process).

Run this inside the forked workflow child (see ``run_workflow_child_inprocess.py``): that makes
it the main thread of a fresh process, so Luigi's signal handler registration
works and each run imports the generated module fresh.
"""

from __future__ import annotations

import importlib
import os
import sys
from typing import TYPE_CHECKING

import luigi

from hera.utils.logging import get_classMethod_logger

# Imported for its side effect: defining the class registers the Luigi event handlers.
from . import luigi_task_events  # noqa: F401

if TYPE_CHECKING:
    from hera.simulations.hermesWorkflowToolkit import hermesWorkflowToolkit


def executeWorkflowFromDB_inprocess(
    workflow_toolkit: hermesWorkflowToolkit,
    nameOrWorkflowFileOrJSONOrResource: str | dict,
    dispatch_id: str | None = None,
    workers: int = 1,
    targetTask: str = "finalnode_xx_0",
) -> str | None:
    """Build a saved workflow and run it in-process with ``luigi.build``.

    Prep (rebuild, build, write, clean target files) is delegated to the toolkit's
    ``prepareWorkflowRunFromDoc``; only the execution differs from
    ``hermesWorkflowToolkit.executeWorkflowFromDB``. Returns the ``dispatch_id``
    used for the run.

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
        # Steps 1-4: rebuild, build, write and clean target files (no execution).
        # Shared with the subprocess path so both build the module the same way.
        pythonFileName, workflowName, targetFilesDir = workflow_toolkit.prepareWorkflowRunFromDoc(doc)
        print(f"[in-process luigi.build] target files for {workflowName} at {targetFilesDir}", flush=True)

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
