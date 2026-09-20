import os
import shutil
import sys
import types

import pytest

# Make the server modules (server.py, api_models.py, ...) importable regardless of
# where pytest is run from.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# server.py calls argparse.parse_args() at import time; give it clean args so it
# doesn't try to parse pytest's argv.
sys.argv = ["server"]


# Source of the generated Luigi module the fake toolkit writes to disk. The child
# imports it and runs finalnode_xx_0 with luigi.build, exactly as it would for a real
# workflow. The task calls the test's on_execute hook (via a module the child inherits
# through fork), then writes its output target so luigi sees it as complete.
_GENERATED_TASK_TEMPLATE = '''\
import os
import luigi
import workflow_test_hooks

WORKFLOW_NAME = {name!r}
TARGET_DIR = {target_dir!r}


class finalnode_xx_0(luigi.Task):
    def output(self):
        return luigi.LocalTarget(os.path.join(TARGET_DIR, "done"))

    def run(self):
        os.makedirs(TARGET_DIR, exist_ok=True)
        callback = getattr(workflow_test_hooks, "on_execute", None)
        if callback is not None:
            callback(WORKFLOW_NAME)
        with self.output().open("w") as handle:
            handle.write("done")
'''


@pytest.fixture
def install_fake_hera(monkeypatch):
    """Install fakes so the runner works without real hera (but with real luigi).

    Returns a factory: ``install_fake_hera(files_directory, on_execute=None)``.

    It fakes two things the workflow child uses:
      - ``hera.toolkitHome`` -> a toolkit exposing ``FilesDirectory`` and
        ``prepareWorkflowRunFromDoc``. The fake prepare writes a real Luigi task
        module (``finalnode_xx_0``) to the files directory, so the child runs it with
        ``luigi.build`` just like a real workflow.
      - ``workflow_test_hooks.on_execute`` -> the test's callback. The generated task
        calls it inside its ``run``, so a test can produce captured output or raise.

    The child runs in a forked process, so anything ``on_execute`` records in memory
    is NOT visible to the parent. Tests must assert on the run's observable result
    instead: the returned chunks / dispatch id, or what ``poll`` reports.
    ``on_execute`` can ``os.write(1/2, ...)`` to produce captured output, and raising
    from it fails the run (the failure surfaces in the run's error / chunks).
    """
    def _install(files_directory, on_execute=None):
        def prepare_workflow_run_from_doc(doc):
            workflow_name = doc.desc["workflowName"]
            target_dir = os.path.join(files_directory, f"{workflow_name}_targetFiles")
            # Match the real prep: clear prior target files so the task re-runs.
            shutil.rmtree(target_dir, ignore_errors=True)
            python_file = os.path.join(files_directory, f"{workflow_name}.py")
            with open(python_file, "w") as handle:
                handle.write(_GENERATED_TASK_TEMPLATE.format(name=workflow_name, target_dir=target_dir))
            return python_file, workflow_name, target_dir

        toolkit = types.SimpleNamespace(
            FilesDirectory=files_directory,
            prepareWorkflowRunFromDoc=prepare_workflow_run_from_doc,
        )

        def get_toolkit(toolkitName, projectName):
            return toolkit

        toolkit_home = types.SimpleNamespace(
            SIMULATIONS_WORKFLOWS="SIMULATIONS_WORKFLOWS",
            getToolkit=get_toolkit,
        )
        fake_hera = types.ModuleType("hera")
        fake_hera.toolkitHome = toolkit_home
        monkeypatch.setitem(sys.modules, "hera", fake_hera)

        # The generated task imports this to reach the test's callback; the forked
        # child inherits sys.modules, so installing it here makes it importable there.
        hooks = types.ModuleType("workflow_test_hooks")
        hooks.on_execute = on_execute
        monkeypatch.setitem(sys.modules, "workflow_test_hooks", hooks)

    return _install
