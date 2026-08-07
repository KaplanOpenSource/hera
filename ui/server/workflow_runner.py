import os
import sys
import tempfile


class WorkflowRunner:
    """Builds and executes saved Hermes workflows, capturing their console output.

    Instantiated once and reused (it will hold run state/config later).
    """

    def run(self, project_name: str, workflow_name: str) -> dict:
        """Build and execute a saved workflow from the DB (local Luigi scheduler).

        Returns ``{"dispatch_id", "output"}`` where output is the run's captured
        stdout+stderr.

        stdout+stderr are captured at the fd level because the Luigi subprocess and
        the node's os.system write to the process file descriptors, not sys.stdout.
        The toolkit's files directory is temporarily added to PYTHONPATH so the
        generated Luigi module is importable by the subprocess. Both are restored
        afterwards.
        """
        from hera import toolkitHome

        workflow_toolkit = toolkitHome.getToolkit(
            toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS,
            projectName=project_name,
        )

        previous_pythonpath = os.environ.get("PYTHONPATH")
        os.environ["PYTHONPATH"] = workflow_toolkit.FilesDirectory + os.pathsep + (previous_pythonpath or "")

        log_fd, log_path = tempfile.mkstemp(suffix=".log")
        saved_stdout_fd = os.dup(1)
        saved_stderr_fd = os.dup(2)
        sys.stdout.flush()
        sys.stderr.flush()
        os.dup2(log_fd, 1)
        os.dup2(log_fd, 2)

        try:
            dispatch_id = workflow_toolkit.executeWorkflowFromDB(workflow_name, scheduler="local")
        finally:
            sys.stdout.flush()
            sys.stderr.flush()
            os.dup2(saved_stdout_fd, 1)
            os.dup2(saved_stderr_fd, 2)
            os.close(saved_stdout_fd)
            os.close(saved_stderr_fd)
            os.close(log_fd)
            if previous_pythonpath is None:
                os.environ.pop("PYTHONPATH", None)
            else:
                os.environ["PYTHONPATH"] = previous_pythonpath

        with open(log_path) as log_file:
            output = log_file.read()
        os.remove(log_path)

        return {"dispatch_id": dispatch_id, "output": output}
