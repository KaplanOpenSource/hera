import os
import sys
import threading


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
        The output is duplicated: it still shows on the server console AND is
        collected into the returned string. The toolkit's files directory is
        temporarily added to PYTHONPATH so the generated Luigi module is importable
        by the subprocess. Both fds and PYTHONPATH are restored afterwards.
        """
        from hera import toolkitHome

        workflow_toolkit = toolkitHome.getToolkit(
            toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS,
            projectName=project_name,
        )

        previous_pythonpath = os.environ.get("PYTHONPATH")
        os.environ["PYTHONPATH"] = workflow_toolkit.FilesDirectory + os.pathsep + (previous_pythonpath or "")

        # Everything written to fds 1/2 goes into the pipe; a reader thread echoes
        # it back to the real console (saved_stdout_fd) and also collects it.
        saved_stdout_fd = os.dup(1)
        saved_stderr_fd = os.dup(2)
        pipe_read_fd, pipe_write_fd = os.pipe()

        collected = []

        def tee():
            while True:
                chunk = os.read(pipe_read_fd, 4096)
                if not chunk:
                    break
                collected.append(chunk)
                os.write(saved_stdout_fd, chunk)

        reader = threading.Thread(target=tee)
        reader.start()

        sys.stdout.flush()
        sys.stderr.flush()
        os.dup2(pipe_write_fd, 1)
        os.dup2(pipe_write_fd, 2)

        try:
            dispatch_id = workflow_toolkit.executeWorkflowFromDB(workflow_name, scheduler="local")
        finally:
            sys.stdout.flush()
            sys.stderr.flush()
            os.dup2(saved_stdout_fd, 1)
            os.dup2(saved_stderr_fd, 2)
            os.close(pipe_write_fd)
            reader.join()
            os.close(pipe_read_fd)
            os.close(saved_stdout_fd)
            os.close(saved_stderr_fd)
            if previous_pythonpath is None:
                os.environ.pop("PYTHONPATH", None)
            else:
                os.environ["PYTHONPATH"] = previous_pythonpath

        output = b"".join(collected).decode(errors="replace")

        return {"dispatch_id": dispatch_id, "output": output}
