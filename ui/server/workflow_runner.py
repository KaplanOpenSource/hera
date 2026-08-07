import os
import sys
import threading


class WorkflowRunner:
    """Builds and executes saved Hermes workflows, capturing their console output.

    Instantiated once and reused (it will hold run state/config later).
    """

    def __init__(self):
        # Serialize runs: fd 1/2 are process-wide, so concurrent runs would cross-capture.
        self._lock = threading.Lock()

    def run(self, project_name: str, workflow_name: str) -> dict:
        """Build and execute a saved workflow from the DB (local Luigi scheduler).

        Returns ``{"dispatch_id", "output"}``. Output is captured at the fd level
        (Luigi subprocess / os.system write to fds, not sys.stdout) and duplicated
        to the console. The files directory is added to PYTHONPATH so the generated
        module is importable. Both fds and PYTHONPATH are restored afterwards.
        """
        from hera import toolkitHome

        with self._lock:
            workflow_toolkit = toolkitHome.getToolkit(
                toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS,
                projectName=project_name,
            )

            previous_pythonpath = os.environ.get("PYTHONPATH")
            os.environ["PYTHONPATH"] = workflow_toolkit.FilesDirectory + os.pathsep + (previous_pythonpath or "")

            # Tee fds 1/2 through a pipe: reader thread echoes to console and collects.
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
