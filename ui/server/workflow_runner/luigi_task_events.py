"""Luigi task event handlers for in-process workflow runs.

Importing this module registers the handlers for every ``luigi.Task`` run in this
process (the ``event_handler`` decorator registers each on class-body execution).
They print ``[luigi-event]`` lines, and START/SUCCESS/FAILURE move the task pointer
so the output router tags each task's output with the task's name.
"""

import luigi

from .task_pointer import BETWEEN, task_pointer


# Every Luigi event line starts with this prefix so the UI log parser can spot
# them and hide them by default, without confusing them for a task's own output.
# Kept in sync with EVENT_PREFIX in the client's classifyLog.ts.
EVENT_PREFIX = "[luigi-event]"

# Task failures seen during the current run, as (task_family, message). The child
# clears this before luigi.build and reads it after, so a failed run can surface the
# real task error instead of luigi.build's bare False.
task_failures = []


def _event(message):
    """Print one Luigi-event line with the shared prefix (flushed for live output)."""
    print(f"{EVENT_PREFIX} {message}", flush=True)


class LuigiTaskEvents:
    """Luigi event handlers, registered at import for every task in this process.

    Each handler is registered when the class body runs (the event_handler decorator
    registers the raw function). ``@staticmethod`` is the outer decorator so Luigi
    calls the plain ``fn(task, ...)`` with no ``self``.
    """

    @staticmethod
    @luigi.Task.event_handler(luigi.Event.START)
    def on_start(task):
        # Point at this task first, so the START line and the run's output land in its bucket.
        task_pointer.current = task.task_family
        _event(f"START {task.task_family}")

    @staticmethod
    @luigi.Task.event_handler(luigi.Event.SUCCESS)
    def on_success(task):
        # Print while still on the task's bucket, then stop pointing at it.
        _event(f"SUCCESS {task.task_family}")
        task_pointer.current = BETWEEN

    @staticmethod
    @luigi.Task.event_handler(luigi.Event.FAILURE)
    def on_failure(task, exception):
        _event(f"FAILURE {task.task_family}: {exception}")
        task_failures.append((task.task_family, str(exception)))
        task_pointer.current = BETWEEN

    @staticmethod
    @luigi.Task.event_handler(luigi.Event.PROCESSING_TIME)
    def on_time(task, seconds):
        _event(f"TIME {task.task_family}: {seconds:.2f}s")

    @staticmethod
    @luigi.Task.event_handler(luigi.Event.BROKEN_TASK)
    def on_broken(task, exception):
        _event(f"BROKEN {task.task_family}: {exception}")
        task_failures.append((task.task_family, str(exception)))
        task_pointer.current = BETWEEN

    @staticmethod
    @luigi.Task.event_handler(luigi.Event.PROGRESS)
    def on_progress(task, progress):
        _event(f"PROGRESS {task.task_family}: {progress}")

    @staticmethod
    @luigi.Task.event_handler(luigi.Event.DEPENDENCY_DISCOVERED)
    def on_dependency_discovered(task, dependency):
        _event(f"DEP DISCOVERED {task.task_family} -> {dependency.task_family}")

    @staticmethod
    @luigi.Task.event_handler(luigi.Event.DEPENDENCY_MISSING)
    def on_dependency_missing(task):
        _event(f"DEP MISSING {task.task_family}")

    @staticmethod
    @luigi.Task.event_handler(luigi.Event.DEPENDENCY_PRESENT)
    def on_dependency_present(task):
        _event(f"DEP PRESENT {task.task_family}")
