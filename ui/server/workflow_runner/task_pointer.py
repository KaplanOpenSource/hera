"""Shared 'current task' pointer for tagging workflow output.

Luigi event handlers set it (START -> task name, SUCCESS/FAILURE -> BETWEEN); the
output router reads it to tag each captured output piece with its task. One per
child process; sequential runs only (parallel workers fork and lose the updates).
"""

# Pointer values for output that does not belong to a specific task.
PREAMBLE = "__preamble__"  # before the first task starts
BETWEEN = "__between__"  # while no task is running


class TaskPointer:
    """Holds the name of the task running right now, for the router to read."""

    def __init__(self):
        self.current = PREAMBLE


# The single shared pointer for this process.
task_pointer = TaskPointer()
