"""Shared 'current chunk' pointer for routing captured workflow output per task.

The idea: instead of parsing the merged log afterwards, we keep a pointer that
says which task is running right now. The Luigi event handlers set the pointer
(START -> the task name, SUCCESS/FAILURE -> "between"); the output router reads
the pointer and drops each captured piece of stdout/stderr into that task's
bucket. One instance lives per workflow child process.

Sequential runs only: with more than one Luigi worker, tasks run in separate
forked processes that each get their own copy of this pointer, so the child's
router would no longer see their updates. Parallel needs per-task capture instead.
"""

# Bucket names for output that does not belong to a specific task.
PREAMBLE = "__preamble__"  # before the first task starts
BETWEEN = "__between__"  # while no task is running


class ChunkState:
    """Holds the current chunk name and an ordered list of output segments.

    Segments are ordered, not keyed by name, so each gap keeps its place: the
    output between task A and task B and the output between task B and task C stay
    as two separate "between" segments in run order, instead of merging into one.
    """

    def __init__(self):
        self.current = PREAMBLE
        # Each entry is [name, list-of-text-pieces], appended in run order.
        self._segments = []

    def append(self, text):
        """Add captured output to the current segment, starting a new one when the
        pointer has moved since the last piece."""
        if not self._segments or self._segments[-1][0] != self.current:
            self._segments.append([self.current, []])
        self._segments[-1][1].append(text)

    def as_list(self):
        """Return [{'name', 'text'}, ...] in run order."""
        result = []
        for name, parts in self._segments:
            result.append({"name": name, "text": "".join(parts)})
        return result


# The single shared pointer/buckets for this process.
state = ChunkState()
