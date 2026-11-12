from __future__ import annotations

import threading

from .workflow_child_result import WorkflowOutput


class WorkflowLogBuilder:
    """Turns the child's tagged output messages into per-task chunks.

    Thread-safe: the queue reader (in ``run``) calls ``add`` while ``poll`` reads
    a ``chunks`` snapshot at the same time.

    Segments are ordered, not keyed by name, so each gap keeps its place: the
    output between task A and task B and the output between task B and task C stay
    as two separate "between" segments in run order, instead of merging into one.
    """

    def __init__(self):
        self._lock = threading.Lock()
        self._segments = []  # ordered [name, [text pieces]] entries

    def add(self, message: WorkflowOutput) -> None:
        """Add one output message to its task's segment."""
        with self._lock:
            if not self._segments or self._segments[-1][0] != message.task:
                self._segments.append([message.task, []])
            self._segments[-1][1].append(message.text)

    def chunks(self) -> list:
        """Per-task segments as ``[{'name', 'text'}, ...]`` in run order."""
        with self._lock:
            return [{"name": name, "text": "".join(parts)} for name, parts in self._segments]
