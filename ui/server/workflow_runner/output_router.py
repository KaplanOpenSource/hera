"""Captures this process's stdout/stderr and streams it to the result queue.

On ``start`` it redirects the real file descriptors 1 and 2 into an internal pipe.
A reader thread drains that pipe and, for every chunk of bytes, puts a
``WorkflowOutput(task, text)`` on the result queue, tagged with the task pointer's
current task. The parent groups these into per-task chunks.

It redirects the real fds, not just Python's ``sys.stdout``, so output from
shelled-out programs (OpenFOAM, the LSM binary) is captured too.

``stop`` restores the original fds, lets the reader drain, and joins it.
"""

from __future__ import annotations

import io
import os
import sys
import threading
from multiprocessing.queues import Queue

from .task_pointer import TaskPointer
from .workflow_child_result import WorkflowMessage, WorkflowOutput


class OutputRouter:
    # Set in start(); declared here so the fd/thread calls type-check cleanly.
    _read_fd: int
    _write_fd: int
    _saved_out: int
    _saved_err: int
    _reader: threading.Thread

    def __init__(self, result_queue: Queue[WorkflowMessage], task_pointer: TaskPointer):
        # Where to send the captured output (the parent reads it off this queue).
        self._queue = result_queue
        # The TaskPointer whose `current` value tags each captured piece.
        self._task_pointer = task_pointer

    def start(self) -> None:
        # The capture pipe: everything written to fd 1/2 lands here.
        self._read_fd, self._write_fd = os.pipe()
        # Save the real stdout/stderr so stop() can put them back.
        self._saved_out = os.dup(1)
        self._saved_err = os.dup(2)
        os.dup2(self._write_fd, 1)
        os.dup2(self._write_fd, 2)
        # Line-buffer stdout/stderr so each print reaches the pipe (and is tagged)
        # while the pointer is still on the right task, instead of block-buffering
        # and flushing everything at the end into whatever task is current then.
        if isinstance(sys.stdout, io.TextIOWrapper):
            sys.stdout.reconfigure(line_buffering=True)
        if isinstance(sys.stderr, io.TextIOWrapper):
            sys.stderr.reconfigure(line_buffering=True)
        self._reader = threading.Thread(target=self._pump, daemon=True)
        self._reader.start()

    def _pump(self) -> None:
        while True:
            chunk = os.read(self._read_fd, 4096)
            if not chunk:
                break
            # Tag by the current pointer and hand the piece to the parent.
            self._queue.put(WorkflowOutput(task=self._task_pointer.current, text=chunk.decode(errors="replace")))

    def stop(self) -> None:
        sys.stdout.flush()
        sys.stderr.flush()
        # Restore the real fds. This drops the last writers on the capture pipe's
        # write end, so the reader hits EOF and exits.
        os.dup2(self._saved_out, 1)
        os.dup2(self._saved_err, 2)
        os.close(self._write_fd)
        self._reader.join()
        os.close(self._read_fd)
        os.close(self._saved_out)
        os.close(self._saved_err)
