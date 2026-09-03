"""Captures this process's stdout/stderr and routes it two ways at once.

On ``start`` it redirects the real file descriptors 1 and 2 into an internal
pipe. A reader thread drains that pipe and, for every chunk of bytes:

    1. appends it to the current bucket in a ChunkState (so output is grouped by
       whichever task is running), and
    2. forwards it unchanged to ``forward_fd`` (the server pipe), so the live
       merged stream the UI already reads is unaffected.

It redirects the real fds, not just Python's ``sys.stdout``, so output from
shelled-out programs (OpenFOAM, the LSM binary) is captured too.

``stop`` restores the original fds, lets the reader drain, and joins it.
"""

import io
import os
import sys
import threading


class OutputRouter:
    # Set in start(); declared here as ints so the fd calls type-check cleanly.
    _forward: int
    _read_fd: int
    _write_fd: int
    _saved_out: int
    _saved_err: int
    _reader: threading.Thread

    def __init__(self, forward_fd: int, state):
        # Where to also send the captured output (the server pipe write end).
        self._forward_fd = forward_fd
        # The ChunkState whose `current` pointer decides the bucket per chunk.
        self._state = state

    def start(self):
        # Our own handle to the forward target, so the caller can close its copy.
        self._forward = os.dup(self._forward_fd)
        # The capture pipe: everything written to fd 1/2 lands here.
        self._read_fd, self._write_fd = os.pipe()
        # Save the real stdout/stderr so stop() can put them back.
        self._saved_out = os.dup(1)
        self._saved_err = os.dup(2)
        os.dup2(self._write_fd, 1)
        os.dup2(self._write_fd, 2)
        # Line-buffer stdout/stderr so each print reaches the pipe (and is bucketed)
        # while the pointer is still on the right task, instead of block-buffering
        # and flushing everything at the end into whatever task is current then.
        # reconfigure exists on TextIOWrapper (the concrete stdout/stderr type).
        if isinstance(sys.stdout, io.TextIOWrapper):
            sys.stdout.reconfigure(line_buffering=True)
        if isinstance(sys.stderr, io.TextIOWrapper):
            sys.stderr.reconfigure(line_buffering=True)
        self._reader = threading.Thread(target=self._pump, daemon=True)
        self._reader.start()

    def _pump(self):
        while True:
            chunk = os.read(self._read_fd, 4096)
            if not chunk:
                break
            # Route by the current pointer, then forward so live output is unchanged.
            self._state.append(chunk.decode(errors="replace"))
            os.write(self._forward, chunk)

    def stop(self):
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
        os.close(self._forward)
