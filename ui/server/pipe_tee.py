import os
import threading


class PipeTee:
    """Captures everything written to a pipe so it can be sent to the UI.

    Hand ``write_fd`` to a child process (as its stdout/stderr). After the writer is
    done and ``close_write`` has been called, ``result`` returns the collected text.
    """

    def __init__(self):
        self._read_fd, self.write_fd = os.pipe()
        self._collected = []
        self._reader = threading.Thread(target=self._pump)
        self._reader.start()

    def _pump(self) -> None:
        while True:
            chunk = os.read(self._read_fd, 4096)
            if not chunk:
                break
            self._collected.append(chunk)

    def snapshot(self) -> str:
        """Return the text collected so far, without waiting for the run to end.

        Copies the buffer first (``list(...)``) so reading stays safe while the
        reader thread appends more chunks. No lock needed.
        """
        return b"".join(list(self._collected)).decode(errors="replace")

    def close_write(self) -> None:
        """Close our copy of the write end so the reader sees EOF once the child exits."""
        os.close(self.write_fd)

    def result(self) -> str:
        """Wait for the reader to drain the pipe and return the captured output."""
        self._reader.join()
        os.close(self._read_fd)
        return b"".join(self._collected).decode(errors="replace")
