from enum import Enum


class RunStatus(str, Enum):
    """Run status reported by poll(); BUSY is returned by start().

    Subclasses ``str`` so members serialize to their plain string value over the
    wire (e.g. ``"running"``) and compare equal to it.
    """

    IDLE = "idle"  # no run started yet on this runner
    RUNNING = "running"
    DONE = "done"
    ERROR = "error"
    BUSY = "busy"
    NOT_FOUND = "not_found"
