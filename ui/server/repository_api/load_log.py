"""Everything hera complains about while a repository loads."""
from __future__ import annotations

import logging
import traceback
from typing import Dict, List

# hera skips a toolkit it cannot build with this line, at info level, and goes on.
SKIPPED = "Skipping key"


class LoadLog(logging.Handler):
    """Listens to every hera logger for the length of one load.

    hera's loader does not raise: a section that fails is logged and the load
    carries on. Worse, hera's own logging config mutes whole families
    (`hera.datalayer`, `hera.measurements`, `hera.riskassessment`) at CRITICAL
    with propagation off, so their complaints reach nothing at all. This opens
    every hera logger and listens on each, then puts the levels back.
    """

    def __init__(self):
        super().__init__(level=logging.INFO)
        self.problems: List[str] = []
        # The level each logger had before, so the load changes nothing lasting.
        self._levels: Dict[str, int] = {}

    def emit(self, record: logging.LogRecord) -> None:
        # The handler sits on every hera logger, so one record arrives again for
        # each logger it propagates through. Mark it and keep only the first.
        if getattr(record, "seenByLoadLog", False):
            return
        record.seenByLoadLog = True
        message = record.getMessage()
        if record.levelno < logging.WARNING and SKIPPED not in message:
            return
        self.problems.append(f"{record.name}: {message}")
        # The exception is the message that matters; never drop it.
        if record.exc_info:
            self.problems.append("".join(traceback.format_exception(*record.exc_info)).strip())

    @staticmethod
    def _heraLoggers() -> List[logging.Logger]:
        """Every hera logger that exists now, the root one of them included."""
        names = ["hera"] + [name for name in logging.Logger.manager.loggerDict
                            if name.startswith("hera.")]
        return [logging.getLogger(name) for name in names]

    def __enter__(self) -> "LoadLog":
        # A logger born during the load inherits the level from these and
        # propagates up to them, so it is heard too.
        for logger in self._heraLoggers():
            self._levels[logger.name] = logger.level
            logger.setLevel(logging.INFO)
            logger.addHandler(self)
        return self

    def __exit__(self, *error) -> None:
        """Puts every level back. Answers nothing, so no error is ever suppressed."""
        for name, level in self._levels.items():
            logger = logging.getLogger(name)
            logger.removeHandler(self)
            logger.setLevel(level)
