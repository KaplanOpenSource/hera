# Logging Configuration

Hera includes a logging system built on Python's standard `logging` module. It uses a JSON-based configuration file stored in `~/.pyhera/log/` and provides helpers to customize log levels, add file handlers, and define formatters.

Logging is initialized automatically when Hera is imported. You only need to call `initialize_logging` explicitly if you want to customize the configuration.

## Importing

```python
from hera.utils import initialize_logging, with_logger, get_logger, add_FileHandler, add_formatter
```

## Basic Setup

### Enable Debug Logging for a Module

```python
from hera.utils import initialize_logging, with_logger

initialize_logging(
    with_logger("hera.utils", level="DEBUG"),
)
```

This sets the `hera.utils` logger to DEBUG while keeping other loggers at their default levels.

### Override Multiple Loggers

```python
initialize_logging(
    with_logger("hera.utils", level="DEBUG"),
    with_logger("hera.datalayer", level="WARNING"),
)
```

## Getting a Logger

Use `get_logger` to obtain a logger instance for your own code:

```python
from hera.utils import get_logger

# Named logger
logger = get_logger(None, "myproject.analysis")
logger.info("Starting analysis...")

# Class-scoped logger (pass an instance)
class MyProcessor:
    def run(self):
        logger = get_logger(self)  # logger name: mymodule.MyProcessor
        logger.debug("Processing started")
```

## Writing Logs to a File

Use `add_FileHandler` to register a file handler, then assign it to a logger:

```python
from hera.utils import initialize_logging, with_logger, add_FileHandler

initialize_logging(
    add_FileHandler("my_file", "experiment.log", mode="w"),
    with_logger("hera.utils", level="DEBUG", handlers=["my_file"]),
)
```

Parameters for `add_FileHandler`:

- `handlerName` -- a name to reference the handler.
- `fileName` -- path to the log file.
- `mode` -- `'w'` to overwrite (default) or `'a'` to append.
- `formatter` -- name of the formatter to use (default: `'default'`).

## Custom Formatters

Define a custom format string and attach it to a handler:

```python
from hera.utils import initialize_logging, add_formatter, add_FileHandler, with_logger

initialize_logging(
    add_formatter("detailed", format="%(asctime)s [%(levelname)s] %(name)s: %(message)s"),
    add_FileHandler("debug_file", "/tmp/hera_debug.log", mode="a", formatter="detailed"),
    with_logger("hera", level="DEBUG", handlers=["debug_file"]),
)
```

The `add_formatter` function accepts:

- `formatterName` -- a name to reference the formatter.
- `format` -- a Python logging format string.
- `datefmt` -- date format (default: `"%Y-%m-%d %H:%M:%S"`).

## Log File Location

By default, Hera stores its log configuration and log files under:

```
~/.pyhera/log/
```

This directory is created automatically the first time logging is initialized.
