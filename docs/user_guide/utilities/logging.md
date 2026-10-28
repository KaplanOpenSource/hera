# Logging Configuration

Hera includes a logging system built on Python's standard `logging` module. It uses a JSON-based configuration file stored in `~/.pyhera/log/` and provides helpers to customize log levels, add file handlers, and define formatters.

Logging is initialized automatically when Hera is imported. You only need to call `initialize_logging` explicitly if you want to customize the configuration.

The default configuration is copied on first use to `~/.pyhera/log/heraLogging.config` (a JSON `logging.config.dictConfig` file), and you can edit that copy. It defines the loggers `hera.datalayer`, `hera.measurements`, `hera.simulations`, `hera.bin` and `hera.riskassessment`, all at level `CRITICAL` and with `propagate: false`, each writing to the console and to a rotating log file in `~/.pyhera/log/` (`datalayer.log`, `measurements.log`, ...). The available formatters are `brief` and `default`.

!!! note "`disable_existing_loggers`"
    `initialize_logging(*overrides, disable_existing_loggers=True)` disables every logger that already exists by default. Since `import hera` already creates loggers, pass `disable_existing_loggers=False` when you call it yourself.

## Importing

```python
from hera.utils import initialize_logging, with_logger, get_logger, add_FileHandler, add_formatter
```

## Basic Setup

### Enable Debug Logging for a Module

```python
from hera.utils import initialize_logging, with_logger

initialize_logging(
    with_logger("hera.datalayer", level="DEBUG"),
    disable_existing_loggers=False,
)
```

This sets the `hera.datalayer` logger to DEBUG while keeping its configured handlers and leaving the other loggers at their default levels. A logger that is not in the default configuration (for example `hera.utils`) has no handlers, so also pass `handlers=["console"]` (or a file handler) to see its messages.

### Override Multiple Loggers

```python
initialize_logging(
    with_logger("hera.utils", level="DEBUG", handlers=["console"]),
    with_logger("hera.datalayer", level="WARNING"),
    disable_existing_loggers=False,
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

`get_logger(instance, name=None)` returns the class logger of `instance` (named `<module>.<ClassName>`) or, when `name` is given, `logging.getLogger(name)`. Related helpers: `getClassLogger(cls)` takes a class, and `get_classMethod_logger(instance, name)` returns `<module>.<ClassName>.<name>`.

## Writing Logs to a File

Use `add_FileHandler` to register a file handler, then assign it to a logger:

```python
from hera.utils import initialize_logging, with_logger, add_FileHandler

initialize_logging(
    add_FileHandler("my_file", "experiment.log", mode="w"),
    with_logger("hera.utils", level="DEBUG", handlers=["my_file"]),
    disable_existing_loggers=False,
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
    disable_existing_loggers=False,
)
```

The `add_formatter` function accepts:

- `formatterName` -- a name to reference the formatter.
- `format` -- a Python logging format string.
- `datefmt` -- date format (default: `"%Y-%m-%d %H:%M:%S"`).

Note that the default loggers (`hera.datalayer`, `hera.measurements`, ...) have `propagate: false`, so a handler attached to the parent `hera` logger does not receive their messages; attach it to those loggers directly.

## Log File Location

By default, Hera stores its log configuration and log files under:

```
~/.pyhera/log/
```

This directory is created automatically the first time logging is initialized (that is, on `import hera`). The path is fixed; it is not configurable through an environment variable.
