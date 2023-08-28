import json
import logging
import logging.config
import os.path
import pathlib
import pdb
from importlib.resources import read_text
from typing import List

# Special Hera logging level
EXECUTION = 15

HERA_DEFAULT_LOG_DIR = pathlib.Path.home() / ".pyhera" / "log"


# This function is named to match the style of the stdlib logging module
# noinspection PyPep8Naming
def getClassLogger(cls: type):
    name = cls.__module__ + "." + cls.__qualname__
    return logging.getLogger(name)


def get_logger(instance, name=None):
    return getClassLogger(instance.__class__) if name is None else logging.getLogger(name)

def get_classMethod_logger(instance, name=None):
    lgname = instance.__class__ if name is None else instance.__class__.__module__ + "." + instance.__class__.__qualname__ + "." + name
    return logging.getLogger(lgname)



def get_default_logging_config(*, disable_existing_loggers: bool = False) -> dict:
    defaultLocalConfig = os.path.join(HERA_DEFAULT_LOG_DIR,'heraLogging.config')
    if not os.path.isfile(defaultLocalConfig):
        defaultConfig = read_text('hera.utils.logging', 'heraLogging.config')
        with open(defaultLocalConfig,'w') as localConfig:
            localConfig.write(defaultConfig)

    with open(defaultLocalConfig,'r') as localConfig:
        config_text = "\n".join(localConfig.readlines())
    config_text = config_text.replace("{hera_log}", str(HERA_DEFAULT_LOG_DIR))
    config = json.loads(config_text)
    assert isinstance(config, dict)
    config['disable_existing_loggers'] = disable_existing_loggers
    return config


def _define_logger_execution():
    def execution(self, message, *args, **kws):
        self.log(EXECUTION, message, *args, **kws)

    logging.Logger.execution = execution
    logging.addLevelName(EXECUTION, 'EXECUTION')
    logging.EXECUTION = EXECUTION


def initialize_logging(*logger_overrides: (str, dict), disable_existing_loggers: bool = True) -> None:
    """
    Initialize logging for the Hera library

    Calling this function configures loggers for the hera library, according to
    the default configuration, with modifications you can add using :py:func:`with_logger` calls.

    The function is called automatically when hera is imported, so you only need
    to call it explicitly if you really want to add such modifications.

    To override e.g. the ``hera.bin`` logger's level to ``DEBUG``, call like so::

        initialize_logging(
            with_logger("hera.bin", level="DEBUG"),
        )

    This will use the rest of the definitions for this logger (handlers, formatters,
    whatever other parameters the logger class takes) from the default configuration,
    if such definitions exist.

    As the ``initialize`` in the name implies, you're expected to call this
    function, if at all, before you start getting logger objects. If you call it
    after some loggers were created, consider passing ``disable_existing_loggers=False``.
    """
    if not os.path.isdir(HERA_DEFAULT_LOG_DIR):
        os.makedirs(HERA_DEFAULT_LOG_DIR,exist_ok=False)
    _define_logger_execution()
    config = get_default_logging_config(disable_existing_loggers=disable_existing_loggers)
    for logger_name, logger_dict in logger_overrides:
        # This says: Use whatever was configured, if any, and update with what was provided
        config['loggers'].setdefault(logger_name, logger_dict).update(logger_dict)
    logging.config.dictConfig(config)


def with_logger(logger_name, level=None, handlers=None, propagate=None) -> (str, dict):
    """Build a dictionary describing a logger, for use with initialize_logging()"""
    logger_dict = dict(level=level, handlers=handlers, propagate=propagate)
    # Remove from dict parameters not supplied; this allows the use here
    # to just override specific settings on existing loggers
    empty = [key for key, value in logger_dict.items() if value is None]
    for key in empty:
        del logger_dict[key]

    return logger_name, logger_dict
