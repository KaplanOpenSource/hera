import logging
import os
import pathlib

import pytest

from .helpers import EXECUTION, getClassLogger, get_default_logging_config, initialize_logging, with_logger

"""
A decorator for tests not to be run by default, because they are likely to mess up logging.
To run them, set the HERA_FULL_LOGGING_TESTS environment variable, e.g.

$ HERA_FULL_LOGGING_TESTS=1 pytest
"""
only_when_full_logging_tests = pytest.mark.skipif(
    os.environ.get('HERA_FULL_LOGGING_TESTS', None) is None,
    reason="HERA_FULL_LOGGING_TESTS not set in environment"
)


class TestClass:
    class InnerTestClass:
        pass


def test_getclasslogger_external():
    logger = getClassLogger(TestClass)
    name_parts = logger.name.rsplit('.', 1)
    assert name_parts[0] == TestClass.__module__
    assert name_parts[-1] == 'TestClass'


def test_getclasslogger_nested():
    logger = getClassLogger(TestClass.InnerTestClass)
    name_parts = logger.name.rsplit('.', 2)
    assert name_parts[0] == TestClass.__module__
    assert name_parts[-2] == 'TestClass'
    assert name_parts[-1] == 'InnerTestClass'


def test_getclasslogger_infunc():

    class InFunc:
        pass

    logger = getClassLogger(InFunc)
    name_parts = logger.name.rsplit('.', 1)
    assert name_parts[-1] == 'InFunc'


def test_get_default_logging_config():
    config = get_default_logging_config()
    assert isinstance(config, dict)
    assert not config.get('disable_existing_loggers', True)
    # Sample configs relying on the '{hera_log}' replacement
    _assert_resolved_path(config['handlers']['bin']['filename'])
    _assert_resolved_path(config['handlers']['simulations.old']['filename'])


def _assert_resolved_path(alleged_path):
    assert isinstance(alleged_path, str)
    assert '{hera_log}' not in alleged_path
    p = pathlib.Path(alleged_path)
    assert p.is_absolute()


def test_with_logger():
    name, shunra = with_logger("shunra", level="EXECUTION", handlers=["bin"], propagate=True)
    assert name == "shunra"
    assert shunra["level"] == "EXECUTION"
    assert shunra["handlers"] == ["bin"]
    assert shunra["propagate"] is True


def test_with_logger_defaults():
    name, zeus = with_logger("zeus", level="DEBUG")
    assert name == "zeus"
    assert zeus["level"] == "DEBUG"
    assert "handlers" not in zeus
    assert "propagate" not in zeus


@only_when_full_logging_tests
def test_initialize_logging_with_logger():
    """
    This test is not run by default because it is likely to mess up logging.
    To run it, set the HERA_FULL_LOGGING_TESTS
    """
    config = get_default_logging_config()
    initialize_logging(
        with_logger("bin", "DEBUG"),
        with_logger("sin", "EXECUTION", propagate=False, handlers=["console"]),
        disable_existing_loggers=True
    )
    bin = logging.getLogger("bin")
    assert config['loggers']['bin']['level'] == "CRITICAL"
    assert bin.level == logging.DEBUG  # Overridden
    assert bin.propagate is False and len(bin.handlers) == 2  # from config

    sin = logging.getLogger("sin")
    assert "sin" not in config["loggers"]
    assert sin.level == EXECUTION and sin.propagate is False
    assert isinstance(sin.handlers[0], logging.StreamHandler)

