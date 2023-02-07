import pytest

from hera.utils.testing.mongo import *  # noqa


pytest_plugins = ["pytester"]  # For use in testing the testing tools


@pytest.fixture()
def cached_data_regression(data_regression):
    """A utility for verifying a piece of data in Mongo stayed the same"""
    check = data_regression.check

    def check_cached(data):
        check(data.asDict())

    data_regression.check = check_cached
    return data_regression
