"""
Tests for testing utilities
"""
import pathlib
import re
import shutil
import warnings

import pytest

from hera.datalayer import getDBObject


original_home = pathlib.Path.home()  # This gets messed up by pytester


@pytest.fixture()
def pytester_with_hera_config(pytester):
    # Required for successful logging initialization in pytester environment
    pytester.mkdir('.pyhera')
    pytester.mkdir('.pyhera/log')
    config_file = original_home / '.pyhera' / 'config.json'
    shutil.copy(config_file, pytester.path / '.pyhera' / 'config.json')
    # Activate the utilities we test here
    pytester.makeconftest("""
        from hera.utils.testing.mongo import *  # noqa
    """)
    # Now go
    return pytester


def test_session_hooks_prefix(pytester_with_hera_config):

    pytester = pytester_with_hera_config
    pytester.makepyfile("""
        from hera.utils.testing.mongo import HERA_TEST_DB_PREFIX
        from hera.datalayer import getDBObject
    
        def test_db_has_prefix():
            db_name = getDBObject('db_name')
            assert db_name.startswith(HERA_TEST_DB_PREFIX)    
    """)

    result = pytester.runpytest()

    result.assert_outcomes(passed=1)


@pytest.mark.parametrize("is_failure", [True, False])
def test_session_hooks_keep_db_only_on_fail(pytester_with_hera_config, is_failure):
    pytester = pytester_with_hera_config
    pytester.makepyfile("""
        import warnings
        from hera.utils.testing.mongo import HERA_TEST_DB_PREFIX
        from hera.datalayer import getDBObject

        def test_tell_db_prefix():
            conn = getDBObject("connection")
            db_name = getDBObject("db_name")
            conn[db_name]['test_collection'].insert_one({'This is a': 'test'})
            warnings.warn(f'HERA TEST DB PREFIX: {HERA_TEST_DB_PREFIX}')
            assert %(expr)s 
    """ % {'expr': not is_failure})  # We use old formatting because the contents uses f-strings and {}

    # Run in a subprocess so we have a separate pid and a separate DB
    result = pytester.runpytest_subprocess("-W", "ignore::DeprecationWarning")
    expect_outcome = "failed" if is_failure else "passed"
    result.assert_outcomes(**{expect_outcome: 1})
    output = result.stdout
    prefix_match = re.search(r"HERA TEST DB PREFIX: (hera@test@\S+)", str(output))
    assert prefix_match
    tested_prefix = prefix_match.group(1)
    # warnings.warn(f"Tested prefix: {tested_prefix}")
    conn = getDBObject('connection')
    databases = conn.list_database_names()
    assert is_failure == any(db.startswith(tested_prefix) for db in databases)
    # warnings.warn(f"Databases after non-removal test: {databases}")

    # Cleanup -- take some care that the prefix is real
    if is_failure:
        assert len(tested_prefix)>10
        from .mongo import drop_databases_by_prefix
        drop_databases_by_prefix(tested_prefix)
        # databases = conn.list_database_names()
        # warnings.warn(f"Databases after non-removal test cleanup: {databases}")


def test_fixture_cleans_database(pytester_with_hera_config):
    pytester = pytester_with_hera_config
    pytester.makepyfile("""
        # Add two tests, each asserting it starts clean and leaving something behind,
        # so that the whole suite fails no matter what order they run in, unless
        # the database is cleaned indeed

        from hera import datalayer
        from hera.utils.testing.mongo import HERA_TEST_DB_PREFIX
        
        def assert_only_system_collections(database):
            "Utility for checking that the database is clean"
            for col_name in database.list_collection_names():
                assert database[col_name].name.startswith("system."), f"Found {col_name} in db {db_name}"
    
        def check_db_starts_clean():
            "Utility for checking that the database stays clean"
            for user in datalayer.getDBNamesFromJSON():
                conn = datalayer.getDBObject("connection", user=user)
                db_name = datalayer.getDBObject("db_name", user=user)
                database = conn[db_name]
                # Make sure the db is empty
                assert_only_system_collections(database)
                # And then add a collection with data
                database['self_col'].insert_one({'This is a': 'test'})
    
        def test_db_starts_clean_1():
            check_db_starts_clean()
        def test_db_starts_clean_2():
            check_db_starts_clean()
    """)

    result = pytester.runpytest()
    result.assert_outcomes(passed=2)
