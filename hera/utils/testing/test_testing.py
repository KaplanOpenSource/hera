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
    return pytester


def test_session_hooks_prefix(pytester_with_hera_config):

    pytester = pytester_with_hera_config
    pytester.makeconftest("""
        from hera.utils.testing.mongo import *  # noqa
    """)
    pytester.makepyfile("""
        from hera.utils.testing.mongo import HERA_TEST_DB_PREFIX
        from hera.datalayer import getDBObject
    
        def test_db_has_prefix():
            db_name = getDBObject('db_name')
            assert db_name.startswith(HERA_TEST_DB_PREFIX)
    
    """)

    result = pytester.runpytest()

    result.assert_outcomes(passed=1)


# TODO: add test-session-leaves-db-in-place-when-tests-fail
# TODO: and while at it, consider code duplications
def test_session_hooks_remove_db_on_success(pytester_with_hera_config):
    pytester = pytester_with_hera_config
    pytester.makeconftest("""
       from hera.utils.testing.mongo import *  # noqa
    """)
    pytester.makepyfile("""
        import warnings
        from hera.utils.testing.mongo import HERA_TEST_DB_PREFIX
        from hera.datalayer import getDBObject

        def test_tell_db_prefix():
            conn = getDBObject("connection")
            db_name = getDBObject("db_name")
            conn[db_name]['test_collection'].insert_one({'This is a': 'test'})
            warnings.warn(f'HERA TEST DB PREFIX: {HERA_TEST_DB_PREFIX}') 
    """)

    # Run in a subprocess so we have a separate pid and a separate DB
    result = pytester.runpytest_subprocess("-W", "ignore::DeprecationWarning")
    result.assert_outcomes(passed=1)
    output = result.stdout
    prefix_match = re.search(r"HERA TEST DB PREFIX: (hera@test@\S+)", str(output))
    assert prefix_match
    tested_prefix = prefix_match.group(1)
    # warnings.warn(f"Tested prefix: {tested_prefix}")
    conn = getDBObject('connection')
    databases = conn.list_database_names()
    assert not any(db.startswith(tested_prefix) for db in databases)
    # warnings.warn(f"Databases after removal test: {databases}")


def test_session_hooks_dont_remove_db_on_failure(pytester_with_hera_config):
    pytester = pytester_with_hera_config
    pytester.makeconftest("""
       from hera.utils.testing.mongo import *
    """)
    pytester.makepyfile("""
        import warnings
        from hera.utils.testing.mongo import HERA_TEST_DB_PREFIX
        from hera.datalayer import getDBObject

        def test_tell_db_prefix():
            conn = getDBObject("connection")
            db_name = getDBObject("db_name")
            conn[db_name]['test_collection'].insert_one({'This is a': 'test'})
            warnings.warn(f'HERA TEST DB PREFIX: {HERA_TEST_DB_PREFIX}')
            assert False, "Failing internally for the bigger test" 
    """)

    # Run in a subprocess so we have a separate pid and a separate DB
    result = pytester.runpytest_subprocess("-W", "ignore::DeprecationWarning")
    result.assert_outcomes(failed=1)
    output = result.stdout
    prefix_match = re.search(r"HERA TEST DB PREFIX: (hera@test@\S+)", str(output))
    assert prefix_match
    tested_prefix = prefix_match.group(1)
    # warnings.warn(f"Tested prefix: {tested_prefix}")
    conn = getDBObject('connection')
    databases = conn.list_database_names()
    assert any(db.startswith(tested_prefix) for db in databases)
    # warnings.warn(f"Databases after non-removal test: {databases}")

    # Cleanup -- take some care that the prefix is real
    assert len(tested_prefix)>10
    from .mongo import drop_databases_by_prefix
    drop_databases_by_prefix(tested_prefix)
    databases = conn.list_database_names()
    # warnings.warn(f"Databases after non-removal test cleanup: {databases}")


def test_fixture_cleans_database(pytester_with_hera_config):
    pytester = pytester_with_hera_config
    pytester.makeconftest("""
        from hera.utils.testing.mongo import *  # noqa
    """)
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
    print("Test output")
    print(result.stdout)
    print("Test error")
    print(result.stderr)

    result.assert_outcomes(passed=2)



