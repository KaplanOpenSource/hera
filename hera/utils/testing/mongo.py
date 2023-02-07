import importlib
from os import getpid

import pytest

from hera import datalayer

# Get a relatively unique name prefix
HERA_TEST_DB_PREFIX = f'hera@test@{getpid()}_'


db_configs = {}


def pytest_sessionstart(session):
    """Called before testing starts"""
    from hera.datalayer.document import getMongoJSON
    db_configs.update(getMongoJSON())
    # Override the created databases; this uses ours instead
    for user, db_dict in db_configs.items():
        # DEBUG warnings.warn(f"Creating test db for {db_dict['dbName']}")
        db_dict["dbName"] = HERA_TEST_DB_PREFIX + db_dict["dbName"]
        datalayer.createDBConnection(user=user, mongoConfig=db_dict)
        # DEBUG warnings.warn(f"Created db {db_dict['dbName']}")

    importlib.reload(datalayer)  # Some initializations there need to be redone too


def pytest_sessionfinish(session, exitstatus):
    """When testing ends successfully, we delete the DB; otherwise, keep it"""
    # exitstatus is 0 (falsey) when everything is ok.
    if not exitstatus:
        for user, db_dict in db_configs.items():
            db_conn = datalayer.getDBObject("connection", user=user)
            db_name = db_dict["dbName"]
            assert db_name.startswith(HERA_TEST_DB_PREFIX)
            db_conn.drop_database(db_name)


@pytest.fixture(autouse=True)
def clean_mongo():
    """Clean all user databases before the test starts"""
    user_dbs = {}
    for user in datalayer.getDBNamesFromJSON():
        user_dbs[user] = datalayer.getDBObject("connection", user=user)

    # Strongly inspired by https://github.com/ClearcodeHQ/pytest-mongo/blob/v2.1.1/src/pytest_mongo/factories.py
    for db_conn in set(user_dbs.values()):
        for db_name in db_conn.list_database_names():
            if db_name.startswith(HERA_TEST_DB_PREFIX):
                database = db_conn[db_name]
                for collection_name in database.list_collection_names():
                    collection = database[collection_name]
                    # Do not delete any of Mongo "system" collections
                    # Note collection_name above but collection.name below
                    if not collection.name.startswith("system."):
                        # warnings.warn(f"Removing collection {collection.name} from db {db_name}")
                        collection.drop()


# For manual use, in case test databases have accumulated
def drop_databases_by_prefix(prefix):
    for user in datalayer.getDBNamesFromJSON():
        db_conn = datalayer.getDBObject("connection", user=user)
        for db_name in db_conn.list_database_names():
            if db_name.startswith(prefix):
                db_conn.drop_database(db_name)


__all__ = ["pytest_sessionstart", "pytest_sessionfinish", "clean_mongo"]
