import importlib
import warnings
from os import getpid

import pytest

from hera import datalayer

# Get a relatively unique name prefix
HERA_TEST_DB_PREFIX = f'hera@test@{getpid()}_'


@pytest.fixture(autouse=True, scope="session")
def use_test_database():
    from hera.datalayer.document import getMongoJSON
    db_configs = getMongoJSON()
    # Override the created databases; this uses ours instead
    for user, db_dict in db_configs.items():
        # DEBUG warnings.warn(f"Creating test db for {db_dict['dbName']}")
        db_dict["dbName"] = HERA_TEST_DB_PREFIX + db_dict["dbName"]
        datalayer.createDBConnection(user=user, mongoConfig=db_dict)
        # DEBUG warnings.warn(f"Created db {db_dict['dbName']}")

    importlib.reload(datalayer)  # Some initializations there need to be redone too

    yield

    for user, db_dict in db_configs.items():
        db_conn = datalayer.getDBObject("connection", user=user)
        db_name = db_dict["dbName"]
        assert db_name.startswith(HERA_TEST_DB_PREFIX)
        db_conn.drop_database(db_name)


@pytest.fixture(autouse=True)
def clean_mongo():

    user_dbs = {}
    for user in datalayer.getDBNamesFromJSON():
        user_dbs[user] = datalayer.getDBObject("connection", user=user)

    yield user_dbs

    for db_conn in set(user_dbs.values()):
        for db_name in db_conn.list_database_names():
            if db_name.startswith(HERA_TEST_DB_PREFIX):
                database = db_conn[db_name]
                for collection_name in database.list_collection_names():
                    collection = database[collection_name]
                    # Do not delete any of Mongo "system" collections
                    # Note collection_name above but collection.name below
                    if not collection.name.startswith("system."):
                        warnings.warn(f"Removing collection {collection.name} from db {db_name}")
                        collection.drop()


@pytest.fixture()
def cached_data_regression(data_regression):

    check = data_regression.check

    def check_cached(data):
        check(data.asDict())

    data_regression.check = check_cached
    return data_regression
