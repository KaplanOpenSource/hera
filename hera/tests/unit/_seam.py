"""Test-only seam that points hera's datalayer at an in-memory MongoDB.

Why a seam and not mocks: ``abstractToolkit`` inherits from ``Project``
(hera/toolkit.py:33), and ``Project.__init__`` reads from the DB, writes to
the DB and creates a directory under ``~/.hera``.  No toolkit can be built
without a database.  Rather than mock every call site, we substitute the one
function that creates the connection, so the real datalayer code — query
building, desc filtering, version resolution, document tagging — runs for
real against an in-memory store.

The substituted function is ``hera.datalayer.document.connectToDatabase``.
``createDBConnection`` calls it by module-global name, so rebinding the
attribute is enough.
"""
import getpass
import json
import os
import pathlib

UNIT_DB_NAME = "heraUnitTest"
UNIT_DB_ALIAS = f"{UNIT_DB_NAME}-alias"

_UNIT_MONGO_CONFIG = dict(
    dbIP="127.0.0.1",
    dbName=UNIT_DB_NAME,
    username="unit",
    password="unit",
)


def write_pyhera_config():
    """Write a placeholder ~/.pyhera/config.json under the current HOME.

    MUST run before the first ``import hera.datalayer``:
    hera/datalayer/__init__.py:7-10 builds collection singletons at import
    time, and ``getDBObject`` raises ``KeyError: user <name> not found``
    when the config has no entry for the current user.

    The placeholder dbIP is never dialled — ``mongoengine.connect()`` is
    lazy, and install() replaces the connection before any query runs.

    The caller is responsible for pointing HOME at a temp directory first;
    see the bootstrap block in conftest.py for why that cannot live here.

    Returns the path of the file written.
    """
    pyhera = pathlib.Path(os.environ["HOME"], ".pyhera")
    pyhera.mkdir(parents=True, exist_ok=True)
    config_path = pyhera / "config.json"
    with open(config_path, "w", encoding="utf-8") as handle:
        json.dump({getpass.getuser(): dict(_UNIT_MONGO_CONFIG)}, handle)
    return config_path


def install():
    """Replace connectToDatabase with a mongomock variant and re-register."""
    import mongomock
    from mongoengine import connect, disconnect

    import hera.datalayer.document as hdoc

    def _mongomock_connect(mongoConfig, alias=None):
        if isinstance(mongoConfig, str):
            mongoConfig = hdoc.parseConnectionString(mongoConfig)
        alias = "%s-alias" % mongoConfig["dbName"] if alias is None else alias
        disconnect(alias)
        return connect(
            alias=alias,
            db=mongoConfig["dbName"],
            mongo_client_class=mongomock.MongoClient,
            uuidRepresentation="standard",
        )

    hdoc.connectToDatabase = _mongomock_connect
    hdoc.createDBConnection(
        connectionName=getpass.getuser(),
        mongoConfig=dict(_UNIT_MONGO_CONFIG),
    )
    _rebind_module_level_collections()


def _rebind_module_level_collections():
    """Rebuild the singletons bound at import time, before the seam existed.

    hera/datalayer/__init__.py:7-10 creates Measurements/Simulations/Cache/All
    eagerly.  abstractcalculator.py:172,190,221 uses ``datalayer.Cache``
    directly, so a stale binding would silently bypass mongomock.
    """
    import hera.datalayer as dl

    dl.Measurements = dl.Measurements_Collection()
    dl.Simulations = dl.Simulations_Collection()
    dl.Cache = dl.Cache_Collection()
    dl.All = dl.AbstractCollection()


def reset():
    """Drop the in-memory database. Called between tests."""
    from mongoengine.connection import get_connection

    get_connection(UNIT_DB_ALIAS).drop_database(UNIT_DB_NAME)
