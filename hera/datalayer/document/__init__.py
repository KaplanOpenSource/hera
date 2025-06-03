from mongoengine import *
from mongoengine.connection import disconnect
import os
import json
import getpass
import logging
from hera.datalayer.document.metadataDocument import MetadataFrame,nonDBMetadataFrame

dbObjects = {}

def getMongoJSON():
    """
    Return the definition of the connection.

    Returns
    -------
        dict
    """
    configFile = os.path.join(os.environ.get('HOME'), '.pyhera', 'config.json')
    if os.path.isfile(configFile):
        with open(configFile, 'r') as jsonFile:
            mongoConfig = json.load(jsonFile)
    else:
        configData = {getpass.getuser(): dict(username='{username}',
                                              password='{password}',
                                              dbIP='{databaseIP}',
                                              dbName='{databaseName}'
                                              )
                      }

        if not os.path.isdir(os.path.dirname(configFile)):
            os.makedirs(os.path.dirname(configFile))

        with open(configFile, 'w') as jsonFile:
            json.dump(configData, jsonFile, indent=4, sort_keys=True)

        errorMessage = "The config file doesn't exist in the default directory.\n" \
                       "A default config data file named '{}' was created. Please fill it and try again.".format(
            configFile)

        raise IOError(errorMessage)
    return mongoConfig

def addOrUpdateDatabase(connectionName, username, password, databaseIP, databaseName):
    """
        Adding or updating the configuration file.

    Parameters
    ----------
    connectionName : str
            The name of the connection
    username : str
            The user name in the database
    password : str
            The password for the database
    databaseIP : str
            The IP of the database server.
    databaseName : str
            The name of the database

    Returns
    -------

    """
    logger = logging.getLogger("hera.datalayer.addOrUpdateDatabase")
    logger.info(f"Creating the connection {connectionName}")
    configFile = os.path.join(os.environ.get('HOME'), '.pyhera', 'config.json')
    if os.path.isfile(configFile):
        with open(configFile, 'r') as jsonFile:
            mongoConfig = json.load(jsonFile)
    else:
        raise FileNotFoundError(f"{configFile} does not exis. Create it first by importing hera and filling up connection names")

    mongoConfig[connectionName] =  dict(username=username, password=password, dbIP=databaseIP, dbName=databaseName)
    logger.debug(f"Creating the connection with the data: {mongoConfig[connectionName]}")

    with open(configFile, 'w') as jsonFile:
        json.dump(mongoConfig, jsonFile, indent=4, sort_keys=True)

    logger.debug(f"Updating the config file {configFile}")

def removeConnection(connectionName):
    """
        Removing the connection name from the configuration file.

    Parameters
    ----------
    connectionName

    Returns
    -------

    """
    logger = logging.getLogger("hera.datalayer.addOrUpdateDatabase")
    logger.info(f"Removing the connection {connectionName}")
    configFile = os.path.join(os.environ.get('HOME'), '.pyhera', 'config.json')
    if os.path.isfile(configFile):
        with open(configFile, 'r') as jsonFile:
            mongoConfig = json.load(jsonFile)
    else:
        raise FileNotFoundError(f"{configFile} does not exis. Create it first by importing hera and filling up connection names")

    deletedConnection = mongoConfig[connectionName]
    del mongoConfig[connectionName]
    logger.info(f"Removing the connection {connectionName}: {deletedConnection}")

    with open(configFile, 'w') as jsonFile:
        json.dump(mongoConfig, jsonFile, indent=4, sort_keys=True)

    print(f"Removing the connection {connectionName}")
    print(json.dumps(deletedConnection,indent=4))
    logger.debug(f"Saved the config file {configFile}")


def getDBNamesFromJSON():
    mongoConfigJSON = getMongoJSON()
    return [x for x in mongoConfigJSON.keys()] 

def getMongoConfigFromJson(connectionName=None):
    """
    Get the mongoConfig of a user from .pyhera/config.json

    :param connectionName:
    :return:
    """
    mongoConfigJSON = getMongoJSON()
    mongoConfig = mongoConfigJSON[getpass.getuser()] if connectionName is None else mongoConfigJSON[connectionName]
    
    return mongoConfig
    ## build the connection to the db.


def connectToDatabase(mongoConfig,alias=None):
    """
    Creates a connection to the database according to the mongoConfig.

    Parameters
    ----------

    mongoConfig: dict
                defines the connection to the DB:

                - dbName : the name of the database.
                - dbIP   : the IP of the database
                - username: the unsername to log in with.
                - password : the user password.

                str
                defines a connection string.

                username:password@dbIP/dbName

    alias: str
            An alternative alias. Used mainly for parallel applications.

    Returns
    -------
        mongodb connection.
    """
    if isinstance(mongoConfig,str):
        mongoConfig = parseConnectionString(mongoConfig)

    alias = '%s-alias' % mongoConfig['dbName'] if alias is None else alias

    disconnect(alias)

    con = connect(alias=alias,
            host=mongoConfig['dbIP'],
            db=mongoConfig['dbName'],
            username=mongoConfig['username'],
            password=mongoConfig['password'],
            authentication_source='admin'
            )
    return con


def createDBConnection(connectionName, mongoConfig, alias=None):
    """
    Creates a connection to the database.
    Creates mongoengine objects and saving them to a global dictionary dbObjects.

    saves DB objects for the user in a DBdict:

        - connection: the connection to the db. has the alias [dbname]-alias
                      or the given alias name.

        - Metadata: the meta data object that holds all the documets.
        - Measurements:  documents of the measurements.old.
        - Cache:      documents of the cache.
        - Simulations:   documents of the simulations.old.

    Parameters
    ----------

    connectionName: str
            The username to register the connection under.
    mongoConfig; dict
            defines the connection to the DB.
            see connectToDatabase for details.
    alias: str
            The name of the connection.
            Used to prevent two connections with the same name.
            if None, use the user name.

    Returns
    -------
        dict.
        return the DBdict.
    """
    dbDict = {}
    if isinstance(mongoConfig,str):
        mongoConfig = parseConnectionString(mongoConfig)


    con = connectToDatabase(mongoConfig=mongoConfig,alias=alias)

    dbName = mongoConfig['dbName']

    new_Metadata = type('Metadata', (DynamicDocument, MetadataFrame), {'meta': {'db_alias': '%s-alias' % dbName,
                                                                                'allow_inheritance': True,
                                                                                'auto_create_indexes': True,
                                                                                'indexes': ['projectName']
                                                                                },
                                                                       '__str__' : lambda self: MetadataFrame.__str__(self)
                                                                       }
                        )

    dbDict['connection'] = con
    dbDict['Metadata'] = new_Metadata

    new_Measurements = type('Measurements', (new_Metadata,), {})
    dbDict['Measurements'] = new_Measurements

    new_Cache = type('Cache', (new_Metadata,), {})
    dbDict['Cache'] = new_Cache

    new_Simulations = type('Simulations', (new_Metadata,), {})
    dbDict['Simulations'] = new_Simulations

    dbObjects[connectionName] = dbDict

    return dbDict



def getDBObject(objectName, connectionName=None):
    """
    Returns the mongoengine object(objectName) of a given user from the global dictionary dbObjects.

    Parameters
    ----------

    objectName: str
        The name of the object to return
    user: str
        Connection name
    :return:
        mongoengine object
    """
    connectionName = getpass.getuser() if connectionName is None else connectionName

    try:
        dbs = dbObjects[connectionName]
    except KeyError:
        allusers = ",".join([x for x in dbObjects.keys()])
        raise KeyError(f"user {connectionName} not found. Must be one of: {allusers}")


    try:
        ret = dbs[objectName]
    except KeyError:
        allobjs = ",".join([x for x in dbs.keys()])
        raise KeyError(f"object {objectName} not found. Must be one of: {allobjs}")

    return ret

def parseConnectionString(conStr):
    username, password = conStr.split("@")[0].split(":")
    dbIP, dbName = conStr.split("@")[1].split("/")
    return dict(username=username,
                       password=password,
                       dbName=dbName,
                       dbIP=dbIP)


# ---------------------default connections--------------------------
for user in getDBNamesFromJSON():
    createDBConnection(connectionName=user,
                       mongoConfig=getMongoConfigFromJson(connectionName=user)
                       )
# -------------------------------------------------------------------

