from mongoengine import *
import os
import json
import getpass
from .metadataDocument import MetadataFrame

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

def getDBNamesFromJSON():
    mongoConfigJSON = getMongoJSON()
    return [x for x in mongoConfigJSON.keys()] 

def getMongoConfigFromJson(user=None):
    """
    Get the mongoConfig of a user from .pyhera/config.json

    :param user:
    :return:
    """
    mongoConfigJSON = getMongoJSON()
    mongoConfig = mongoConfigJSON[getpass.getuser()] if user is None else mongoConfigJSON[user]
    
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

    alias: str
            An alternative alias. Used mainly for parallel applications.

    Returns
    -------
        mongodb connection.
    """
    alias = '%s-alias' % mongoConfig['dbName'] if alias is None else alias

    con = connect(alias=alias,
            host=mongoConfig['dbIP'],
            db=mongoConfig['dbName'],
            username=mongoConfig['username'],
            password=mongoConfig['password'],
            authentication_source='admin'
            )
    return con


def createDBConnection(user, mongoConfig,alias=None):
    """
    Creates a connection to the database.
    Creates mongoengine objects and saving them to a global dictionary dbObjects.

    saves DB objects for the user in a DBdict:

        - connection: the connection to the db. has the alias [dbname]-alias
                      or the given alias name.

        - Metadata: the meta data object that holds all the documets.
        - Measurements:  documents of the measurements.
        - Cache:      documents of the cache.
        - Simulations:   documents of the simulations.

    Parameters
    ----------

    user: str
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

    con = connectToDatabase(mongoConfig=mongoConfig,alias=alias)

    dbName = mongoConfig['dbName']

    new_Metadata = type('Metadata', (DynamicDocument, MetadataFrame), {'meta': {'db_alias': '%s-alias' % dbName,
                                                                                'allow_inheritance': True,
                                                                                'auto_create_indexes': True,
                                                                                'indexes': [('geometery', '2dsphere')]
                                                                                }
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

    dbObjects[user] = dbDict

    return dbDict



def getDBObject(objectName, user=None):
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
    user = getpass.getuser() if user is None else user

    try:
        dbs = dbObjects[user]
    except KeyError:
        allusers = ",".join([x for x in dbObjects.key()])
        raise KeyError(f"user {user} not found. Must be one of: {allusers}")


    try:
        ret = dbs[objectName]
    except KeyError:
        allobjs = ",".join([x for x in dbs.keys()])
        raise KeyError(f"object {objectName} not found. Must be one of: {allobjs}")

    return ret



# ---------------------default connections--------------------------
for user in getDBNamesFromJSON():
    createDBConnection(user=user,
                       mongoConfig=getMongoConfigFromJson(user=user)
                       )
# -------------------------------------------------------------------

