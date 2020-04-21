from mongoengine import *
import os
import json
import getpass
from .metadataDocument import MetadataFrame

dbObjects = {}

def getMongoJSON():
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

    :param mongoConfig: dict
                defines the connection to the DB:

                - dbName : the name of the database.
                - dbIP   : the IP of the database
                - username: the unsername to log in with.
                - password : the user password.

    :param alias: str
            An alternative alias. Used mainly for parallel applications.

    :return:
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
        - Analysis:      documents of the analysis.
        - Simulations:   documents of the simulations.

    :param user:
            The username to register the connection under.
    :param mongoConfig; dict
            defines the connection to the DB.
            see connectToDatabase for details.
    :return:
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

    new_Analysis = type('Analysis', (new_Metadata,), {})
    dbDict['Analysis'] = new_Analysis

    new_Simulations = type('Simulations', (new_Metadata,), {})
    dbDict['Simulations'] = new_Simulations

    dbObjects[user] = dbDict

    return dbDict


# ---------------------default connections--------------------------
for user in getDBNamesFromJSON():
    createDBConnection(user=user,
                       mongoConfig=getMongoConfigFromJson(user=user)
                       )
# -------------------------------------------------------------------


def getDBObject(objectName, user=None):
    """
    Returns the mongoengine object(objectName) of a given user from the global dictionary dbObjects.

    :param objectName:
    :param user:
    :return:
    """
    user = getpass.getuser() if user is None else user
    return dbObjects[user][objectName]

