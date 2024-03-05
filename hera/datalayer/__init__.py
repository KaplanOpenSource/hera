from .document import getDBObject, createDBConnection, getMongoConfigFromJson,getDBNamesFromJSON,removeConnection,addOrUpdateDatabase,getMongoJSON
from .collection import AbstractCollection,Measurements_Collection,Simulations_Collection,Cache_Collection
from .project import getProjectList,Project,createProjectDirectory
from .document import nonDBMetadataFrame
from .datahandler import datatypes

Measurements = Measurements_Collection()
Simulations = Simulations_Collection()
Cache = Cache_Collection()
All = AbstractCollection()

