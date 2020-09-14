from .document import getDBObject, createDBConnection, getMongoConfigFromJson,getDBNamesFromJSON
from .collection import AbstractCollection,Measurements_Collection,Simulations_Collection,Cache_Collection
from .project import getProjectList
from .project import Project, ProjectMultiDB, ProjectMultiDBPublic

from .utils import dictToMongoQuery
from .datahandler import datatypes


Measurements = Measurements_Collection()
Simulations = Simulations_Collection()
Cache = Cache_Collection()
All = AbstractCollection()

