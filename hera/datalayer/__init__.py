from .document import getDBObject, createDBConnection, getMongoConfigFromJson
from .collection import AbstractCollection,Measurements_Collection,Simulations_Collection,Analysis_Collection
from .project import getProjectList
from .project import Project

from .utils import dictToMongoQuery
from .datahandler import datatypes

Measurements = Measurements_Collection()
Simulations = Simulations_Collection()
Analysis = Analysis_Collection()

All = AbstractCollection()
