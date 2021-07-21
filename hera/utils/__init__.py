import logging
from unum import Unum

from .query import andClause,dictToMongoQuery
from .ConvertJSONtoConf import ConvertJSONtoConf



tonumber = lambda x,theunit: x.asNumber(theunit) if isinstance(x,Unum) else x
tounit   = lambda x,theunit: x.asUnit(theunit) if isinstance(x,Unum) else x*theunit

toMeteorologicalAngle = lambda mathematical_angle: (270-mathematical_angle) if ((270-mathematical_angle) >= 0) else (630-mathematical_angle)
toMathematicalAngle  = toMeteorologicalAngle
toAzimuthAngle = lambda mathematical_angle: (90 - mathematical_angle) if ((90 - mathematical_angle) >= 0) else (450 - mathematical_angle)

#############

class loggedObject:

    _logger = None

    @property
    def logger(self):
        return self._logger

    def __init__(self):
        name = ".".join(str(self.__class__)[8:-2].split(".")[1:])
        self._logger = logging.getLogger(name)

