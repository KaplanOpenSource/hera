import logging
from unum import Unum

def andClause(excludeFields=[], **kwargs):

    L = []
    for key, value in kwargs.items():
        if key in excludeFields:
            continue

        if isinstance(value, list):
            conditionStr = "%s in %s"
        elif isinstance(value, str):
            conditionStr = "%s == '%s'"
        elif isinstance(value, dict):
            conditionStr = "%s " + value['operator'] + " %s"
            value = value['value']
        else:
            conditionStr = "%s == %s"

        L.append(conditionStr % (key, value))

    return " and ".join(L)


toMeteorologicalAngle = lambda mathematical_angle: (270-mathematical_angle) if ((270-mathematical_angle) >= 0) else (630-mathematical_angle)
toMathematicalAngle  = toMeteorologicalAngle
toAzimuthAngle = lambda mathematical_angle: (90 - mathematical_angle) if ((90 - mathematical_angle) >= 0) else (450 - mathematical_angle)


from .query import andClause,dictToMongoQuery
from .unum import tounit,tonumber,tounum


#############

class loggedObject:

    _logger = None

    @property
    def logger(self):
        return self._logger

    def __init__(self):
        name = ".".join(str(self.__class__)[8:-2].split(".")[1:])
        self._logger = logging.getLogger(name)

