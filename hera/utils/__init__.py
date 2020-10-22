from unum.units import *
from unum import Unum
from unum import NameConflictError

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



tonumber = lambda x,theunit: x.asNumber(theunit) if isinstance(x,Unum) else x
tounum = lambda x,theunit: x.asUnit(theunit) if isinstance(x,Unum) else x*theunit


def toMeteorologicalAngle(mathematical_angle):
	return (270-mathematical_angle) if ((270-mathematical_angle)   >= 0)   else (630-mathematical_angle)
def toMatematicalAngle(meteorological_angle):
	return (270-meteorological_angle) if ((270-meteorological_angle) >= 0) else (630 - meteorological_angle)

def toAzimuthAngle(mathematical_angle):
	return (90-mathematical_angle) if ((90-mathematical_angle) >= 0) else (450 - mathematical_angle)


try:
	dosage = Unum.unit('dosage',mg*min/m**3,"Exposure dosage")
# except NameConflictError:
except NameConflictError:
	pass
