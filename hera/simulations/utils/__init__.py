from unum import Unum
from unum.units import *
from unum import NameConflictError

def toUnum(x,unit):
    return x.asUnit(unit) if isinstance(x,Unum) else x*unit

def toNumber(x,unit):
    return x.asNumber(unit) if isinstance(x,Unum) else x


########################## Units for the model.
try:
	atm   = Unum.unit('atm',1.01325*bar,'atmosphere')
	mbar  = Unum.unit('mbar',bar/1000,'millibar')

	mmHg  = Unum.unit('mmHg',atm/760.,'mmHg = 1 torr')
	torr  = Unum.unit('torr',atm/760.,'torr = 1 mmHg')

	dyne  = Unum.unit('dyne',1e-5*N,'dyne')
	poise = Unum.unit('poise',g/cm/s,'poise')
	cpoise = Unum.unit('cpoise',poise/10.,'centipoise')
except NameConflictError:
	pass

####

from .matplotlibCountour import toGeopandas
