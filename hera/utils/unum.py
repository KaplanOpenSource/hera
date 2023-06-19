from unum import Unum
from unum.units import *
from unum import NameConflictError

tonumber = lambda x,theunit: x.asNumber(theunit) if isinstance(x,Unum) else x
tounit   = lambda x,theunit: x.asUnit(theunit) if isinstance(x,Unum) else x*theunit
tounum   = tounit

try:
    atm   = Unum.unit('atm',1.01325*bar,'atmosphere')
except NameConflictError as e:
    atm = Unum.unit('atm',conv=None,name='atmosphere')

try:
    mbar  = Unum.unit('mbar',bar/1000,'millibar')
except NameConflictError as e:
    mbar = Unum.unit('mbar',conv=None,name='millibar')

try:
    mmHg  = Unum.unit('mmHg',atm/760.,'mmHg = 1 torr')
except NameConflictError as e:
    mmHg = Unum.unit('mmHg',conv=None,name='mmHg = 1 torr')

try:
    torr  = Unum.unit('torr',atm/760.,'torr = 1 mmHg')
except NameConflictError as e:
    torr = Unum.unit('torr',conv=None,name='torr = 1 mmHg')

try:
    dyne  = Unum.unit('dyne',1e-5*N,'dyne')
except NameConflictError as e:
    dyne = Unum.unit('dyne',conv=None,name='dyne')

try:
    poise = Unum.unit('poise',g/cm/s,'poise')
except NameConflictError as e:
    poise = Unum.unit('poise',conv=None,name='poise')

try:
    cpoise = Unum.unit('cpoise',poise/10.,'centipoise')
except NameConflictError as e:
    cpoise = Unum.unit('cpoise',conv=None,name='centipoise')

try:
    mL = Unum.unit('mL', 1e-3 * L)
except NameConflictError:
    pass

try:
    uL = Unum.unit('uL', 1e-6 * L)
except NameConflictError:
    pass

try:
    nL = Unum.unit('nL', 1e-9 * L)
except NameConflictError:
    pass



def unumToStr(obj):
    if isinstance(obj, Unum):
        ret = str(obj).replace(" [", "*").replace("]", "")
    else:
        ret = str(obj)
    return ret

def strToUnum(value):
    try:
        ret = eval(str(value))
    except:
        ret = value
    return ret
