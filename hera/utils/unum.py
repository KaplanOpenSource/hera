from unum import Unum

tonumber = lambda x,theunit: x.asNumber(theunit) if isinstance(x,Unum) else x
tounit   = lambda x,theunit: x.asUnit(theunit) if isinstance(x,Unum) else x*theunit
tounum   = tounit
