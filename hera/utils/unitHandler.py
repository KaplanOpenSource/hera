from unum import Unum
from unum.units import *
from unum import NameConflictError
import re

from pint import UnitRegistry
from pint import Quantity
from pint.errors import UndefinedUnitError, DimensionalityError

from deprecated import deprecated

def tonumber(x,theunit):
    if isinstance(x,Unum):
        ret = x.asNumber(theunit)
    elif isinstance(x,Quantity):
        ret = x.to(theunit).magnitude
    else:
        ret = x

    return ret

def tounit(x,theunit):
    if isinstance(x,Unum):
        ret = x.asUnit(theunit)
    elif isinstance(x,Quantity):
        ret = x.to(theunit)
    else:
        ret = Quantity(x,theunit)
    return ret

# tonumber = lambda x, theunit: x.asNumber(theunit) if isinstance(x, Unum) else x
# tounit = lambda x, theunit: x.asUnit(theunit) if isinstance(x, Unum) else x * theunit
tounum = tounit

# Initialize the Pint unit registry
ureg = UnitRegistry()
ureg.define('mmH2O = atm / 10197.162129779')
ureg.define('dunam = 1000*m**2')

try:
    atm = Unum.unit('atm', 1.01325 * bar, 'atmosphere')
except NameConflictError as e:
    atm = Unum.unit('atm', conv=None, name='atmosphere')

try:
    mbar = Unum.unit('mbar', bar / 1000, 'millibar')
except NameConflictError as e:
    mbar = Unum.unit('mbar', conv=None, name='millibar')

try:
    mmHg = Unum.unit('mmHg', atm / 760., 'mmHg = 1 torr')
except NameConflictError as e:
    mmHg = Unum.unit('mmHg', conv=None, name='mmHg = 1 torr')

try:
    mmH2O = Unum.unit('mmH2O', atm / 10197.162129779, 'mmH2O = 0.0000980665bar')
except NameConflictError as e:
    mmH2O = Unum.unit('mmH2O', conv=None, name='mmH2O = 1 torr')

try:
    torr = Unum.unit('torr', atm / 760., 'torr = 1 mmHg')
except NameConflictError as e:
    torr = Unum.unit('torr', conv=None, name='torr = 1 mmHg')

try:
    dyne = Unum.unit('dyne', 1e-5 * N, 'dyne')
except NameConflictError as e:
    dyne = Unum.unit('dyne', conv=None, name='dyne')

try:
    poise = Unum.unit('poise', g / cm / s, 'poise')
except NameConflictError as e:
    poise = Unum.unit('poise', conv=None, name='poise')

try:
    cpoise = Unum.unit('cpoise', poise / 10., 'centipoise')
except NameConflictError as e:
    cpoise = Unum.unit('cpoise', conv=None, name='centipoise')

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


@deprecated(reason="Doesn't work for some cases For example kg/s/m. Move to Pint")
def convert_unum_units_to_eval_str(unit_str):
    """
    Converts a string like 'kg/m3' to 'kg/m**3' for eval-safe unit expression.
    It assumes the input string is a valid unum-style unit string.

    Args:
        unit_str (str): A unit string in unum format (e.g., 'kg/m3').

    Returns:
        str: A string formatted for eval (e.g., 'kg/m**3').
    """

    # Convert shorthand exponents like m3 to m**3 (only when preceded by a letter)
    def replace_exponents(match):
        unit = match.group(1)
        exponent = match.group(2)
        return f"{unit}**{exponent}"

    # Regex matches: a letter or unit symbol followed by an integer exponent
    pattern = re.compile(r'([a-zA-Z]+)(-?\d+)')
    # Apply the transformation only to denominator (after '/') or any place it appears
    unit_str = re.sub(pattern, replace_exponents, unit_str)

    return unit_str


@deprecated(reason="Doesn't work for some cases For example kg/s/m")
def unumToStr(obj):
    if isinstance(obj, Unum):
        objStr = convert_unum_units_to_eval_str(str(obj))
        ret = objStr.replace(" [", "*").replace("]", "")
    else:
        ret = str(obj)
    return ret


@deprecated(reason="Doesn't work for some cases For example kg/s/m")
def strToUnum(value):
    if isinstance(value, Unum):
        ret = value
    else:
        try:
            ret = eval(str(value))
        except:
            ret = value
    return ret


# Map pint unit names to unum.unit unit objects
PINT_TO_UNUM_MAP = {
    'meter': m,
    'metre': m,
    "hour" : h,
    "centimeter" : cm,
    "mmHg": mmHg,
    "mbar": mbar,
    'second': s,
    "milligram": mg,
    "minute": min,
    "atm": atm,
    "milliliter": mL,
    'kilogram': kg,
    'ampere': A,
    'kelvin': K,
    'mole': mol,
    'candela': cd,
    'radian': rad,
    'steradian': sr,
    'hertz': Hz,
    'newton': N,
    'pascal': Pa,
    'joule': J,
    'watt': W,
    'coulomb': C,
    'volt': V,
    'farad': F,
    'ohm': ohm,
    'siemens': S,
    'weber': Wb,
    'tesla': T,
    'henry': H,
    'lumen': lm,
    'lux': lx,
    'becquerel': Bq,
    'gray': Gy,
    'sievert': Sv,
    'katal': kat,
}


def extractUnumUnitsFromPint(pint_quantity):
    """
    Converts the pint unit string 'meters**1*second**-2' to unum str m**1*s**-2
    Parameters
    ----------
    pint_units_str : str
        The str of the unum.

    Returns
    -------

    """
    units = pint_quantity._units

    unum_unit = 1 * m / m  # unitless number

    for unit_name, power in units.items():
        if unit_name not in PINT_TO_UNUM_MAP:
            raise ValueError(f"Unit '{unit_name}' not mapped to Unum unit object.")
        unum_base = PINT_TO_UNUM_MAP[unit_name]
        unum_unit *= unum_base ** power

    return unum_unit


def pintToUnum(pint_quantity):
    """
        Conver pint quantity ot Unum.
    Parameters
    ----------
    pint_quantity (pint) : The value to convert.

    Returns
    -------
        unum.
    """
    magnitude = pint_quantity.magnitude
    units = pint_quantity._units

    unum_unit = 1 * m / m  # unitless number

    for unit_name, power in units.items():
        if unit_name not in PINT_TO_UNUM_MAP:
            raise ValueError(f"Unit '{unit_name}' not mapped to Unum unit object.")
        unum_base = PINT_TO_UNUM_MAP[unit_name]
        unum_unit *= unum_base ** power

    return magnitude * unum_unit


def unumToPint(unum_obj, value=1.0):
    """
    Convert a Unum object to a Pint Quantity.

    Args:
        unum_obj (Unum): A unum object (e.g., kg / m**3).
        value (float): Optional numerical value.

    Returns:
        pint.Quantity: A pint quantity with the same unit.
    """
    # Get the string representation of the unum object
    unit_str = str(unum_obj)
    pint_str = convert_unum_units_to_eval_str(unit_str)
    return value * ureg.parse_expression(pint_str)


def unumToBaseUnits(unum_obj):
    """
        Convert the unum object to MKS units.
    Parameters
    ----------
    unum_obj (Unum): A unum object (e.g., kg / m**3).

    Returns
    -------
    Unum
    """
    pint_obj = unumToPint(unum_obj)
    standardize = pint_obj.to_base_units()
    return pintToUnum(standardize)
