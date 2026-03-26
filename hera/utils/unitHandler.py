"""
Unit handling module for Hera.

**Pint is the preferred unit system.** Use ``ureg`` for all new code.
Unum is supported as a fallback for backward compatibility with older code.

Usage::

    from hera.utils.unitHandler import ureg, Quantity

    speed = 5.0 * ureg.m / ureg.s
    area = 10 * ureg.dunam          # custom: 1 dunam = 1000 m²
"""

from hera.utils.logging import get_logger
from pint import Unit, UnitRegistry
from pint import Quantity
from pint.errors import UndefinedUnitError, DimensionalityError
from deprecated import deprecated

# ---------------------------------------------------------------------------
# Pint unit registry (preferred)
# ---------------------------------------------------------------------------

ureg = UnitRegistry()
ureg.define('mmH2O = atm / 10197.162129779')
ureg.define('dunam = 1000*m**2')

# Pint-based unit constants for convenience
celsius = ureg.degC
K = ureg.kelvin

# ---------------------------------------------------------------------------
# Unum support (backward compatibility only)
# ---------------------------------------------------------------------------

try:
    from unum import Unum
    from unum import NameConflictError
    import unum.units as _unum_units_module
    import re

    # Export unum unit symbols (m, kg, s, etc.) for backward compatibility.
    # NOTE: celsius and K are NOT overridden — pint versions are kept.
    _unum_units = {k: v for k, v in vars(_unum_units_module).items() if not k.startswith('_')}
    _unum_compat = {k: v for k, v in _unum_units.items() if k not in ('celsius', 'K')}
    globals().update(_unum_compat)

    # Unum celsius/K available as _unum_celsius/_unum_K for
    # backward-compatible .asNumber() calls in old code paths.
    _unum_celsius = _unum_units.get('celsius', None)
    _unum_K = _unum_units.get('K', None)

    # Override celsius/K with unum versions since existing code like
    # temperature.asNumber(celsius) needs unum units.
    # New code should use ureg.degC and ureg.kelvin instead.
    if _unum_celsius is not None:
        celsius = _unum_celsius
    if _unum_K is not None:
        K = _unum_K

    unumSupport = True
except ImportError:
    class Unum:
        """Placeholder when unum is not installed."""
        pass
    _unum_units = {}
    _unum_celsius = None
    _unum_K = None
    unumSupport = False


# ---------------------------------------------------------------------------
# Conversion helpers
# ---------------------------------------------------------------------------

def tonumber(x, theunit):
    """
    Convert a value with units to a plain number in the given unit.

    Parameters
    ----------
    x : Unum, Quantity, or numeric
        The value to convert.
    theunit : unit
        The target unit.

    Returns
    -------
    float
    """
    if isinstance(x, Unum) and unumSupport:
        ret = x.asNumber(theunit)
    elif isinstance(x, Quantity):
        ret = x.to(theunit).magnitude
    else:
        ret = x
    return ret


def tounit(x, theunit):
    """
    Convert a value to the given unit, handling both pint and unum.

    Parameters
    ----------
    x : Unum, Quantity, or numeric
        The value to convert.
    theunit : Unit, Unum, or str
        The target unit.

    Returns
    -------
    Quantity or Unum
    """
    logger = get_logger(None, "hera.utils.tounit")
    if isinstance(x, Unum) and unumSupport:
        logger.warning("Please prefer using pint (ureg) for units")
        if isinstance(theunit, Unit):
            ret = unumToPint(x).to(theunit)
        elif isinstance(theunit, Unum):
            ret = x.asUnit(theunit)
        else:
            logger.error(f"can't convert {x} to units {theunit}")
            ret = x
    elif isinstance(x, Quantity):
        if isinstance(theunit, Unum) and unumSupport:
            logger.warning("Please prefer using pint (ureg) for units")
            ret = pintToUnum(x).asUnit(theunit)
        elif isinstance(theunit, (Unit, str)):
            ret = x.to(theunit)
        else:
            logger.error(f"can't convert {x} to units {theunit}")
            ret = x
    else:
        ret = Quantity(x, theunit)
    return ret


tounum = tounit


# ---------------------------------------------------------------------------
# Unum custom units and conversion helpers (only when unum is installed)
# ---------------------------------------------------------------------------

if unumSupport:
    # Access unum unit objects from the private dict
    _m = _unum_units['m']
    _s = _unum_units['s']
    _g = _unum_units['g']
    _kg = _unum_units['kg']
    _cm = _unum_units['cm']
    _N = _unum_units['N']
    _bar = _unum_units['bar']
    _L = _unum_units['L']

    # Define custom unum units
    try:
        _unum_atm = Unum.unit('atm', 1.01325 * _bar, 'atmosphere')
    except NameConflictError:
        _unum_atm = Unum.unit('atm', conv=None, name='atmosphere')

    try:
        _unum_mbar = Unum.unit('mbar', _bar / 1000, 'millibar')
    except NameConflictError:
        _unum_mbar = Unum.unit('mbar', conv=None, name='millibar')

    try:
        _unum_mmHg = Unum.unit('mmHg', _unum_atm / 760., 'mmHg = 1 torr')
    except NameConflictError:
        _unum_mmHg = Unum.unit('mmHg', conv=None, name='mmHg = 1 torr')

    try:
        _unum_mmH2O = Unum.unit('mmH2O', _unum_atm / 10197.162129779, 'mmH2O')
    except NameConflictError:
        _unum_mmH2O = Unum.unit('mmH2O', conv=None, name='mmH2O')

    try:
        _unum_torr = Unum.unit('torr', _unum_atm / 760., 'torr = 1 mmHg')
    except NameConflictError:
        _unum_torr = Unum.unit('torr', conv=None, name='torr = 1 mmHg')

    try:
        _unum_dyne = Unum.unit('dyne', 1e-5 * _N, 'dyne')
    except NameConflictError:
        _unum_dyne = Unum.unit('dyne', conv=None, name='dyne')

    try:
        _unum_poise = Unum.unit('poise', _g / _cm / _s, 'poise')
    except NameConflictError:
        _unum_poise = Unum.unit('poise', conv=None, name='poise')

    try:
        _unum_cpoise = Unum.unit('cpoise', _unum_poise / 10., 'centipoise')
    except NameConflictError:
        _unum_cpoise = Unum.unit('cpoise', conv=None, name='centipoise')

    try:
        Unum.unit('mL', 1e-3 * _L)
    except NameConflictError:
        pass

    try:
        Unum.unit('uL', 1e-6 * _L)
    except NameConflictError:
        pass

    try:
        Unum.unit('nL', 1e-9 * _L)
    except NameConflictError:
        pass

    # Map pint unit names to unum unit objects (for conversion)
    PINT_TO_UNUM_MAP = {
        'meter': _m,
        'metre': _m,
        "hour": _unum_units.get('h', _unum_units.get('hour', None)),
        "centimeter": _cm,
        "mmHg": _unum_mmHg,
        "mbar": _unum_mbar,
        'second': _s,
        "milligram": _unum_units.get('mg', None),
        "minute": _unum_units.get('min', None),
        "atm": _unum_atm,
        "milliliter": _unum_units.get('mL', None),
        'kilogram': _kg,
        'ampere': _unum_units.get('A', None),
        'kelvin': _unum_units.get('K', None),
        'mole': _unum_units.get('mol', None),
        'candela': _unum_units.get('cd', None),
        'radian': _unum_units.get('rad', None),
        'steradian': _unum_units.get('sr', None),
        'hertz': _unum_units.get('Hz', None),
        'newton': _N,
        'pascal': _unum_units.get('Pa', None),
        'joule': _unum_units.get('J', None),
        'watt': _unum_units.get('W', None),
        'coulomb': _unum_units.get('C', None),
        'volt': _unum_units.get('V', None),
        'farad': _unum_units.get('F', None),
        'ohm': _unum_units.get('ohm', None),
        'siemens': _unum_units.get('S', None),
        'weber': _unum_units.get('Wb', None),
        'tesla': _unum_units.get('T', None),
        'henry': _unum_units.get('H', None),
        'lumen': _unum_units.get('lm', None),
        'lux': _unum_units.get('lx', None),
        'becquerel': _unum_units.get('Bq', None),
        'gray': _unum_units.get('Gy', None),
        'sievert': _unum_units.get('Sv', None),
        'katal': _unum_units.get('kat', None),
    }
    # Remove None entries
    PINT_TO_UNUM_MAP = {k: v for k, v in PINT_TO_UNUM_MAP.items() if v is not None}

    @deprecated(reason="Doesn't work for some cases. Move to Pint")
    def convert_unum_units_to_eval_str(unit_str):
        """Convert a unum-style unit string to an eval-safe expression."""
        def replace_exponents(match):
            """Replace unit-exponent pairs with Python power syntax."""
            unit = match.group(1)
            exponent = match.group(2)
            return f"{unit}**{exponent}"
        pattern = re.compile(r'([a-zA-Z]+)(-?\d+)')
        unit_str = re.sub(pattern, replace_exponents, unit_str)
        return unit_str

    @deprecated(reason="Doesn't work for some cases")
    def unumToStr(obj):
        """Convert a Unum object to a string representation."""
        if isinstance(obj, Unum):
            objStr = convert_unum_units_to_eval_str(str(obj))
            ret = objStr.replace(" [", "*").replace("]", "")
        else:
            ret = str(obj)
        return ret

    @deprecated(reason="Doesn't work for some cases")
    def strToUnum(value):
        """Convert a string to a Unum object."""
        if isinstance(value, Unum):
            ret = value
        else:
            try:
                ret = eval(str(value))
            except:
                ret = value
        return ret

    def extractUnumUnitsFromPint(pint_quantity):
        """Extract unum unit equivalent from a pint Quantity."""
        units = pint_quantity._units
        unum_unit = 1 * _m / _m  # unitless
        for unit_name, power in units.items():
            if unit_name not in PINT_TO_UNUM_MAP:
                raise ValueError(f"Unit '{unit_name}' not mapped to Unum.")
            unum_unit *= PINT_TO_UNUM_MAP[unit_name] ** power
        return unum_unit

    def pintToUnum(pint_quantity):
        """
        Convert a pint Quantity to a Unum object.

        Parameters
        ----------
        pint_quantity : Quantity
            The pint value to convert.

        Returns
        -------
        Unum
        """
        magnitude = pint_quantity.magnitude
        units = pint_quantity._units
        unum_unit = 1 * _m / _m
        for unit_name, power in units.items():
            if unit_name not in PINT_TO_UNUM_MAP:
                raise ValueError(f"Unit '{unit_name}' not mapped to Unum.")
            unum_unit *= PINT_TO_UNUM_MAP[unit_name] ** power
        return magnitude * unum_unit

    def unumToPint(unum_obj, value=1.0):
        """
        Convert a Unum object to a pint Quantity.

        Parameters
        ----------
        unum_obj : Unum
            The unum value to convert.
        value : float
            Optional numerical value.

        Returns
        -------
        Quantity
        """
        if isinstance(unum_obj, Quantity):
            return unum_obj
        unit_str = str(unum_obj)
        pint_str = convert_unum_units_to_eval_str(unit_str)
        return value * ureg.parse_expression(pint_str)

    def unumToBaseUnits(unum_obj):
        """Convert a Unum object to MKS base units."""
        pint_obj = unumToPint(unum_obj)
        standardize = pint_obj.to_base_units()
        return pintToUnum(standardize)
