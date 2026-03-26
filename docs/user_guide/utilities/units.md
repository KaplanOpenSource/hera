# Unit Handling

Hera uses [pint](https://pint.readthedocs.io/) as its preferred unit system. A shared `UnitRegistry` called `ureg` is available throughout the library. It comes pre-configured with custom units relevant to environmental and hydrological work, such as `dunam` (1 dunam = 1000 m^2) and `mmH2O`.

Legacy code may still use unum, but all new code should use pint via `ureg`.

## Importing

```python
from hera.utils.unitHandler import ureg, Quantity
# Or from the top-level namespace:
from hera.utils import ureg
```

## Creating Quantities

Attach units to numbers by multiplying with `ureg.<unit>`:

```python
speed = 5.0 * ureg.m / ureg.s
area = 10 * ureg.dunam
pressure = 50 * ureg.mmH2O

print(speed)       # 5.0 meter / second
print(area)        # 10 dunam
print(pressure)    # 50 mmH2O
```

## Converting Between Units

```python
# Convert to different units
speed_kmh = speed.to(ureg.km / ureg.hour)
print(speed_kmh)  # 18.0 kilometer / hour

# Convert to base (MKS) units
area_m2 = area.to_base_units()
print(area_m2)  # 10000.0 meter ** 2
```

## Extracting the Numeric Value

Use `tonumber` to strip units and get a plain float in a specific unit:

```python
from hera.utils.unitHandler import tonumber

value = tonumber(speed, ureg.m / ureg.s)
print(value)  # 5.0
print(type(value))  # <class 'float'>
```

This works transparently with both pint `Quantity` objects and legacy unum objects, as well as plain numbers (which pass through unchanged).

## Converting with `tounit`

Use `tounit` to ensure a value carries the desired unit, regardless of what it started as:

```python
from hera.utils.unitHandler import tounit

# Plain number gets units attached
q = tounit(100, ureg.m)
print(q)  # 100 meter

# Existing quantity gets converted
q2 = tounit(speed, ureg.km / ureg.hour)
print(q2)  # 18.0 kilometer / hour
```

## Using Units with Toolkits

When storing configurations that contain units (for example, in a MongoDB-backed toolkit), combine `ureg` quantities with the JSON utilities:

```python
from hera.utils import ureg, ConfigurationToJSON, JSONToConfiguration

config = {
    "windSpeed": 3.5 * ureg.m / ureg.s,
    "domainSize": 2 * ureg.dunam,
}

# Serialize to JSON-safe dict (with split magnitude/units)
json_safe = ConfigurationToJSON(config, splitUnits=True)
# {'windSpeed': {'magnitude': 3.5, 'units': '3.5 meter / second', 'value': '3.5 meter / second'},
#  'domainSize': {'magnitude': 2, 'units': '2 dunam', 'value': '2 dunam'}}

# Restore back to pint Quantity objects
restored = JSONToConfiguration(json_safe)
print(restored["windSpeed"])  # 3.5 meter / second
```

For the full details on saving and loading configurations with units — including `splitUnits` for unit-aware database queries, `JSONVariations` for parameter sweeps, and `compareJSONS` for diffing configurations — see **[JSON Utilities](json.md)**.

## Custom Units Reference

| Unit | Definition | Typical Use |
|------|-----------|-------------|
| `dunam` | 1000 m^2 | Land area (common in Israel) |
| `mmH2O` | atm / 10197.162129779 | Water column pressure |
