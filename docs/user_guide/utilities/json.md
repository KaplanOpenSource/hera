# JSON Utilities

The `hera.utils.jsonutils` module provides functions for serializing configurations that contain physical units, comparing JSON structures side by side, generating parameter sweeps, and manipulating values by JSON path.

## Importing

```python
from hera.utils import (
    ConfigurationToJSON,
    JSONToConfiguration,
    JSONVariations,
    compareJSONS,
    setJSONPath,
)
```

## ConfigurationToJSON -- Serialize Configs with Units

Convert a Python dict (which may contain pint `Quantity` values) into a JSON-safe dict where every value is a string or plain number.

```python
from hera.utils import ureg, ConfigurationToJSON

config = {
    "solver": {
        "velocity": 10.0 * ureg.m / ureg.s,
        "iterations": 500,
    }
}

# Basic serialization -- units become strings
json_dict = ConfigurationToJSON(config)
# {'solver': {'velocity': '10.0 meter / second', 'iterations': 500}}
```

### Splitting Units for Queryable Storage

When `splitUnits=True`, each quantity is expanded into a sub-dict with `magnitude`, `units`, and `value` keys (`magnitude` and `units` describe the quantity after optional standardization; `value` is the original quantity string and is omitted when `keepOriginalUnits=False`). This is useful when you need to query by numeric magnitude (e.g., in MongoDB with `$lt` / `$gt`):

```python
json_dict = ConfigurationToJSON(config, splitUnits=True)
# {'solver': {
#     'velocity': {
#         'magnitude': 10.0,
#         'units': '10.0 meter / second',
#         'value': '10.0 meter / second'
#     },
#     'iterations': 500
# }}
```

Set `standardize=True` to convert all quantities to MKS base units before serializing. Legacy `unum` values are converted to pint first.

## JSONToConfiguration -- Restore Units

Convert a serialized dict back into a configuration with live pint `Quantity` objects:

```python
from hera.utils import JSONToConfiguration

restored = JSONToConfiguration(json_dict)
print(restored["solver"]["velocity"])  # 10.0 meter / second
```

Strings that pint can parse (such as `'10.0 meter / second'`) become `Quantity` objects; other strings, numbers and booleans are returned unchanged. A `{magnitude, units, value}` sub-dict produced by `splitUnits=True` is turned back into a single `Quantity`.

Options:

- `returnStandardize=True` -- for split-unit sub-dicts, rebuild the quantity from `units` (the standardized form) instead of `value`.
- `returnUnum=True` -- return unum objects instead of pint (legacy support; needs `unum`).

## stripConfigurationUnits -- Keep Only Magnitudes

`stripConfigurationUnits(config, returnStandardize=False, ignoreStandardization=[])` replaces every quantity with its magnitude. With `returnStandardize=True` the quantities are converted to MKS first, except under keys listed in `ignoreStandardization`:

```python
from hera.utils import ureg, stripConfigurationUnits

stripConfigurationUnits({"speed": 36 * ureg.km / ureg.hour}, returnStandardize=True)
# {'speed': 10.0}
```

## JSONVariations -- Parameter Sweeps

Generate the Cartesian product of parameter variations applied to a base configuration. Each "variation group" is a dict mapping JSON paths to lists of values. Parameters within one group change together (all lists in a group must have the same length, otherwise `ValueError` is raised); the Cartesian product is taken across groups. `JSONVariations` is a generator that yields a new deep copy of the base for each combination.

```python
from hera.utils import JSONVariations

base = {
    "mesh": {"resolution": 100, "layers": 5},
    "physics": {"viscosity": 0.01},
}

variations = [
    # Group 1: resolution and layers change together
    {"mesh.resolution": [100, 200], "mesh.layers": [5, 10]},
    # Group 2: viscosity varies independently
    {"physics.viscosity": [0.01, 0.001, 0.0001]},
]

for variant in JSONVariations(base, variations):
    res = variant["mesh"]["resolution"]
    visc = variant["physics"]["viscosity"]
    print(f"resolution={res}, viscosity={visc}")

# resolution=100, viscosity=0.01
# resolution=100, viscosity=0.001
# resolution=100, viscosity=0.0001
# resolution=200, viscosity=0.01
# resolution=200, viscosity=0.001
# resolution=200, viscosity=0.0001
```

## compareJSONS -- Side-by-Side Comparison

Compare multiple JSON configurations and see only the fields that differ:

```python
from hera.utils import compareJSONS

config_a = {"solver": {"dt": 0.1, "maxIter": 100}}
config_b = {"solver": {"dt": 0.2, "maxIter": 100}}

diff = compareJSONS(run1=config_a, run2=config_b)
print(diff)
# Shows a DataFrame with solver.dt differing between run1 and run2
```

```
datasetName            run1  run2
parameterNameFullPath
solver.dt               0.1   0.2
```

The JSONs are flattened with `convertJSONtoPandas` (list items become `name_0`, `name_1`, ...) and compared with [`compareDataframeConfigurations`](dataframes.md). Pass `longFormat=True` to get long-format output. Use `changeDotToUnderscore=True` to make the parameter names compatible with `pandas.query()`.

## diffJSONS -- JSON Patch Between Two Documents

For exactly two documents, `diffJSONS(source, target)` returns the RFC 6902 JSON Patch that turns `source` into `target` (an empty list when they are equal):

```python
from hera.utils import diffJSONS

diffJSONS({"a": 1}, {"a": 2})
# [{'op': 'replace', 'path': '/a', 'value': 2}]
```

`diffJSONS` imports the `jsonpatch` package, which is not pinned in `requirements.txt` or `setup.py` (it usually arrives as a transitive dependency); if it is missing, run `pip install jsonpatch`.

## loadJSON and convertJSONtoPandas

- `loadJSON(jsonData)` accepts a file path, a JSON string, a file-like object, a dict or a list and returns the loaded object.
- `convertJSONtoPandas(jsonData, nameColumn="parameterNameFullPath", valueColumn="value")` flattens a JSON (path, string or dict) into a two-column DataFrame of JSON paths and values.

## setJSONPath -- Modify Nested Values by Path

Set values deep inside a nested dict using dot-separated JSON paths:

```python
from hera.utils import setJSONPath

base = {"a": {"b": 1}, "c": {"d": 2}}
updated = setJSONPath(base, {"a.b": 99})
print(updated)
# {'a': {'b': 99}, 'c': {'d': 2}}
```

The original dict is not modified unless you pass `inPlace=True`. Paths may be given with or without the leading `$.`. Each path must already exist in `base`; a missing path raises `IndexError`.
