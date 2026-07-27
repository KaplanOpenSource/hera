"""
Lazy package init for hera.utils.

Logging helpers stay eager (cheap and needed during `import hera`). Every other
submodule — unitHandler (pint), jsonutils (pandas + unum transitively), query
(pandas), matplotlibCountour (shapely + geopandas), angle, zipUtils — is loaded
on first attribute access via PEP 562 `__getattr__`. This keeps `import hera`
from pulling pint/pandas/shapely/geopandas when the caller does not need them.
"""

import importlib

from hera.utils.logging import (
    with_logger,
    initialize_logging,
    get_logger,
    getClassLogger,
    get_classMethod_logger,
    add_formatter,
    add_FileHandler,
)

_LAZY_SUBMODULES = (
    "unitHandler",
    "jsonutils",
    "query",
    "matplotlibCountour",
    "angle",
    "zipUtils",
)

# Public names resolvable via __getattr__. Listed explicitly so `from hera.utils
# import *` works without forcing any of the lazy submodules to load at star-
# import time (Python calls getattr for each name in __all__, which then routes
# through __getattr__ and only imports the submodule that actually owns it).
# Only names that are always present (regardless of optional deps like unum)
# are listed here. Optional symbols (e.g. pintToUnum, NameConflictError) are
# still reachable via `from hera.utils import <name>` when the optional
# dependency is installed — they just don't star-import.
__all__ = [
    # logging (already eager above)
    "with_logger", "initialize_logging", "get_logger",
    "getClassLogger", "get_classMethod_logger",
    "add_formatter", "add_FileHandler",
    # unitHandler (always present — unum-only names omitted)
    "ureg", "celsius", "K",
    "Unit", "UnitRegistry", "Quantity",
    "UndefinedUnitError", "DimensionalityError",
    "Unum", "unumSupport",
    "tonumber", "tounit", "tounum",
    # jsonutils
    "compareJSONS", "ConfigurationToJSON", "JSONToConfiguration",
    "stripConfigurationUnits", "loadJSON",
    "processJSONToPandas", "convertJSONtoPandas",
    "setJSONPath", "JSONVariations", "JSONvariationItem",
    # query
    "dictToMongoQuery",
    # matplotlibCountour
    "standardize_polygon", "toGeopandas",
    # angle
    "toMeteorologicalAngle", "toMathematicalAngle", "toAzimuthAngle",
    # zipUtils
    "zip_items", "list_json_files_in_zip",
]

_NOT_FOUND = object()


def __getattr__(name):
    if name.startswith("_"):
        raise AttributeError(f"module 'hera.utils' has no attribute {name!r}")
    for submod_name in _LAZY_SUBMODULES:
        try:
            mod = importlib.import_module(f"hera.utils.{submod_name}")
        except ImportError:
            continue
        value = getattr(mod, name, _NOT_FOUND)
        if value is not _NOT_FOUND:
            globals()[name] = value
            return value
    raise AttributeError(f"module 'hera.utils' has no attribute {name!r}")


def __dir__():
    return sorted(set(__all__) | set(globals()))
