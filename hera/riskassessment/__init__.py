"""
Lazy package init for hera.riskassessment.

Every public name here pulls in pint (via hera.utils.unitHandler) and geopandas,
which costs ~3s. The CLI only needs argparse wiring at import time, so the
symbols are resolved on first attribute access via PEP 562 __getattr__ —
the same pattern used in hera.utils.
"""

import importlib

_LAZY_NAMES = {
    "Agent": (".agents.Agents", "Agent"),
    "thresholdGeoDataFrame": (".agents.effects.thresholdGeoDataFrame", "thresholdGeoDataFrame"),
    "getRiskAreaAlgorithm": (".analysis.riskAreas", "getRiskAreaAlgorithm"),
    "casualtiesPlot": (".presentation.casualtiesFigs", "casualtiesPlot"),
    "ProtectionPolicy": (".protectionpolicy.ProtectionPolicy", "ProtectionPolicy"),
    "RiskToolkit": (".riskToolkit", "RiskToolkit"),
}

__all__ = list(_LAZY_NAMES)


def __getattr__(name):
    entry = _LAZY_NAMES.get(name)
    if entry is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module, attribute = entry
    value = getattr(importlib.import_module(module, __name__), attribute)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(_LAZY_NAMES))
