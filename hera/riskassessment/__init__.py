import sys
import importlib

# All heavy symbols are deferred to first attribute access so that
# `from hera.riskassessment import CLI` is instant.
_LAZY_SYMBOLS = {
    "Agent":                  ".agents.Agents",
    "RiskToolkit":            ".riskToolkit",
    "ProtectionPolicy":       ".protectionpolicy.ProtectionPolicy",
    "thresholdGeoDataFrame":  ".agents.effects.thresholdGeoDataFrame",
    "getRiskAreaAlgorithm":   ".analysis.riskAreas",
    "casualtiesPlot":         ".presentation.casualtiesFigs",
}

_singletons = {}  # AgentHome, casualtiesPlots — created on first access


def __getattr__(name):
    if name in _LAZY_SYMBOLS:
        mod = importlib.import_module(_LAZY_SYMBOLS[name], package=__name__)
        value = getattr(mod, name)
        globals()[name] = value
        return value
    if name == "AgentHome":
        RiskToolkit = __getattr__("RiskToolkit")
        val = RiskToolkit("")
        globals()["AgentHome"] = val
        return val
    if name == "casualtiesPlots":
        casualtiesPlot = __getattr__("casualtiesPlot")
        val = casualtiesPlot()
        globals()["casualtiesPlots"] = val
        return val
    raise AttributeError(f"module 'hera.riskassessment' has no attribute {name!r}")
