from .agents.Agents import Agent
from .riskToolkit import RiskToolkit
from .protectionpolicy.ProtectionPolicy import  ProtectionPolicy
from .agents.effects.thresholdGeoDataFrame import thresholdGeoDataFrame

from .analysis.riskAreas import getRiskAreaAlgorithm

from .presentation.casualtiesFigs import casualtiesPlot



AgentHome       = RiskToolkit("")
casualtiesPlots = casualtiesPlot()

