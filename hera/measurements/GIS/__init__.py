from .datalayer.datalayer import GIS_datalayer
from .analytics.Plotting import Plotting
from .utils.STL import convert
from .analytics.dataManipulations import dataManipulations
from .demography.population import population
import sys
version = sys.version_info[0]
if version==2:
    from .synthetic.synthetic import synthetic
dataManipulations = dataManipulations()