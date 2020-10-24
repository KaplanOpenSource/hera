"""
==========================
plot contourf by seasons
==========================
"""
from hera.measurements.meteorology.lowfreqdata import lowfreqDataLayer,seasonplots
import seaborn # import for documentation purose
seaborn.set() # for documentation purose
import os # import for documentation purose
import matplotlib.pyplot as plt



###############################################################################
# First, we read data from the data base:
#
# [For documentation purposes, we will read a sample data (using getDocFromFile function) and create a document-like object].

stationName= 'AVNE_ETAN'
station=lowfreqDataLayer.getStationDataFromFile(pathToData=stationName,
                                                     fileFormat=lowfreqDataLayer.JSONIMS,
                                                     storeDB =False)

data=station.getData()


###############################################################################
#  Now we can plot:

plotField='TD'
seasonplots.plotProbContourf_bySeason(data,plotField)

###############################################################################
# The user can override the default ax settings, colors, labeles and more.
# See the full documentation of the plotProbContourf function for details.
#
# Here is an example:

ax_functions_properties=dict(set_ylim=[0,45],
                   set_yticks=[x for x in range(0, 46, 5)],
                   set_ylabel= 'Temperature [°c]'
                  )
withLabels=False
Cmapname='coolwarm'

seasonplots.plotProbContourf_bySeason(data,plotField,ax_functions_properties=ax_functions_properties,withLabels=withLabels,Cmapname=Cmapname)

