"""
=====================
plot single contourf
=====================
"""

from hera.measurements.meteorology.lowfreqdata import lowfreqDataLayer,dailyplots
import seaborn # import for documentation purose
seaborn.set() # for documentation purose
import os # import for documentation purose
import matplotlib.pyplot as plt

from hera.measurements.meteorology.lowfreqdata.analysis import seasonsdict

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

plotField='WD'
dailyplots.plotProbContourf(data,plotField)

###############################################################################
# The user can override the default ax settings, colors, labeles and more.
# See the full documentation of the plotProbContourf function for details.
#
# Here is an example:

ax_functions_properties=dict(set_ylim=[0,360],
                             set_yticks=[x for x in range(0, 361, 30)]
                            )
dailyplots.plotProbContourf(data,plotField,ax_functions_properties=ax_functions_properties)

plt.show()