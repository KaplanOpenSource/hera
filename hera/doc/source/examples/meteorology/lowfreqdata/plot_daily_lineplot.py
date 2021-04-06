"""
==========================
plot single daily lineplot
==========================
"""

from hera.measurements.meteorology.lowfreqdata import lowfreqDataLayer,dailyplots
import seaborn # import for documentation purose
seaborn.set() # for documentation purose
import os # import for documentation purose
import matplotlib.pyplot as plt

###############################################################################
# First, we read data from the data base:
#
# For documentation purposes, we will read a sample data (using getDocFromFile function)
# and create a document-like object.

stationName= 'AVNE_ETAN'
station=lowfreqDataLayer.getStationDataFromFile(pathToData=stationName,
                                                     fileFormat=lowfreqDataLayer.JSONIMS,
                                                     storeDB =False)

data=station.getData()

###############################################################################
# Then, select a column from the data to plot
#
# Let's print the available columns and see what the column we selected stands for:
#
#
print(data.columns)
#
#
###############################################################################
#  Now we can plot:

plotField='RH'
datetetoplot='2015-02-26'


dailyplots.dateLinePlot(data,plotField=plotField,date=datetetoplot)

##############################################################################
# The user can override the default ax settings, colors, labeles and more.
# See the full documentation of the plotProbContourf function for details.
#
# Here is an example:

ax_functions_properties=dict(set_ylim=[0,100],
                   set_yticks=[x for x in range(0, 101, 5)]
                  )
line_properties=dict(color='k')

dailyplots.dateLinePlot(data,
                                  plotField=plotField,
                                  date=datetetoplot,
                                  ax_functions_properties=ax_functions_properties,
                                  line_properties=line_properties)

