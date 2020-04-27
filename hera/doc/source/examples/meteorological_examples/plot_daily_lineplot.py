"""
=====================
plot single daily lineplot
=====================
"""
from hera import meteo
import seaborn
seaborn.set()

###############################################################################
# First, we read data from the data base:

projectName='IMS_Data'
stationName= 'AVNE ETAN'
station=meteo.IMS_datalayer.getDocFromDB(projectName=projectName, StationName=stationName)
data=station[0].getData()

###############################################################################
# Then, select a column from the data to plot
#
# Let's print the available columns and see what the column we selected stands for:

print(data.columns)
print( 'The selected column stands for %s' %station[0].desc['RH'])

###############################################################################
#  Now we can plot:

plotField='RH'
datetetoplot='2015-02-26'

meteo.IMS_dailyplots.dateLinePlot(data,plotField=plotField,date=datetetoplot)

###############################################################################
# The user can override the default ax settings, colors, labeles and more.
# See the full documentation of the plotProbContourf function for details.
#
# Here is an example:

ax_functions_properties=dict(set_ylim=[0,100],
                   set_yticks=[x for x in range(0, 101, 5)]
                  )
line_properties=dict(color='k')

meteo.IMS_dailyplots.dateLinePlot(data,plotField=plotField,date=datetetoplot,ax_functions_properties=ax_functions_properties,line_properties=line_properties)
