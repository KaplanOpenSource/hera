"""
=====================
plot contourf by seasons
=====================
"""
from hera import meteo
import seaborn
seaborn.set()

###############################################################################
# First, we read data from the data base:

projectName='IMS_Data'
stationName= 'AVNE ETAN'
station=meteo.IMS_getDocFromDB(projectName=projectName, StationName=stationName)
data=station[0].getData()
print(data)

###############################################################################
# Then, select a column from the data to plot
#
# Let's print the available columns and see what the column we selected stands for:


print(data.columns)
print( 'The selected column stands for %s' %station[0].desc['TD'])

###############################################################################
#  Now we can plot:

plotfield='TD'

meteo.seasonplots.plotProbContourf_bySeason(data,plotfield)

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

meteo.seasonplots.plotProbContourf_bySeason(data,plotfield,ax_functions_properties=ax_functions_properties,withLabels=withLabels,Cmapname=Cmapname)

