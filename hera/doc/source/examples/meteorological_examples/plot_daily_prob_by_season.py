"""
==========================
plot contourf by seasons
==========================
"""
from hera import meteo
import seaborn
import os
seaborn.set()

###############################################################################
# First, we read data from the data base:

projectName='IMS_Data'
stationName= 'AVNE ETAN'
station=meteo.IMS_datalayer.getDocFromDB(projectName=projectName, StationName=stationName)
data=station[0].getData()
print(data)

###############################################################################
# for documentation, we will read an example data from directory:

station=meteo.IMS_datalayer.getDocFromFile(path=os.path.join("documentationData","AVNE_ETAN"))
data=station[0].getData()
###############################################################################
# Then, select a column from the data to plot
#
# Let's print the available columns and see what the column we selected stands for:


print(data.columns)
#print( 'The selected column stands for %s' %station[0].desc['TD'])

###############################################################################
#  Now we can plot:

plotField='TD'
meteo.IMS_seasonplots.plotProbContourf_bySeason(data,plotField)

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

meteo.IMS_seasonplots.plotProbContourf_bySeason(data,plotField,ax_functions_properties=ax_functions_properties,withLabels=withLabels,Cmapname=Cmapname)

