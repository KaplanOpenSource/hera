"""
=====================
plot single contourf
=====================
"""
from hera import meteo

###############################################################################
# First, we read data from the data base:

projectName='IMS_Data'
stationName= 'AVNE ETAN'
station=meteo.IMS_datalayer.getDocFromDB(projectName=projectName, StationName=stationName)
data=station[0].getData()
print(data)

###############################################################################
# Then, select a column from the data to plot
#
# Let's print the available columns and see what the column we selected stands for:

print(data.columns)
print( 'The selected column stands for %s' %station[0].desc['WD'])

###############################################################################
#  Now we can plot:

plotField='WD'
meteo.IMS_dailyplots.plotProbContourf(data,plotField)

###############################################################################
# The user can override the default ax settings, colors, labeles and more.
# See the full documentation of the plotProbContourf function for details.
#
# Here is an example:

ax_functions_properties=dict(set_ylim=[0,360],
                             set_yticks=[x for x in range(0, 361, 30)]
                            )
meteo.IMS_dailyplots.plotProbContourf(data,plotField,ax_functions_properties=ax_functions_properties)

