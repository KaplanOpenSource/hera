"""
==============================================
Plot a variable at a fixed height along slice
==============================================
"""

#######################
# This example shows how to plot a variable at a fixed height over terrain along the distance downwind,
# in addition to the terrain height along the distance.
#
# In order to do so, first of all the data should be arranged.
# We load data of a slice of the domain from the database (look for example in simulation OpenFOAM to see how do add the data to the database).
# The type of such data is OFsimulation. In the presented project, the desired data is the second data of that kind in the database.
#
# The data may be loaded from the database, as demonstrated below.

from hera import datalayer
#data = datalayer.Measurements.getDocuments(projectName="Example", type="OFsimulation", filter="Slice1)[0].getData().compute()

###########################
# For the documentation, we would load it directly from the resource.

import dask.dataframe
data = dask.dataframe.read_parquet("documentationData/Slice1.parquet").compute()

#######################
# Now, we use a function of hera openfoam in order to arrange the data.
# The arrangement adds columns that hold the velocity magnitude, the height ofer the terrain
# and the distance downwind from the domain's corner.

# from hera import openfoam
# op = openfoam.dataManipulations(projectName="Example")
# data = op.arrangeSlice(data=data, ydir=False)
#
# #######################
# # ydir = False means that the wind component in the y direction is negative.
# #
# # The plot is plotted using a module named plotting.
# # We will plot the velocity at height 30 m above the terrain.
#
# plotting = openfoam.Plotting()
# plotting.variableAgainstDistance(data=data, variable="Velocity", height=30)
#
# #######################
# # The default style of the plot is a lineplot. A scatter plot can be plotted instead:
#
# plotting.variableAgainstDistance(data=data, variable="Velocity", height=30, style="scatter")
#
# #######################
# # The colors and labels can be changed from the default values.
# # In addition, signed distances may be added as in the next example.
# # The colors of the signed lines are default colors that may be changed using the parameter signedColors.
#
# plotting.variableAgainstDistance(data=data, variable="Velocity", height=30, signedDists=[500,1500,2500,3500],
#                                  colors=["orange", "green"], labels=["Distance (m)", "Velocity (m/s)", "Terrain (m)"])
#
