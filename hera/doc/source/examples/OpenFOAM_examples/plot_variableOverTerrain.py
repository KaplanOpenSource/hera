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
# The data may be loaded from the database, as demonstrated below.

from hera.simulations.openfoam.postProcess.meshUtilities.datalayer import datalayer
OFdatalayer = datalayer(projectName="Documentation")
#data = datalayer.getDocuments(filter="Slice1")[0].getData().compute()

###########################
# For the documentation, we would load it directly from the resource.

import dask.dataframe
data = dask.dataframe.read_parquet("documentationData/Slice1.parquet").compute()

#######################
# Now, we use a function of the analysis layer in order to arrange the data.
# It adds columns that hold the velocity magnitude, the height over the terrain
# and the distance downwind from the domain's corner.

data = OFdatalayer.analysis.arrangeSlice(data=data, ydir=False)

#######################
# ydir = False means that the wind component in the y direction is negative.
#
# The plot is plotted using the resentation layer.
# We will plot the velocity at height 30 m above the terrain.

OFdatalayer.presentation.variableAgainstDistance(data=data, variable="Velocity", height=30)

#######################
# The default style of the plot is a lineplot. A scatter plot can be plotted instead:

OFdatalayer.presentation.variableAgainstDistance(data=data, variable="Velocity", height=30, style="scatter")

#######################
# The colors and labels can be changed from the default values.
# In addition, signed distances may be added as in the next example.
# The colors of the signed lines are default colors that may be changed using the parameter signedColors.

OFdatalayer.presentation.variableAgainstDistance(data=data, variable="Velocity", height=30, signedDists=[500,1500,2500,3500],
                                 colors=["orange", "green"], labels=["Distance (m)", "Velocity (m/s)", "Terrain (m)"])

