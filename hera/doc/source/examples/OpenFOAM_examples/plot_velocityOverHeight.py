"""
==============================================
Plot the velocity over height.
==============================================
"""

#######################
# This example shows how to plot the change of the velocity with height in different points,
# in addition to the terrain height along the distance.
#
# In order to do so, first of all the data should be arranged.
# We load data of a slice of the domain from the database (look for example in simulation OpenFOAM to see how do add the data to the database).
# The type of such data is OFsimulation. In the presented project, the desired data is the second data of that kind in the database.
#
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
# The arrangement adds columns that hold the velocity magnitude, the height ofer the terrain
# and the distance downwind from the domain's corner.

data = OFdatalayer.analysis.arrangeSlice(data=data, ydir=False)

#######################
# ydir = False means that the wind component in the y direction is negative.
#
# The plot is plotted using the presentation layer.
# The colors and labels can be changed from the default values by addressing parameters called "colors" and "labels".

d = OFdatalayer.presentation.UinLocations(data=data, points=[1000,2000,3000])