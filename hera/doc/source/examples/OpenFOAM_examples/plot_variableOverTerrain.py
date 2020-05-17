"""
==============================================
Plot a variable at a fixed height over terrain
==============================================
"""

#######################
# This example shows how to plot a variable at a fixed height over terrain along the distance downwind,
# in addition to the terrain height along the distance.
#
# In order to do so, first of all the data should be arranged.
# We load data of a slice of the domain from the database (look for example in simulation OpenFOAM to see how do add the data to the database).
# The type of such data is OFsimulation. In the presented project, the desired data is the first data of that kind in the database.

from hera import datalayer
data = datalayer.Measurements.getDocuments(projectName="Examples", type="OFsimulation")[0].getData().compute()

#######################
# Now, we use a function of hera openfoam in order to arrange the data.
# The arrangement adds columns that hold the velocity magnitude, the height ofer the terrain
# and the distance downwind from the domain's corner.

from hera import openfoam
op = openfoam.dataManipulations()
data = op.arrangeSlice(data=data, ydir=False)

#######################
# ydir = False means that the wind component in the y direction is negative.
# The goal is to plot the velocity at height of 10 m above the terrain.
# Therefore, we select part of the data.

dataForPlot = data.query("heightOverTerrain > 9.8 and heightOverTerrain < 10.2")

#######################
# The plot is plotted using a module named plotting.

plotting = openfoam.Plotting(dataForPlot)
plotting.variableAgainstDistance(variable="Velocity")

#######################
# The colors and labels can be changed from the default values.
# In addition, signed distances may be added as in the next example.
# The colors of the signed lines are default colors that may be changed using the parameter signedColors.

plotting.variableAgainstDistance(variable="Velocity", signedDists=[1000,3000,5000,7000], colors=["orange", "green"], labels=["Distance (m)", "Velocity (m/s)", "Terrain (m)"])

