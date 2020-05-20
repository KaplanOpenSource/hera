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

from hera import datalayer
data = datalayer.Measurements.getDocuments(projectName="Example", type="OFsimulation")[1].getData().compute()

#######################
# Now, we use a function of hera openfoam in order to arrange the data.
# The arrangement adds columns that hold the velocity magnitude, the height ofer the terrain
# and the distance downwind from the domain's corner.

from hera import openfoam
op = openfoam.dataManipulations(projectName="Example")
data = op.arrangeSlice(data=data, ydir=False)

#######################
# ydir = False means that the wind component in the y direction is negative.
#
# We have to find points in which there is enough data.
# In order to do so, we may use the next function:

optionalLocations = op.findDetaiedLocations(data)

#########################
# The function returns a list of specific distances downwinds which have data of velocities in many different heights.
# One may choose any values from the list. We chose three values randomly.
#
# The plot is plotted using the plotting model.
# The colors and labels can be changed from the default values by addressing parameters called "colors" and "labels".

plotting = openfoam.Plotting()
plotting.UinLocations(data=data, points=[optionalLocations[9], optionalLocations[56], optionalLocations[80]])