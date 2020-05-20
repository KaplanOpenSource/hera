"""
==============================================
Plot a variable at a fixed height over terrain
==============================================
"""

#######################
# This example shows how to plot a variable at a fixed height over the whole terrain,
# next to a plot of the terrain. We would plot the velocity at height 20 m.
# In order to plot it, we would load data that was made in advance by the function makeClipHeightData
# of openfoam.dataManipulations, as demonstrated in the openfoam tutorial.
# The loaded data is a pandas dataframe that holds x and y coordinates and the velocity at height 20 m.
#
# We would also use the GIS module in order to load the GIS data of the region.

from hera import datalayer
from hera import GIS
import matplotlib.pyplot as plt

GISdatalayer = GIS.GIS_datalayer(projectName="Example", FilesDirectory="")
simulationData = datalayer.Measurements.getDocuments(projectName="Example")[11].getData(usePandas=True)
GISdata = GISdatalayer.getGISDocuments()[0].getData()

#######################
# We would use matplotlib tricontourf in order to plot the velocity.
# The plot is done this way:

levels=[0,1,1.5,1.75,2,2.25,2.5] # velocity levels for the tricontour
fig, ax = plt.subplots(1,2, figsize=(10,15))
CS = ax[1].tricontourf(simulationData.x, simulationData.y, simulationData.Velocity, levels=levels)
cbax1 = fig.add_axes([0.91, 0.3, 0.03, 0.39])
cbax2 = fig.add_axes([0.44, 0.3, 0.03, 0.39])
ax[1].set_aspect('equal', adjustable='box')
GISdata.plot(column="HEIGHT", ax=ax[0])
sm = plt.cm.ScalarMappable(norm=plt.Normalize(vmin=min(GISdata.HEIGHT), vmax=max(GISdata.HEIGHT)))
fig.colorbar(CS, cax=cbax1)
fig.colorbar(sm, cax=cbax2)
plt.subplots_adjust(wspace=0.5)
ax[0].set_title("Terrain Height Above Sea (m)")
ax[1].set_title("Velocity at 20 m \nAbove Terrain (m/s)")
ax[0].set_yticks([737000,738000,739000])
ax[1].set_yticks([737000,738000,739000])
ax[0].set_xticks([197000,198000,199000])
ax[1].set_xticks([197000,198000,199000])
ax[0].set_xlabel("ITM")
ax[0].set_ylabel("ITM")
ax[1].set_xlabel("ITM")
ax[1].set_ylabel("ITM")