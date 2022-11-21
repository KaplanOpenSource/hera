import pandas
import matplotlib.pyplot as plt
import math

class presentation():

    _datalayer = None
    _analysis = None

    @property
    def datalayer(self):
        return self._datalayer

    @property
    def analysis(self):
        return self._analysis

    def __init__(self, dataLayer, Analysis):

        self._datalayer =  dataLayer
        self._analysis = Analysis

    def variableAgainstDistance(self, data, variable, height, style="plot", colors=["red", "blue"], signedColors=["blue", "orange", "green", "red"],
                                signedDists=None, labels=None, topography=True, ax=None):
        """
        Plots the values of a variable and the terrain height in a slice along the distance downwind.
        Params:
            data: The data of the slice (pandas DataFrame)
            variable: The name of the column of the variable (string)
            colors: The colors of the plot, default red for the variable and blue for the terrain (list of two strings)
            signedDists: Default is None. Else, a list of distances that should be signed using a dashed line. (list of floats)
            signedColors: A list of colors to use for the signed distances (list of strings)
            labels: A list of labels for the two y axis and the x axis, default is ["Distance Downwind", variable, "Terrain"] (list of three strings)
            topography: If true, plots the terrain height.
            ax: ax in which to plot, default is None.
        Return: ax

        """
        if ax is None:
            fig, ax = plt.subplots()
        else:
            plt.sca(ax)

        dataHeight = self.analysis.regularizeTimeSteps(data=data, fieldList=[variable], coord1="distance",coordinateType="constantHeight",heightsList=[height]).reset_index()
        data = data.sort_values(by="distance").reset_index().drop(0)
        dataHeight = dataHeight.sort_values(by="distance").drop(0)
        labels = labels if labels is not None else ["Distance Downwind (m)", variable, "Terrain (m)"]

        getattr(ax, style)(dataHeight.distance, dataHeight[variable], color=colors[0])
        ax.tick_params(axis='y', labelcolor=colors[0])
        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1], color=colors[0])
        if topography:
            ax2 = ax.twinx()
            ax2.plot(data.distance, data.terrain, color=colors[1])
            ax2.tick_params(axis='y', labelcolor=colors[1])
            ax2.set_ylabel(labels[2], color=colors[1])
        if signedDists is not None:
            for i in range(len(signedDists)):
                ax2.plot([signedDists[i],signedDists[i]], [data.z.min(), data.z.max()], color=signedColors[i], linestyle="--")
        return ax

    def UinLocations(self, data, points, style="plot", colors=["blue", "red"], labels=["Distance Downwind (m)", "Velocity (m/s)", "Height (m)"],ax=None):

        if ax is None:
            fig, ax = plt.subplots()
        else:
            plt.sca(ax)
        data = data.sort_values(by=["distance", "heightOverTerrain"])
        ax.plot(data.distance, data.terrain, zorder=10, color=colors[0])
        ax.set_ylim(data.z.min(), data.z.max())
        xticks = [i for i in range(int(data.Velocity.max()+2))]
        ax.set_ylabel(labels[2])
        ax.set_xlabel(labels[0])
        for point in points:
            data = data.append(pandas.DataFrame([{"distance":point}]))
        data = data.sort_values(by="distance").interpolate(by="distance")
        terrain = data.query("distance in @points")
        dataPoints = self.analysis.regularizeTimeSteps(data=data, fieldList=["Velocity"], coord1="distance",coordinateType="Sigma", terrain=terrain, n=100).reset_index()
        for point in points:
            axins = ax.inset_axes([point, data.loc[data.distance==point].terrain.mean(), data.distance.max()/10,
                                   data.z.max()-data.loc[data.distance==point].terrain.mean()], transform=ax.transData)
            getattr(axins, style)(dataPoints.loc[dataPoints.distance==point].sort_values(by=["distance", "sigma"]).Velocity,
                    dataPoints.loc[dataPoints.distance==point].sort_values(by=["distance", "sigma"]).z, color=colors[1], zorder=0)
            axins.set_ylim(data.loc[data.distance==point].terrain.mean(), dataPoints.loc[dataPoints.distance==point].z.max())
            axins.get_yaxis().set_visible(False)
            axins.xaxis.set_ticks_position("top")
            axins.xaxis.set_label_position("top")
            axins.set_xlim(0, data.Velocity.max()+0.5)
            axins.set_xticks(xticks)
            axins.set_xlabel(labels[1])
            ax.scatter(point, data.loc[data.distance==point].terrain.mean(), color=colors[1], zorder=15)
        return ax