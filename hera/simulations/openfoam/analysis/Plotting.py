from .... import GIS
import matplotlib.pyplot as plt

class Plotting():

    _data = None

    def __init__(self, data):

        self._data = data


    def variableAgainstDistance(self, variable, colors=["red", "blue"], signedColors=["blue", "orange", "green", "red"],
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

        data = self._data.sort_values(by="distance")
        labels = labels if labels is not None else ["Distance Downwind (m)", variable, "Terrain (m)"]

        ax.plot(data.distance, data[variable], color=colors[0])
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

    def variableAgainstHeight(self, variable, nOfPoints=4, labels=None, ax=None):

        if ax is None:
            fig, ax = plt.subplots()
        else:
            plt.sca(ax)

        labels = labels if labels is not None else ["Height (m)", variable]

        optional = []
        dists = []
        for d in self._data.distance.drop_duplicates():
           # if len(data.query("distance==@d")) > 20:
           if len(self._data.query("distance==@d and heightOverTerrain<10")) > 3 and len(
                   self._data.query("distance>@d-0.5 and distance<@d+0.5 and heightOverTerrain>10")) > 10:
               optional.append(d)
        delta = int(len(optional)/(nOfPoints+1))
        optional.sort()
        data = self._data.sort_values(by="heightOverTerrain")
        for i in range(nOfPoints):
            dist = optional[delta*(i+1)]
            dists.append(dist)
            ax.plot(data.query("distance>@dist-0.5 and distance<@dist+0.5").heightOverTerrain, data.query("distance>@dist-0.5 and distance<@dist+0.5")[variable], label=int(dist))
        ax.legend()
        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])

        return dists

    def UinLocations(self, points, colors=["blue", "red"], labels=["Distance Downwind (m)", "Velocity (m/s)", "Height (m)"],ax=None):

        if ax is None:
            fig, ax = plt.subplots()
        else:
            plt.sca(ax)
        data = self._data.sort_values(by=["distance", "heightOverTerrain"])
        ax.plot(data.distance, data.terrain, zorder=10, color=colors[0])
        ax.set_ylim(data.z.min(), data.z.max())
        xticks = [i for i in range(int(data.Velocity.max()+2))]
        ax.set_ylabel(labels[2])
        ax.set_xlabel(labels[0])
        for point in points:
            axins = ax.inset_axes([point, data.loc[data.distance==point].terrain.mean(), data.distance.max()/10,
                                   data.z.max()-data.loc[data.distance==point].terrain.mean()], transform=ax.transData)
            axins.plot(data.query("distance>@point-0.5 and distance<@point+0.5").sort_values(by=["distance", "heightOverTerrain"]).Velocity,
                    data.query("distance>@point-0.5 and distance<@point+0.5").sort_values(by=["distance", "heightOverTerrain"]).z, color=colors[1], zorder=0)
            axins.set_ylim(data.loc[data.distance==point].terrain.mean(), data.query("distance>@point-0.5 and distance<@point+0.5").z.max())
            axins.get_yaxis().set_visible(False)
            axins.xaxis.set_ticks_position("top")
            axins.xaxis.set_label_position("top")
            axins.set_xlim(0, data.Velocity.max()+0.5)
            axins.set_xticks(xticks)
            axins.set_xlabel(labels[1])
            ax.scatter(point, data.loc[data.distance==point].terrain.mean(), color=colors[1], zorder=15)