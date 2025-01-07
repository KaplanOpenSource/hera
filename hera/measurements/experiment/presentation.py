import matplotlib.pyplot as plt
import seaborn
import matplotlib.colors as mcolors
import numpy
from hera.measurements.GIS.utils import WSG84,ITM,convertCRS
import pandas as pd
import numpy as np
from hera import toolkitHome
import jinja2
import os

class experimentPresentation:
    """
        A basic implementation of the presentation layer for the experiment.

        Has 3 different plots:

        - setup plots/tables:
                Plotting or tabulating the location and the positions of devices.

        - technical plots/tables:
                Plotting or tabulating the technical aspects of the data (frequency and ect).

        - data   plots/tables:
                Plotting or tabulating the data itself.

        Hence all the procedures are:

        [setup|technical|data][Plot|Table]<the name>.

    """
    NAMES_SHORT = "short"  # just plot the number
    NAMES_LONG = "long"  # plot the entire name.

    ########################
    ###
    ###  Init and private functions
    ###
    #########################

    _saveFigures = None
    _savePath = None

    _datalayer = None
    _analysis = None

    _cmap = None

    @property
    def cmap(self):
        return self._cmap
    @property
    def saveFigures(self):
        return self._saveFigures

    @saveFigures.setter
    def saveFigures(self, value):
        self._saveFigures = value

    @property
    def savePath(self):
        return os.path.abspath(self._savePath)

    @savePath.setter
    def savePath(self, value):
        self._savePath = value

    @property
    def datalayer(self):
        return self._datalayer

    @property
    def analysis(self):
        return self._analysis

    def __init__(self,datalayer,analysis):
        self._datalayer = datalayer
        self._analysis  = analysis
        colors = [
            [1., 0.01176471, 0.01176471],
            [1., 0.01176471, 0.01176471],
            [1., 0.64705884, 0.11764706],
            [1., 0.64705884, 0.11764706],
            [0.15294118, 0.68235296, 0.3764706],
            [0.15294118, 0.68235296, 0.3764706]
        ]
        self._cmap_levels = [0, 0.15, 0.25, 0.75, 0.85, 1]
        self._cmap = self._get_continuous_cmap(colors, self._cmap_levels)


    ########################
    ###
    ###  Setup plots
    ###
    #########################


    def setupPlotOrigin(self, ax=None, s=50):
        if ax is None:
            fig, ax = plt.subplots(1, 1)

        plt.scatter(0, 0, s=s, marker='X', c='xkcd:lavender')

        return ax

    def setupPlotImage(self, imageName, ax=None, xlabel=True, ylabel=True, withGrid = True,plt_kwargs=dict(),majorLocator=10):
        """
                Plots an image of the experiment.
                If axes is not supplied, creates the figure and the axes. 

        Parameters
        ----------
        imageName: str
                The name of the image
        ax: Axis
                Figure axes.

        xlabel: bool
                If true adds the label x [m] to x-axis
        ylabel: bool
                If true adds the label y [m] to y-axis

        withGrid : bool
                If true plot the grid.

        plt_kwargs: dict
                a dict of dicts with matplotlib.pyplot function -> args.
                The code tranvse and activates all.

                For example, in order to import xlim and ylim

                    > plt_kwargs=dict(xlim=(1,2),ylim=(4,5))

        Returns
        -------
            Axes
        """

        if ax is None:
            fig, ax = plt.subplots(1, 1)
        else:
            plt.sca(ax)

        metadata = self.datalayer.getImageMetadata(imageName)
        img = self.datalayer.getImage(imageName)

        plt.imshow(img, extent=[metadata['left'],
                                metadata['right'],
                                metadata['lower'],
                                metadata['upper']]
                   )

        for fnc,values in plt_kwargs.items():
            plt_func = getattr(plt,fnc)
            plt_func(values)


        if xlabel:
            plt.xlabel("x [m]")

        if ylabel:
            plt.ylabel("y [m]")

        if withGrid :
            plt.grid(b=True, which='major', linestyle='-')
            plt.grid(b=True, which='minor', linestyle='--')

            yLocatorSpacing = (metadata['upper']-metadata['lower'])//majorLocator
            xLocatorSpacing = (metadata['right'] - metadata['left']) // majorLocator

            locatorSpacing = numpy.min([yLocatorSpacing,xLocatorSpacing])

            ax.xaxis.set_major_locator(plt.MultipleLocator(locatorSpacing))
            # ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
            ax.yaxis.set_major_locator(plt.MultipleLocator(locatorSpacing))
            # ax.yaxis.set_minor_locator(plt.MultipleLocator(1))
            #plt.minorticks_on()

        plt.axis("equal")

        return ax

    def _plotEntityLocationScatter(self, entityTypeName, trialSet, trialName, status, floorName, ax=None,
                                   plotNameMode=None, scatter_kw=dict()):
        """
            Plots a scatter plot of the entity.
            the scatter_kw is used to send more arguments to the scater function inside.

        :param entityTypeName: str
                The type of the entity to plot.
        :param trialSet: str
                The name of the trialset
        :param trialName: str
                The name of the trial.

        :param status: str
                Is it design or deploy.
        :param floorName: is it Concourse or platform.

        plotName: None, NAME_SHORT, NAME_LONG
            if not None, plot the names of the devices (can be short, long)

        scatter_kw: dict
            Parameters to pass to the scatter.

        :param ax: axis.
            The axes to draw on.
        :return:
            axis
        """

        data = self.datalayer.experimentSetup.trialSet[trialSet][trialName].entitiesTable(status).query("locationName==@floorName and entityType==@entityTypeName")

        if ax is None:
            fig, ax = plt.subplots(1, 1)

        if self._entityMarkers.get(entityTypeName) is not None:
            scatterProperties = self._entityMarkers.get(entityTypeName)['scatter']

            scatter_kw_final = dict(scatterProperties['attrs'])
            scatter_kw_final.update(scatter_kw)

            if scatterProperties['color_handler'] is not None:
                scatter_kw_final = scatterProperties['color_handler'](data, scatter_kw_final)
        else:
            scatter_kw_final = dict()

        plt.scatter(data.longitude, data.latitude, **scatter_kw_final)

        if plotNameMode is not None:
            for indx, itr in data.iterrows():
                if plotNameMode == self.NAMES_SHORT:
                    pltname = itr.entityName.split(" ")[-1]

                elif plotNameMode == self.NAMES_LONG:
                    pltname = itr.entityName
                else:
                    raise ValueError(
                        f"plotNameMode can be either {self.NAMES_SHORT} or {self.NAMES_LONG}. Got {plotNameMode}")

                textShiftDict = {FLOOR_PLATFORM: 1, FLOOR_CONCOURSE: 2}

                plt.text(itr.longitude - 1, itr.latitude - textShiftDict[floorName], s=pltname)

        return ax

    def _plotEntityLocationNames(self, entityTypeName, trialSet, trialName, status, floorName, ax=None,
                                   plotNameMode=None, text_kw=dict()):
        """

        :param entityTypeName: str
                The type of the entity to plot.
        :param trialSet: str
                The name of the trialset
        :param trialName: str
                The name of the trial.

        :param status: str
                Is it design or deploy.
        :param floorName: is it Concourse or platform.

        plotName: None, NAME_SHORT, NAME_LONG
            if not None, plot the names of the devices (can be short, long)

        scatter_kw: dict
            Parameters to pass to the scatter.

        :param ax: axis.
            The axes to draw on.
        :return:
            axis
        """
        entities = self.datalayer.experimentSetup.trialSet[trialSet][trialName].entitiesTable(status)
        if len(entities) > 0:
            data = entities.query("locationName==@floorName and entityType==@entityTypeName")
        else:
            data = pandas.DataFrame()


        if ax is None:
            fig, ax = plt.subplots(1, 1)


        if plotNameMode is not None:
            for indx, itr in data.iterrows():
                if plotNameMode == self.NAMES_SHORT:
                    pltname = itr.entityName.split(" ")[-1]

                elif plotNameMode == self.NAMES_LONG:
                    pltname = itr.entityName
                else:
                    raise ValueError(
                        f"plotNameMode can be either {self.NAMES_SHORT} or {self.NAMES_LONG}. Got {plotNameMode}")

                text_kw.setdefault("size",15)
                text_kw.setdefault("bbox",dict(boxstyle="square",ec=(1., 0.5, 0.5),fc=(1., 0.8, 0.8)))
                plt.text(itr.longitude, itr.latitude, s=pltname,**text_kw)

        return ax



    def _scatter_height_color(self, data, attrMap):
        """
            Adds the heights in data as color in attrmap.
            This map will be sent to scatter.

        :param data:
        :param attrMap:
        :return:
        """
        hts = [float(x) for x in data.height.values]

        if len(hts)==0:

            attrMap['vmin'] = 0
            attrMap['vmax'] = 0

        else:
            attrMap['vmin'] = numpy.min(hts)
            attrMap['vmax'] = numpy.max(hts)

        return attrMap

    ########################
    ###
    ###  Technical plots
    ###
    #########################



    def _get_continuous_cmap(self, cmap_list, float_list=None):
        cdict = dict()

        if float_list is None:
            float_list = list(numpy.linspace(0, 1, len(cmap_list)))

        for num, col in enumerate(['red', 'green', 'blue']):
            col_list = [[float_list[i], cmap_list[i][num], cmap_list[i][num]] for i in range(len(float_list))]
            cdict[col] = col_list
            cmp = mcolors.LinearSegmentedColormap('technical', segmentdata=cdict, N=256)
        return cmp

    def _splitName(self,x):
        if len(x.split(" ")) > 1:
            return int(x.split(" ")[1])
        else:
            return x.split(" ")[0]

    def plotDeviceTypeFunctionality(self,
                                    deviceType,
                                    trialName,
                                    trialSetName=None,
                                    samplingWindow="1min",
                                    equalSquares=False,
                                    ax=None):


        """
            Plots a heatmap of the actual normalized frequency.
            Fills in devices that did not get data with 0.

            Parameters
            ----------
            deviceType : string
                    The type of the device to present.

            trialName : string
                    The name of the trial to show.

            trialSetName : string
                    The name of the trial set.

            samplingWindow : string
                    A time string (of pandas). i.e. '1min' and ect.

            equalSquares : bool
                    If true, then force each square in the heatmap to be of equal size.

            ax : axis
                If None, build a new figure.

            completeTimeSeries : bool
                    If true, and the release has end and start times, then add all the missing time slots
                    with the requested interval.
                    default: True

            completeDevices : bool=True
                    If True, add all the missing devices that were related to the trial.
                    default: True


        """
        # get all the data from the DB. Right now, getting the raw. Maybe to accelerate we will
        # add a flag to try to obtain the data from a cached database.
        trialSetName = self.datalayer.defaultTrialSet if trialSetName is None else trialSetName
        analysisLayer = self.analysis

        pvt = analysisLayer.getDeviceTypeTransmissionFrequencyOfTrial(deviceType,
                                                                      trialName,
                                                                      trialSetName,
                                                                      samplingWindow=samplingWindow,
                                                                      normalize=True,
                                                                      completeTimeSeries=True,
                                                                      completeDevices=True,
                                                                      wideFormat=True)



        pvt = pvt.assign(time=pvt.index)
        pvt = pvt.assign(time=pvt.time.apply(lambda x: x.strftime("%H:%M"))).set_index("time").sort_index(
            ascending=False)
        pvt= pvt.fillna(0)
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(18, 12))
        else:
            plt.sca(ax)

        # Sorting the columns according to the device ID
        experimentSetup = self.datalayer.setup

        allDevices = sorted(list(self.datalayer.entityType[deviceType].keys()), key=lambda x: self._splitName(x))
        cmap = self.cmap
        xticklabels = [f" {self._splitName(x)} " for x in allDevices]
        tmp = seaborn.heatmap(pvt[allDevices],
                               cmap=cmap,
                               vmin=0,
                               vmax=1,
                               cbar_kws=dict(ticks=self._cmap_levels),
                               xticklabels=xticklabels,
                              square=equalSquares,
                               linecolor=[0.39607844, 0.34509805, 0.34509805],
                              # linewidths=0.01)
                              )
        tmp = plt.xlabel("Device name")
        tmp = plt.ylabel("Time")

        # if self.presentation.saveFigures:
        #     figname = os.path.join(self.presentation.savePath, "technical",
        #                            f"{trialName}_{deviceType}_DataFrequency_{trialSetName}_{samplingWindow}.png")
        #     plt.savefig(figname)

        return ax, pvt


    def plotNDIRFrequencyDistribution(self,
                                      trialName,
                                      trialSetName,
                                      ax=None):
        """
            Calculates the cumulative histogram data of the normalized frequency

        Parameters
        ----------
        deviceType
        trialNameOrList
        trialSetName
        samplingWindow
        ax


        Returns
        -------

        """
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        else:
            plt.sca(ax)

        analysisLayer = self.datalayer.analysisLayer
        long = analysisLayer.getDeviceTypeTransmissionFrequency(NDIR, trialName, trialSetName, normalize=True,
                                                                wideFormat=False)

        long.groupby("deviceName").mean().hist(density=True, cumulative=True, ax=ax)
        plt.title("")
        plt.xlim(0, 1)
        plt.xlabel("Normalized frequency")
        plt.ylabel("Fraction of devices")

        if self.presentation.saveFigures:
            figname = os.path.join(self.presentation.savePath, "technical",
                                   f"{trialName}_NDIR_DeviceFrequencyDistribution_{trialSetName}.png")
            plt.savefig(figname)


    def plotMessageFrequencyDistribution(self,
                                         deviceType,
                                         trialName,
                                         trialSetName,
                                         ax=None):
        """
            Calculates the cumulative histogram data of the normalized frequency

        Parameters
        ----------
        deviceType
        trialNameOrList
        trialSetName
        samplingWindow
        ax


        Returns
        -------

        """
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        else:
            plt.sca(ax)

        analysisLayer = self.datalayer.analysisLayer
        long = analysisLayer.getDeviceTypeTransmissionFrequency(deviceType, trialName, trialSetName, normalize=True,
                                                                wideFormat=False)

        long['frequency'].hist(density=True, cumulative=True, ax=ax)
        plt.title("")
        plt.xlim(0, 1)
        plt.xlabel("Normalized frequency")
        plt.ylabel("Fraction of messages")

        if self.presentation.saveFigures:
            figname = os.path.join(self.presentation.savePath, "technical",
                                   f"{trialName}_NDIR_MessagesFrequencyDistribution_{trialSetName}.png")
            plt.savefig(figname)



    def devicesProperties(self):
        devices = []
        for trialSet in self.datalayer.setup['trialSets']:
            trialSetName = trialSet['name']
            for trial in trialSet['trials']:
                trialName = trial['name']
                date = trial['createdDate']
                for device in trial['devicesOnTrial']:
                    d = {}
                    d['type'] = device['deviceTypeName']
                    d['name'] = device['deviceItemName']
                    d['latitutde'] = device['location']['coordinates'][0]
                    d['longitute'] = device['location']['coordinates'][1]
                    d['trialSetName'] = trialSetName
                    d['trialName'] = trialName
                    d['trialCreatedDate'] = date
                    devices.append(d)

        return devices

    def plot_devices(self,trialSetName,trialName,deviceType,ax=None,plot_kwargs=None,toolkitDataSource=None,display=True):
        """
        Plot map of devices type places in a specific trial set and trial.

        Parameters
        ----------
        trialSetName : str
            Trial Set Name.
        trialName: str
            Trial Name.
        deviceType: str
            Device type name.
        plotkwargs: dict
            Parameters for matplotlib.pyplot subplot.

        Returns
        -------
            fig
            ax
        """
        tiles_tk = toolkitHome.getToolkit(toolkitHome.GIS_TILES,projectName=self.datalayer.projectName)

        devices_df = self.datalayer.trialSet[trialSetName][trialName].entitiesTable.copy()
        devices_df = devices_df[devices_df['deviceTypeName']==deviceType]
        devices_df[['ITM_Latitude', 'ITM_Longitude']] = devices_df.apply(self.datalayer._process_row, axis=1)

        minx,miny,maxx,maxy = self.datalayer.get_devices_image_coordinates(trialSetName,trialName,deviceType)

        region = dict(minx=minx, maxx=maxx, maxy=maxy, miny=miny, zoomlevel=17, inputCRS=ITM, tileServer=toolkitDataSource)
        img = tiles_tk.getImageFromCorners(**region)

        if ax is None:
            plot_kwargs = plot_kwargs or {}
            fig, ax = plt.subplots(1, 1, **plot_kwargs)
        else:
            fig = ax.figure

        plot = tiles_tk.presentation.plot(img, ax=ax, display=True)
        extent = plot.get_extent()
        x_min, x_max, y_min, y_max = extent
        ax.imshow(plot.get_array(), extent=extent, origin='lower', cmap='gray')
        d = {}
        for row in devices_df.itertuples():
            x = row.ITM_Latitude
            y = row.ITM_Longitude
            try:
                d[row.stationName] += 1
            except:
                d[row.stationName] = 1

            num_of_devices_in_station = d[row.stationName]
            delta = num_of_devices_in_station * 0.02

            ax.scatter(x, y, color='red', marker='o', s=50)  # 's' controls size
            ax.text(x, y + (y_max - y_min) * delta, f"{row.deviceItemName}", color='red', fontsize=20, ha='center',
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.8))

        if display:
            plt.show()
        else:
            plt.close(fig)

        return fig, ax

    def generate_latex_folder(self,latex_template,folder_path):
        """
        Save folder for overleaf website upload to transform to PDF.

        Parameters
        ----------
        latex_template: str
            Latex Template.
        folder_path: str
            Path to save folder

        Returns
        -------
        """
        data = {}
        data['trialSets'] = []
        os.makedirs(folder_path, exist_ok=True)
        for trialSet in self.datalayer.setup['trialSets']:
            trialSet_dict = {}
            trialSet_dict['trialSet_name'] = trialSet['name']
            trialSet_dict['trials'] = []
            for trial in trialSet['trials']:
                trial_dict = {}
                trial_dict['trial_name'] = trial['name']
                devices_df = self.datalayer.trialSet['Measurements']['Measurements'].entitiesTable
                trial_dict['devices'] = []
                for device_name in devices_df['deviceTypeName'].unique():
                    device_dict = {}
                    fig , _ = self.plot_devices(trialSetName=trialSet['name'],trialName=trial['name'],device=device_name,display=False)
                    image_path = os.path.join(folder_path,f"{device_name}.png")
                    fig.savefig(image_path)
                    device_dict['device_name'] = device_name
                    device_dict['map_image_path'] = f"{device_name}.png"
                    device_dict['locations_table'] = []
                    device_df = devices_df[devices_df['deviceTypeName'] == device_name]
                    for row in device_df.itertuples():
                        location = {}
                        location["latitude"] = row.Latitude
                        location['longitude'] = row.Longitude
                        location['device_name'] = str(row.deviceItemName).replace("_", " ")
                        location['station'] = str(row.stationName).replace("_", " ")
                        device_dict['locations_table'].append(location)
                    trial_dict['devices'].append(device_dict)
                trialSet_dict['trials'].append(trial_dict)
            data['trialSets'].append(trialSet_dict)


        template = jinja2.Template(latex_template)
        latex_content = template.render(trialSets=data["trialSets"])
        tex_path = os.path.join(folder_path,f"latex_document.tex")
        with open(tex_path, "w", encoding="utf-8") as file:
            file.write(latex_content)
        print(f"LaTeX document generated at: {tex_path}")