import pandas
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import numpy as np


class PlotFields(object):

    def __init__(self, data = None):
        self._Data = data

    def plotField(self, fieldName, beginDate, endDate, resampleList=[],movingList=[], yLabel=None, ax=None, legend=False, plotRaw=False, save=False, saveTo=None):
        """
        Plot a single field in the requested date range. Adds averaged fields according to user input.

        :param fieldName:
        :param beginDate:
        :param endDate:
        :param averageList: A list of average periods. For example ['5min','10min'] adds five and ten minutes averages.
        :param yLabel:
        :param axe:
        :return:
        """
        if ax is None:
            plt.figure(figsize=[20, 10])
        else:
            plt.sca(ax)

        dataToPlot = self._Data[fieldName].loc[beginDate:endDate]
        if plotRaw:
            dataToPlot.plot(label='raw')

        for averageI in resampleList:
            dataToPlot.resample(averageI).mean().plot(label='resampled %s' % averageI)

        for averageI in movingList:
            dataToPlot.rolling(averageI).mean().plot(label='moving mean %s' % averageI)

        if yLabel is not None:
            plt.ylabel(yLabel)

        plt.xlim(dataToPlot.index.min(),dataToPlot.index.max())

        if legend:
            plt.legend()

        if save:
            plt.savefig(saveTo)

    def plotFields(self, fieldNames, beginDate, endDate, resampleList=[], movingList=[], yLabels=None, axes=None, legends=None, plotRaw=False, save=False, saveTo=None):
        """

        :param fieldNames: A list of the field names to be plotted.
        :param beginDate:
        :param endDate:
        :param averageList: A list of average periods. For example ['5min','10min'] adds five and ten minutes averages.
        :param yLabels: A list of the y labels for the plots.
        :return:
        """
        plotsNum = len(fieldNames)

        fig, axes = plt.subplots(1, plotsNum, figsize=[20*plotsNum, 10])

        for i in range(plotsNum):
            self.plotField(fieldName=fieldNames[i],
                           beginDate=beginDate,
                           endDate=endDate,
                           resampleList=resampleList,
                           movingList=movingList,
                           ax=axes[i],
                           yLabel=yLabels[i] if yLabels is not None else None,
                           legend=legends[i] if legends is not None else None,
                           plotRaw=plotRaw)

        if save:
            fig.savefig(saveTo)

    def test(self, fieldNames, beginDate, endDate, resampleList=[], movingList=[], yLabels=None, axes=None, legends=None, plotRaw=False):
        print('fieldNames: %s' % fieldNames)
        print('beginDate: %s' % beginDate)
        print('endDate: %s' % endDate)
        print('resampleList: %s' % resampleList)
        print('movingList: %s' % movingList)
        print('yLabels: %s' % yLabels)
        print('axes: %s' % axes)
        print('legends: %s' % legends)

    def plot_frequency_bar(self, data, category, grouped_by = "hour", axis = None, ax_props = None, cat_list = None, **kwargs):
        """

            plot_hist2d_windVec2(data=..,y="WD",month=3,hour=1)
                        query field is month,hour

                        querystr = ["month==3","hour==1"] ==> "month==3 and hour==1"

            plot_hist2d_windVec2(data=..,y="WD",month=3)
                        query field is month
                        querystr = ["month==3"]


            :param data:
            :param category:
            :param grouped_by:
            :param axis:
            :param ax_props: dict with axis functions and parameters. { funcname : {parameters}}
            :param cat_list:

            :param kwargs:
                    all the rest are used as query (use the and predicate on all the fields).
                    month=1
                    day = 2


            :return:
            """

        # kwargs_used_fields = ["axis", "fig", "bins", "levels"]
        # queryfield = [x for x in kwargs.key() if x not in kwargs_used_fields]
        # querystr = " and ".join(["{fieldname}=={fieldvalue}".format(fieldname=fieldname, fieldvalue=kwargs[fieldname]) for fieldname in queryfield])
        #
        # month_data = data.assign(year=data.index.year) \
        #     .assign(month=data.index.month) \
        #     .assign(day=data.index.day) \
        #     .assign(hour=data.index.hour) \
        #     .query(querystr)

        if axis is None:
            fig, axis = plt.subplots()
        else:
            plt.sca(axis)

        data_grouped = data.groupby([grouped_by, category])[category].count().unstack(category).fillna(0)
        if not cat_list is None:
            for cat in cat_list:
                if cat not in data_grouped.columns:
                    data_grouped[cat] = 0
            data_grouped = data_grouped[cat_list]

        data_grouped = data_grouped.div(data_grouped.sum(axis=1), axis=0)

        data_grouped.plot(kind='bar', stacked=True, ax = axis)

        if not ax_props is None:
            for func in ax_props:
                getattr(axis, func)(**ax_props[func])

    def plot_scatter_func(self, x, y, data = None, func_dict = None, axis = None, ax_props = None, SelfCorrelation = None, **kwargs):
        """

            plot_hist2d_windVec2(data=..,y="WD",month=3,hour=1)
                        query field is month,hour

                        querystr = ["month==3","hour==1"] ==> "month==3 and hour==1"

            plot_hist2d_windVec2(data=..,y="WD",month=3)
                        query field is month
                        querystr = ["month==3"]

            :param x:
            :param y:
            :param data:
            If data is None (default) then y,x are interpreted as vector data; else interpreted as names of columns in data.

            :param func_dict: dictionary of functions to plot against the data. Must be of the format {..., index:(lambda function, legend of function), ...}
            :param axis:
            :param ax_props: dict with axis functions and parameters. { funcname : {parameters}}
            :param SelfCorrelation: whether or not to plot the self correlation test of the data (scatter plot of data in reverse order of y values)

            :param kwargs:
                    all the rest are used as query (use the and predicate on all the fields).
                    month=1
                    day = 2


            :return:
            """

        # kwargs_used_fields = ["axis", "fig", "bins", "levels"]
        # queryfield = [x for x in kwargs.key() if x not in kwargs_used_fields]
        # querystr = " and ".join(["{fieldname}=={fieldvalue}".format(fieldname=fieldname, fieldvalue=kwargs[fieldname]) for fieldname in queryfield])
        #
        # month_data = data.assign(year=data.index.year) \
        #     .assign(month=data.index.month) \
        #     .assign(day=data.index.day) \
        #     .assign(hour=data.index.hour) \
        #     .query(querystr)

        if axis is None:
            fig, axis = plt.subplots()
        else:
            plt.sca(axis)

        if data is None:
            xdata = x; ydata = y
        else:
            xdata = data[x]; ydata = data[y]

        sns.scatterplot(x=xdata, y=ydata, ax = axis, zorder=2)
        if SelfCorrelation:
            sns.scatterplot(x=np.array(xdata), y=np.array(ydata)[::-1], ax = axis, zorder=1)

        if not func_dict is None:
            # z = np.arange(min(xdata) - 1, max(xdata) + 1, 0.01)
            # for fun in func_dict:
            #     axis.plot(z, np.array([func_dict[fun][0](zi) for zi in z]),
            #              label=func_dict[fun][1])
            z1 = np.append(np.arange(min(xdata) - 1, 0, 0.01), -0.00001)
            z2 = np.arange(0.00001, max(xdata) + 1, 0.01)
            for fun in func_dict:
                lines = axis.plot(z1, np.array([func_dict[fun][0](zi) for zi in z1]),
                         label=func_dict[fun][1])
                axis.plot(z2, np.array([func_dict[fun][0](zi) for zi in z2]), color = lines[0].get_color())

        if not ax_props is None:
            for func in ax_props:
                getattr(axis, func)(**ax_props[func])


if __name__ == "__main__":
    dirPath = '/mnt/public/Lidor/HaifaData/hdf'
    station = 'Check_Post'
    type = 'Sonic'
    year = '2014'
    month = '09'
    height = 'h9m'
    day = '18'
    matplotlib.style.use('fivethirtyeight')

    fullPath = '%s/%s/%s/%s/%s_%s/TOA5_Raw_%s_%s_%s_%s.hdf' % (dirPath, station, type, height, year, month, type, year, month, day)

    data = pandas.read_hdf(fullPath, key='data')

    startTime = data.index[0]
    endTime = data.index[10000]
    varsToPlot = ['u', 'w']

    a = PlotFields(data=data)
    a.plotFields(fieldNames=varsToPlot, beginDate=startTime, endDate=endTime, resampleList=['1min', '2min'], yLabels=varsToPlot)

    plt.show()

