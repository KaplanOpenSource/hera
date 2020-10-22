import matplotlib.pyplot as plt
import numpy
import seaborn
import dask.dataframe
from itertools import product

from .analysis import calcHourlyDist
from .analysis import seasonsdict

class Plots(object):
    """
    This class
    """

    _contourvalsdict = None
    _plotfieldaxfuncdict =None
    _axfuncdict =None
    _labelsdict=None
    _scatterdict=None

    def __init__(self):

        self._contourvalsdict=dict(under_value=0.1,
                                   contourskip=2,
                                   contourfnum=10,
                                   max_value=1.0
                                   )

        self._plotfieldaxfuncdict=dict(WD=dict(#set_ylim=[0, 360],
                                              #set_yticks=[x for x in range(0, 361, 30)],
                                              set_ylabel= 'Wind Direction [° from N]'
                                              ),
                                      WS=dict(#set_ylim=[0,20],
                                          set_ylabel='Wind Speed [m/s]',
                                          #set_yticks= [x for x in range(0, 21, 4)],
                                          #set_yticklabels= [str(x) if x % 2 == 0 else "" for x in range(0, 21,4)]
                                          ),
                                      RH=dict(#set_ylim=[0,100],
                                          #set_yticks= [x for x in range(0, 101,5)],
                                          #set_yticklabels= [str(x) if x % 10 == 0 else "" for x in range(0, 101,5)],
                                          set_ylabel= 'Relative Humidity [%]')
                                      )

        self._axfuncdict=dict(set_xlim= [0, 24],
                              set_xticks=[x for x in range(0,25)],
                              set_xticklabels= [str(x) if x % 2 == 0 else "" for x in range(0, 25)],
                              set_xlabel= 'Time [Hours]'
                              )

        self._labelsdict=dict(levels=numpy.round(numpy.linspace(self._contourvalsdict['under_value'],
                                                                self._contourvalsdict['max_value'],
                                                                self._contourvalsdict['contourfnum']), 2)[::self._contourvalsdict['contourskip']],
                              inline=True,
                              fontsize=8,
                              fmt='%1.2f'
                              )

        self._scatterdict=dict(zorder=1,
                               color='k',
                               marker='.',
                               size=0.5,
                               edgecolors='k',
                               alpha=0.55,
                               legend=False
                               )

    def _getCountourDict(self, params):

        """
        Make a dictionary with the properties for contour plot
        :param params: parametes from the user
        :return: dict
        """

        return dict(zorder=3,
                    levels=numpy.round(numpy.linspace(params['under_value'],
                                                      params['max_value'],
                                                      params['contourfnum']), 2)[::self._contourvalsdict['contourskip']],
                    linewidths=0.5,
                    colors='k')

    def _getContourfDict(self, params,Cmapname='jet'):

        """
        Make a dictionary with the properties for contourf plot

        :param params: parametes from the user
        :return: dict
        """

        return  dict(zorder=2,
                     cmap=self._getcmap(name=Cmapname,under=True,undercolor='0.9',over=False,overcolor=None,alpha=0.05),
                     levels=numpy.round(numpy.linspace(params['under_value'],
                                                       params['max_value'],
                                                       params['contourfnum']), 2),
                     extend='min'
                     )

    def _getcmap(self,under=False,undercolor='0.9',over=False,overcolor='0.9',alpha=0.05,name='jet'):

        """
        Creates a colormap object with under/over range properties

        Parameters
        ----------

        name : string.
        The name of requested colormap (for ex: 'jet')
        under : boolean, default False.
         Whether or not to set a low out-of-range color
        undercolor : String.
         The requested under color name
        over : boolean, default False.
         Whether or not to set a high out-of-range color
        overcolor : String.
         The requested over color name
        alpha : float, default 0.05.
         Out-of-range color transparency

        Returns
        -------
        cmap : colormap object
        """

        cmap=plt.get_cmap(name)
        if under is not False:
            cmap.set_under(color=undercolor,alpha=alpha)
        if over is not False:
            cmap.set_over(color=overcolor, alpha=alpha)

        return cmap

class SeasonalPlots(Plots):


    def __init__(self):

        super().__init__()

    def plotProbContourf_bySeason(self, data, plotField, levels = None, scatter = True, withLabels = True, colorbar=True,
                                  Cmapname='jet',ax=None,scatter_properties = dict(),contour_values = dict(),
                                  contour_properties = dict(),contourf_properties = dict(), labels_properties = dict(),
                                  ax_functions_properties = dict(),normalization = 'max_normalized', figsize=[15, 10]):

        """
        applying plotProbContourf by season. see plotProbContourf for full documentation

        Parameters
        ----------

        data : dask or pandas dataframe
        plotField :
         The data column name to plot
        levels : list, optional.
        scatter : boolean. default is True
        withLabels : boolean. default is True
        colorbar : boolean. default is True
        ax : optional
        scatter_properties : dict, optional
        contour_values : dict, optional
        contour_properties : dict, optional
        contourf_properties : dict, optional
        labels_properties : dict, optional
        ax_functions_properties : dict, optional
            A dict with axes functions to add/replace the default functions applied on the ax
        normalization:string, default 'max_normalized'
        figsize : list, default is [15,10]

        returns
        -------

        ax :

        """

        if ax is None:
            fig, ax = plt.subplots(2,2,figsize=figsize)
        else:
            plt.sca(ax)

        curdata=data
        curdata = curdata.assign(curdate=curdata.index)
        curdata = curdata.assign(monthonly=curdata.curdate.dt.month)

        axPositionList = [x for x in product(range(ax.shape[0]), range(ax.shape[1]))]

        for axPosition, season in zip(axPositionList, seasonsdict):
            qstring = 'monthonly in %s' % seasonsdict.get(season)['monthes']
            seasondata=curdata.query(qstring)

            CS,CFS,ax_i=dailyplots.plotProbContourf(self,
                                                    seasondata,
                                                    plotField,
                                                    ax=ax[axPosition[0], axPosition[1]],
                                                    colorbar=False,
                                                    Cmapname=Cmapname,
                                                    levels = levels,
                                                    scatter = scatter,
                                                    withLabels = withLabels,
                                                    scatter_properties = scatter_properties,
                                                    contour_values = contour_values,
                                                    contour_properties = contour_properties,
                                                    contourf_properties = contourf_properties,
                                                    labels_properties = labels_properties,
                                                    ax_functions_properties= ax_functions_properties,
                                                    normalization = normalization)

            ax_i.set_title('%s %s' % (season,seasonsdict.get(season)['strmonthes']))

        plt.tight_layout(rect=[0, 0, 1, 0.95])
        if colorbar==True:
            plt.colorbar(ax=ax,ticks=CFS.levels)

        return ax

class DailyPlots(Plots):

    _linedict=None

    def __init__(self):

        super().__init__()

        self._linedict=dict(zorder=4,
                            color='magenta',
                            linestyle='-',
                            linewidth=3
                            # size=0.5,
                            # edgecolors='k',
                            # alpha=0.55,
                            # legend=False
                            )

    def plotScatter(self,data,plotField,ax=None,scatter_properties=dict(),ax_functions_properties=dict()):

        """

        Make a scatter plot for selected plotField vs daily time (in 24 hours)

        Parameters
        ----------
        data : dask or pandas dataframe
        plotField : string.
            The data column name to plot
        ax : Axes object, optional
            Axes object to plot in
        scatter_properties : dict, optional
            dict with parameters to add/update the default dict passed to seaborn.scatterplot
        ax_functions_properties : dict, optional
            A dict with axes functions to add/replace the default functions applied on the ax

        returns
        -------
        ax : The axes object
        """

        if ax is None:
            fig, ax = plt.subplots()
        else:

            plt.sca(ax)

        scatter_props = dict(self._scatterdict)
        scatter_props.update(scatter_properties)

        # curdata = data.dropna(subset=[plotField])
        curdata = data.copy()
        curdata[plotField] = curdata[plotField].where(curdata[plotField] > -9000)

        # curdata = curdata.query("%s > -9990" % plotField)
        curdata = curdata.assign(curdate=curdata.index)
        curdata = curdata.assign(houronly=curdata.curdate.dt.hour + curdata.curdate.dt.minute / 60.)


        ax= seaborn.scatterplot(curdata['houronly'], curdata[plotField], ax=ax, **scatter_props)

        ax_func_props = dict(self._plotfieldaxfuncdict.get(plotField, dict()))
        ax_func_props.update(self._axfuncdict)
        ax_func_props.update(ax_functions_properties)

        for func in ax_func_props:
            getattr(ax, func)(ax_func_props[func])

        return ax

    def dateLinePlot(self,data,plotField,date,legend=True,ax=None,line_properties=dict(),ax_functions_properties=dict()):

        """
        Make a line plot for selected plotField vs daily time (in 24 hours)

        Parameters
        ----------
        data : dask or pandas dataframe.
            The data to plot
        plotField : string
            The data column name to plot
        date : string
            The date to plot (in format of 'YYYY-MM-DD').
        legend : boolean
            default True
        ax : Axes object, optional
            Axes object to plot in
        line_properties : dict' optional

        ax_functions_properties : dict, optional
            A dict with axes functions to add/replace the default functions applied on the ax

        returns
        -------

        ax :
        line :

        """


        if ax is None:
            fig, ax = plt.subplots()
        else:
            plt.sca(ax)

        line_props = dict(self._linedict)
        line_props.update(line_properties)

        curdata = data.copy()
        curdata[plotField] = curdata[plotField].where(curdata[plotField] > -9000)

        # curdata = data.dropna(subset=[plotField])
        # curdata = curdata.query("%s > -9990" % plotField)
        curdata = curdata.assign(curdate=curdata.index)
        curdata = curdata.assign(dateonly=curdata.curdate.dt.date.astype(str))
        curdata = curdata.assign(houronly=curdata.curdate.dt.hour + curdata.curdate.dt.minute / 60.)

        qstring = "dateonly == '%s'" % date
        if isinstance(curdata, dask.dataframe.core.DataFrame):
            dailydata = curdata.query(qstring).compute()
        else:
            dailydata = curdata.query(qstring)

        # ax= seaborn.lineplot(dailydata['houronly'], dailydata[plotField], ax=ax, **line_props)
        line = plt.plot(dailydata['houronly'], dailydata[plotField], axes=ax, label=date, **line_props)
        if legend==True:
            ax.legend()

        ax_func_props = dict(self._plotfieldaxfuncdict.get(plotField, dict()))
        ax_func_props.update(self._axfuncdict)
        ax_func_props.update(ax_functions_properties)

        for func in ax_func_props:
            getattr(ax, func)(ax_func_props[func])

        return ax, line

    def plotProbContourf(self, data, plotField, levels=None, scatter=True, withLabels=True, colorbar=True,Cmapname='jet', ax=None, scatter_properties=dict(),
                         contour_values=dict(), contour_properties=dict(), contourf_properties=dict(), labels_properties=dict(),
                         ax_functions_properties=dict(), normalization='max_normalized'):

        """
        Make a probability contour plot for selected plotField vs daily time (in 24 hours)

        Parameters
        ----------

        data : dask or pandas dataframe.
            The data to plot
        plotField : string.
            The data column name to plot
        levels : list, optional.
            Overrides the defaults values of countour and countourf levels
        scatter : boolean. default is True
            Whether or not to add a scatterplot of the data
        withLabels : boolean. default is True
            Whether or not to add a labels to the contour lines
        colorbar : boolean. default is True
            Whether or not to add a colorbar
        Cmapname : string, default 'jet'
            The name of requested colormap
        ax : Axes object, optional
            Axes object to plot in
        scatter_properties : dict, optional
            dict with parameters to add/update the default dict passed to seaborn.scatterplot
        contour_values : dict, optional
            A dict with parameters to add/update the default contour and contourf *values*
            passed to _getCountourDict and _getCountourfDict
        contour_properties : dict, optional
            A dict with contour parameters to add/replace the default contour parmeters passed to plt.contour
        contourf_properties : dict, optional
            A dict with contourf parameters to add/replace the default contour parmeters passed to plt.contourf
        labels_properties : dict, optional
            A dict with labels parameters to add/replace the default labels parmeters passed to ax.clabel
        ax_functions_properties : dict, optional
            A dict with axes functions to add/replace the default functions applied on the ax
        normalization : string, default 'max_normalized'
            The method of data normalization. see hera.analytics.statistics.calcDist2d for more details

        Returns
        ------

        CS : contour set
        CFS : contourf set
        ax : The axes object
        """

        # Compute histogram #

        x_hist, y_hist, M_hist = calcHourlyDist(data, plotField, normalization=normalization)


        # Read and update contour and contourf properties #

        conrourvals_props=dict(self._contourvalsdict)
        conrourvals_props.update(contour_values)

        contour_props = self._getCountourDict(conrourvals_props)
        contourf_props = self._getContourfDict(conrourvals_props,Cmapname)


        if normalization=='y_normalized':
            contour_props.update(dict(levels= numpy.linspace(conrourvals_props['under_value'],
                                                             numpy.ceil(M_hist.max() * 100) / 100,
                                                             conrourvals_props['contourfnum'])[::conrourvals_props['contourskip']]))

            contourf_props.update(dict(levels=numpy.round(numpy.linspace(conrourvals_props['under_value'],
                                                                         numpy.ceil(M_hist.max() * 100) / 100,
                                                                         conrourvals_props['contourfnum']), 2)))

        contour_props.update(contour_properties)
        contourf_props.update(contourf_properties)

        # Read and update axes properties #

        ax_func_props = dict(self._plotfieldaxfuncdict.get(plotField, dict()))
        ax_func_props.update(self._axfuncdict)
        ax_func_props.update(ax_functions_properties)


        countour_levels = contour_props.pop("levels")
        countourf_levels = contourf_props.pop("levels")

        # Read and update levels and labels properties #

        labels_props = dict(self._labelsdict)
        if normalization == 'y_normalized':
            labels_props.update(dict(levels= numpy.linspace(conrourvals_props['under_value'],
                                                            numpy.ceil(M_hist.max() * 100) / 100,
                                                            conrourvals_props['contourfnum'])[::conrourvals_props['contourskip']]))

        if levels is not None:
            countour_levels = levels
            countourf_levels = levels
            del labels_props['levels']

        # Plot scatter, conrour and contourf #

        if ax is None:
            fig, ax = plt.subplots()
        else:
            plt.sca(ax)

        if scatter==True:
            ax = dailyplots.plotScatter(self,data,plotField,ax=ax,scatter_properties=scatter_properties)
            # ax = self.plotScatter(data,plotField,ax=ax,scatter_properties=scatter_properties)

        CS = plt.contour(x_hist, y_hist, M_hist,levels=countour_levels,**contour_props)

        if withLabels==True:

            labels_props.update(labels_properties)
            labels_list = labels_props.pop("levels",CS.levels)

            ax.clabel(CS, labels_list,**labels_props)

        CFS= plt.contourf(x_hist, y_hist, M_hist,levels=countourf_levels ,**contourf_props)

        # Apply ax functions from ax_func_props

        for func in ax_func_props:
            getattr(ax, func)(ax_func_props[func])

        if colorbar==True:
            plt.colorbar(ax=ax, ticks=countourf_levels)

        return CS,CFS,ax



seasonplots = SeasonalPlots()
dailyplots  = DailyPlots()

