import pandas
from .turbulencestatistics import singlePointTurbulenceStatistics
from .abstractcalculator import AbstractCalculator

class RawdataAnalysis:

    _project = None

    @property
    def project(self):
        return self._project

    def __init__(self,project):
        self._project = project

    # def TrhMeanDataFromDB(self,
    #                       projectName,
    #                       samplingWindow,
    #                       stationName,
    #                       height,
    #                       start,
    #                       end,
    #                       inMemory=None,
    #                       compute_mode='not_from_db_and_not_save',
    #                    saveProperties=None,
    #                    usePandas=False, **kwargs):
    #     """
    #     This method loads the raw data that corresponds to the requirements (projectName, station, instrument..., with the
    #     assumption that instrument = "TRH", "Tct_TRH" etc.), calculates the mean temperature with the desirable sampling
    #     window and returns a pandas dataframe of the result. It also updates the cache depending on the save mode.
    #
    #
    #     Parameters
    #     ----------
    #     projectName : str
    #         The name of the project.
    #
    #     samplingWindow : str
    #         The desirable sampling window.
    #
    #     start : str/pandas.Timestamp
    #         Datetime of the begin.
    #
    #     end : str/pandas.Timestamp
    #         Datetime of the end.
    #
    #     compute_mode : str, positional, default 'not_from_db_and_not_save'. Other options are 'from_db_and_not_save',
    #                     'not_from_db_and_save', or 'from_db_and_save'
    #         Instruction of how to compute and whether or not to save the data into the cache.
    #
    #     saveProperties : dict, positional, default None
    #         Saving properties for the cache, goes into the abstract calculator.
    #
    #     usePandas : bool, positional, default False
    #         A flag of whether or not to use pandas.
    #
    #     inMemory : boolean
    #         Default value is None.
    #
    #     kwargs :
    #         Other query arguments.
    #
    #     Returns
    #     -------
    #     A pandas dataframe of the resulting averaged data
    #     """
    #
    #     if type(start) is str:
    #         start = pandas.Timestamp(start)
    #
    #     if type(end) is str:
    #         end = pandas.Timestamp(end)
    #
    #     identifier = {'projectName': self.project.getProjectName(),
    #                   'samplingWindow': samplingWindow,
    #                   'station': None,
    #                   'instrument': None,
    #                   'height': None,
    #                   'start': start,
    #                   'end': end
    #                   }
    #     identifier.update(kwargs)
    #
    #     projectData = self.project.getMetadata()[['height', 'instrument', 'station']].drop_duplicates()
    #
    #     if identifier['station'] is not None:
    #         stationData = projectData.query("station=='%s'" % identifier['station']).iloc[0]
    #         identifier['buildingHeight'] = stationData.get('buildingHeight', None)
    #         identifier['averagedHeight'] = stationData.get('averagedHeight', None)
    #
    #     calculator = AbstractCalculator(rawData=rawData, metadata=projectData, identifier=identifier)
    #
    #     if saveProperties is not None:
    #         calculator.set_saveProperties(**saveProperties)
    #
    #     if calculator._InMemoryAvgRef is None:
    #         calculator._InMemoryAvgRef = inMemory
    #
    #     calculator._TemporaryData = rawData.resample(samplingWindow).mean()
    #     calculator._CalculatedParams += [[col, {}] for col in rawData.columns]
    #
    #     return calculator.compute(mode=compute_mode)

    def singlePointTurbulenceStatisticsFromDB(self,
                                              stationName,
                                              height,
                                              samplingWindow,
                                              start,
                                              end,
                                              inmemory=False,
                                              isMissingData=False, **kwargs):
        """
        This method loads the raw data that corresponds to the requirements (projectName, station, instrument.. ) and
        creates a turbulence calculator with the desirable sampling window.


        Parameters
        ----------
        projectName : str
            The name of the project.

        samplingWindow : str
            The desirable sampling window.

        start : str/pandas.Timestamp
            Datetime of the begin.

        end : str/pandas.Timestamp
            Datetime of the end.

        inmemory : bool, positional, default False
            A flag of whether or not to use pandas.

        isMissingData : bool, positional, default False
            A flag if there is a missing data to compute accordingly.

        kwargs :
            Other query arguments.

        Returns
        -------
        singlePointTurbulenceStatistics
            A turbulence calculator of the loaded raw data.
        """


        identifier = {'projectName': self.project.getProjectName(),
                      'samplingWindow': samplingWindow,
                      'station': None,
                      'instrument': None,
                      'height': None,
                      'start': start,
                      'end': end
                      }
        identifier.update(kwargs)

        projectData = self.project.getMetadata()[['height', 'instrument', 'station']].drop_duplicates()
        rawData     = self.project.getSonicData(stationName=stationName,
                                                height=height,
                                                inmemory=inmemory,
                                                start=start,
                                                end=end,
                                                **kwargs)


        if identifier['station'] is not None:
            stationData = projectData.query("station=='%s'" % identifier['station']).iloc[0]
            identifier['buildingHeight'] = stationData.get('buildingHeight', None)
            identifier['averagedHeight'] = stationData.get('averagedHeight', None)

        return singlePointTurbulenceStatistics(rawData = rawData, metadata=projectData, identifier=identifier, isMissingData=isMissingData)


    def singlePointStatisticsFromSonicData(self, data, samplingWindow, isMissingData=False):
        """
        This method returns turbulence calculator from a given data and sampling window.

        Parameters
        ----------

        data : pandas.DataFrame/dask.dataframe
            The raw data for the calculations.

        samplingWindow : str
            The desirable sampling window.

        isMissingData : bool, optional, default False
            A flag if there is a missing data to compute accordingly.

        Returns
        -------
        singlePointTurbulenceStatistics
            A turbulence calculator of the given data.
        """
        identifier = {'samplingWindow': samplingWindow}

        return singlePointTurbulenceStatistics(rawData=data, metadata={}, identifier=identifier, isMissingData=isMissingData)


    def singlePointStatistics(self, data=None, projectName=None, **kwargs):
        if data is not None:
            return self.singlePointStatisticsFromSonicData(data=data, **kwargs)
        elif projectName is not None:
            return self.singlePointTurbulenceStatisticsFromDB(projectName=projectName, **kwargs)
        else:
            raise ValueError("'data' argument or 'projectName' argument must be delivered")


