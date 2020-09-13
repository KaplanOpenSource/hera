import pandas
from .singlepointturbulencestatistics import singlePointTurbulenceStatistics

class RawdataAnalysis:

    _project = None

    @property
    def project(self):
        return self._project

    def __init__(self,project):
        self._project = project



    def singlePointTurbulenceStatisticsFromDB(self,
                                              samplingWindow,
                                              start,
                                              end,
                                              usePandas=False,
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

        usePandas : bool, positional, default False
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
