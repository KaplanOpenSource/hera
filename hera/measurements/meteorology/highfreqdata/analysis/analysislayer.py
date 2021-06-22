import pandas
import dask
from .turbulencestatistics import singlePointTurbulenceStatistics
from ....experiment.experiment import experimentToolKit

class RawdataAnalysis:

    _datalayer = None

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self, datalayer):
        self._datalayer = datalayer


    def singlePointTurbulenceStatistics(self,
                                        deviceNameOrData,
                                        samplingWindow,
                                        start,
                                        end,
                                        height,
                                        buildingHeight,
                                        averagedHeight,
                                        inmemory=False,
                                        isMissingData=False,
                                        **kwargs):
        """
        This method loads the raw data that corresponds to the requirements (projectName, station, instrument.. ) and
        creates a turbulence calculator with the desirable sampling window.


        Parameters
        ----------
        deviceNameOrData : str / pandas.DataFrame / dask.Dataframe / None
            the data to process.

            if str, queries the database with the deviceName as input.

             if NOne, query the database only with kwargs.

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
            Other query arguments for the database.

        Returns
        -------
        singlePointTurbulenceStatistics
            A turbulence calculator of the loaded raw data.
        """

        identifier = {'projectName': self.datalayer.projectName,
                      'samplingWindow': samplingWindow,
                      'height': height,
                      'buildingHeight':buildingHeight,
                      'averagedHeight':averagedHeight,
                      'start': start,
                      'end': end
                      }
        identifier.update(kwargs)

        if isinstance(deviceNameOrData,pandas.DataFrame) or isinstance(deviceNameOrData,dask.dataframe.DataFrame):

            rawData=deviceNameOrData

        else:
<<<<<<< HEAD
            raise ValueError("deviceNameOrData must be a dask/pandas dataframe")


        return singlePointTurbulenceStatistics(rawData = rawData, metadata=identifier, isMissingData=isMissingData)

    #
    # def singlePointStatisticsFromSonicData(self, data, samplingWindow, isMissingData=False):
    #     """
    #     This method returns turbulence calculator from a given data and sampling window.
    #
    #     Parameters
    #     ----------
    #
    #     data : pandas.DataFrame/dask.dataframe
    #         The raw data for the calculations.
    #
    #     samplingWindow : str
    #         The desirable sampling window.
    #
    #     isMissingData : bool, optional, default False
    #         A flag if there is a missing data to compute accordingly.
    #
    #     Returns
    #     -------
    #     singlePointTurbulenceStatistics
    #         A turbulence calculator of the given data.
    #     """
    #     identifier = {'samplingWindow': samplingWindow}
    #
    #     return singlePointTurbulenceStatistics(rawData=data, metadata={}, identifier=identifier, isMissingData=isMissingData)
    #
    #
    # def singlePointStatistics(self, data=None, projectName=None, **kwargs):
    #     if data is not None:
    #         return self.singlePointStatisticsFromSonicData(data=data, **kwargs)
    #     elif projectName is not None:
    #         return self.singlePointTurbulenceStatisticsFromDB(projectName=projectName, **kwargs)
    #     else:
    #         raise ValueError("'data' argument or 'projectName' argument must be delivered")
    #
    #
=======
            raise ValueError("'data' argument or 'projectName' argument must be delivered")


















def getSinglePointTurbulenceStatisticsFromDB(project,
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


    identifier = {'projectName': project.getProjectName(),
                  'samplingWindow': samplingWindow,
                  'station': None,
                  'instrument': None,
                  'height': None,
                  'start': start,
                  'end': end
                  }
    identifier.update(kwargs)

    projectData = project.getMetadata()[['height', 'instrument', 'station']].drop_duplicates()
    rawData     = project.getSonicData(stationName=stationName,
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


def getSinglePointStgetStisticsFromSonicData(self, data, samplingWindow, isMissingData=False):
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


def getSinglePointStatistics(self, data=None, project=None, **kwargs):
    if data is not None:
        return self.singlePointStatisticsFromSonicData(data=data, **kwargs)
    elif project is not None:
        return self.singlePointTurbulenceStatisticsFromDB(project, **kwargs)
    else:
        raise ValueError("'data' argument or 'projectName' argument must be delivered")





>>>>>>> c8d9ea77f73852c66b87efa24554ffa2b9b9e543
