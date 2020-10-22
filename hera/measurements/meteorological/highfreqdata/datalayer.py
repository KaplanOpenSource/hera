import pandas
import dask

from ....datalayer import ProjectMultiDBPublic
from .analysis.analysislayer import RawdataAnalysis


class RawData(ProjectMultiDBPublic):

    _dataType = None

    @property
    def dataType(self):
        return self._dataType


    def __init__(self, projectName,dataType,databaseNameList=None,useAll=True): # databaseNameList=None, useAll=False

        super().__init__(projectName=projectName,
                         publicProjectName=f"{dataType}_highfreq",
                         databaseNameList=databaseNameList,
                         useAll=useAll)


    def _getRawData(self,stationName,height=None,start=None,end=None,inmemory=False,**kwargs):
        """
            Returns raw data of the requested type.

        Parameters
        ----------

        stationName: str
                The name of the station

        dataType: str
                The type of data to retrieve.
                Either 'rawsonic' or 'temperature'.

        height: float
                The height of the sonic from the mast base.

        start:  pandas.Timestamp or str
                Start time for date range.

        end:  pandas.Timestamp or str
                End time for date range.

        Returns
        -------

        dask.dataframe or pandas.Dataframe.
            The raw data.
        """
        if type(start) is str:
            start = pandas.Timestamp(start)

        if type(end) is str:
            end = pandas.Timestamp(end)

        # shoudl add a document type='SonicData' or something.
        qry = dict(type=self.dataType,
                   stationName=stationName,
                   height=height)

        kwargs.update(qry)


        docList = self.getMeasurementsDocuments(**kwargs)
        rawData = dask.dataframe.concat([doc.getData() for doc in docList])[start:end]

        return rawData.compute() if inmemory else rawData


class RawSonic(RawData):


    _analysis = None

    @property
    def analysis(self):
        return self._analysis


    def __init__(self, projectName, databaseNameList=None, useAll=False):

        super().__init__(projectName=projectName,
                         publicProjectName="RawSonic",
                         databaseNameList= databaseNameList,
                         useAll = useAll)

        self._analysis = RawdataAnalysis(self)



    def getData(self,stationName,height=None,start=None,end=None,inmemory=False,**kwargs):
        """
            Return the sonic raw data either in dask or loaded to memeory (pandas).

        Parameters
        ----------

        stationName: str
                The name of the station
        height: float
                The height of the sonic from the mast base.

        start:  pandas.Timestamp or str
                Start time for date range.

        end:  pandas.Timestamp or str
                End time for date range.

        kwargs : dict,
            Additional query

        Returns
        -------

        dask.dataframe or pandas.Dataframe.
            The raw data.

        """

        self._getRawData(stationName=stationName,
                         dataType='RawSonic',
                         height=height,
                         start=start,
                         end=end,
                         inmemory=inmemory, **kwargs)




class Temperature(RawData):



    def __init__(self, projectName, databaseNameList=None, useAll=False):

        super().__init__(projectName=projectName,
                         publicProjectName="RawSonic",
                         databaseNameList= databaseNameList,
                         useAll = useAll)


    def getData(self,stationName,height=None,start=None,end=None,inmemory=False,**kwargs):
            """
                Return the temperature raw data either in dask or loaded to memeory (pandas).


                Parameters
                ----------

                stationName: str
                        The name of the station
                height: float
                        The height of the sonic from the mast base.

                start:  pandas.Timestamp or str
                        Start time for date range.

                end:  pandas.Timestamp or str
                        End time for date range.

                kwargs : dict,
                    Additional query

            """
            self._getRawData(stationName=stationName,
                             dataType='Temperature',
                             height=height,
                             start=start,
                             end=end,
                             inmemory=inmemory, **kwargs)


