import pandas
import dask.dataframe
from ..analytics.turbulencecalculator import TurbulenceCalculator
from ..analytics.abstractcalculator import AbstractCalculator, getSaveData, SaveDataHandler
from .... import datalayer



def getTrhMeanData(projectName, samplingWindow, start, end, compute_mode='not_from_db_and_not_save', saveProperties=None,
                       usePandas=False, inMemory = None, **kwargs):
    """
    This method loads the raw data that corresponds to the requirements (projectName, station, instrument..., with the
    assumption that instrument = "TRH", "Tct_TRH" etc.), calculates the mean temperature with the desirable sampling
    window and returns a pandas dataframe of the result. It also updates the cache depending on the save mode.


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

    compute_mode : str, positional, default 'not_from_db_and_not_save'. Other options are 'from_db_and_not_save',
                    'not_from_db_and_save', or 'from_db_and_save'
        Instruction of how to compute and whether or not to save the data into the cache.

    saveProperties : dict, positional, default None
        Saving properties for the cache, goes into the abstract calculator.

    usePandas : bool, positional, default False
        A flag of whether or not to use pandas.

    inMemory : boolean
        Default value is None.

    kwargs :
        Other query arguments.

    Returns
    -------
    A pandas dataframe of the resulting averaged data
    """

    if type(start) is str:
        start = pandas.Timestamp(start)

    if type(end) is str:
        end = pandas.Timestamp(end)


    docList = datalayer.Measurements.getDocuments(projectName = projectName, **kwargs)
    dataList = [doc.getData(usePandas=usePandas) for doc in docList]

    rawData = pandas.concat(dataList) if usePandas else dask.dataframe.concat(dataList)
    rawData = rawData[start:end]

    identifier = {'projectName': projectName,
                  'samplingWindow': samplingWindow,
                  'station': None,
                  'instrument': None,
                  'height': None,
                  'start': start,
                  'end': end
                  }
    identifier.update(kwargs)

    projectData = datalayer.Project(projectName=projectName).getMetadata()[['height', 'instrument', 'station']].drop_duplicates()

    if identifier['station'] is not None:
        stationData = projectData.query("station=='%s'" % identifier['station']).iloc[0]
        identifier['buildingHeight'] = stationData.get('buildingHeight', None)
        identifier['averagedHeight'] = stationData.get('averagedHeight', None)

    calculator = AbstractCalculator(rawData=rawData, metadata=projectData, identifier=identifier)

    if saveProperties is not None:
        calculator.set_saveProperties(**saveProperties)

    if calculator._InMemoryAvgRef is None:
        calculator._InMemoryAvgRef = inMemory

    calculator._TemporaryData = rawData.resample(samplingWindow).mean()
    calculator._CalculatedParams +=  [[col, {}] for col in rawData.columns]

    return calculator.compute(mode = compute_mode)


def getTrhMeanData_old(projectName, samplingWindow, start, end, from_db = False, saveProperties=None, usePandas=False, **kwargs):
    """
    This method loads the raw data that corresponds to the requirements (projectName, station, instrument..., with the
    assumption that instrument = "TRH", "Tct_TRH" etc.), calculates the mean temperature with the desirable sampling 
    window and returns a pandas dataframe of the result. It also updates the cache depending on the save mode.


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

    save : bool, positional, default True
        A flag of whether or not to save the results into the cache.

    usePandas : bool, positional, default False
        A flag of whether or not to use pandas.

    isMissingData : bool, positional, default False
        A flag if there is a missing data to compute accordingly.

    kwargs :
        Other query arguments.

    Returns
    -------
    TurbulenceCalculator
        A turbulence calculator of the loaded raw data.
    """

    if type(start) is str:
        start = pandas.Timestamp(start)

    if type(end) is str:
        end = pandas.Timestamp(end)


    docList = datalayer.Measurements.getDocuments(projectName = projectName, **kwargs)
    dataList = [doc.getData(usePandas=usePandas) for doc in docList]

    rawData = pandas.concat(dataList) if usePandas else dask.dataframe.concat(dataList)
    rawData = rawData[start:end]

    identifier = {'projectName': projectName,
                  'samplingWindow': samplingWindow,
                  'station': None,
                  'instrument': None,
                  'height': None,
                  'start': start,
                  'end': end
                  }
    identifier.update(kwargs)

    projectData = datalayer.Project(projectName=projectName).getMetadata()[['height', 'instrument', 'station']].drop_duplicates()

    if identifier['station'] is not None:
        stationData = projectData.query("station=='%s'" % identifier['station']).iloc[0]
        identifier['buildingHeight'] = stationData.get('buildingHeight', None)
        identifier['averagedHeight'] = stationData.get('averagedHeight', None)

    meanData = rawData.resample(samplingWindow).mean()
    if ~usePandas:
        meanData = meanData.compute()

    # if saveProperties is not None:
    #     query = dict(projectName=Identifier['projectName'],
    #                  start=Identifier['start'],
    #                  end=Identifier['end'],
    #                  samplingWindow=SamplingWindow,
    #                  station=Identifier['station'],
    #                  instrument=Identifier['instrument'],
    #                  height=Identifier['height']
    #                  )
    #
    #     doc = {}
    #     doc['projectName'] = query.pop('projectName')
    #     doc['dataFormat'] = saveProperties['dataFormat']
    #     doc['type'] = 'meteorological'
    #     doc['desc'] = query
    #     doc['desc']['start'] = identifier['start']
    #     doc['desc']['end'] = identifier['end']
    #     doc['desc']['samplingWindow'] = samplingWindow
    #     doc['desc']['params'] = list(meanData.columns)
    #     doc['resource'] = getSaveData(data=meanData, **saveProperties)
    #     datalayer.Cache.addDocument(**doc)


    return meanData


def getTurbulenceCalculatorFromDB(projectName, samplingWindow, start, end, usePandas=False, isMissingData=False, **kwargs):
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
    TurbulenceCalculator
        A turbulence calculator of the loaded raw data.
    """

    if type(start) is str:
        start = pandas.Timestamp(start)

    if type(end) is str:
        end = pandas.Timestamp(end)


    docList = datalayer.Measurements.getDocuments(projectName = projectName, **kwargs)
    dataList = [doc.getData(usePandas=usePandas) for doc in docList]

    rawData = pandas.concat(dataList) if usePandas else dask.dataframe.concat(dataList)
    rawData = rawData[start:end]

    identifier = {'projectName': projectName,
                  'samplingWindow': samplingWindow,
                  'station': None,
                  'instrument': None,
                  'height': None,
                  'start': start,
                  'end': end
                  }
    identifier.update(kwargs)

    projectData = datalayer.Project(projectName=projectName).getMetadata()[['height', 'instrument', 'station']].drop_duplicates()

    if identifier['station'] is not None:
        stationData = projectData.query("station=='%s'" % identifier['station']).iloc[0]
        identifier['buildingHeight'] = stationData.get('buildingHeight', None)
        identifier['averagedHeight'] = stationData.get('averagedHeight', None)

    return TurbulenceCalculator(rawData = rawData, metadata=projectData, identifier=identifier, isMissingData=isMissingData)


def getTurbulenceCalculatorFromData(data, samplingWindow, isMissingData=False):
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
    TurbulenceCalculator
        A turbulence calculator of the given data.
    """
    identifier = {'samplingWindow': samplingWindow
                  }

    return TurbulenceCalculator(rawData=data, metadata={}, identifier=identifier, isMissingData=isMissingData)


def getTurbulenceCalculator(data=None, projectName=None, **kwargs):
    if data is not None:
        return getTurbulenceCalculatorFromData(data=data, **kwargs)
    elif projectName is not None:
        return getTurbulenceCalculatorFromDB(projectName=projectName, **kwargs)
    else:
        raise ValueError("'data' argument or 'projectName' argument must be delivered")
