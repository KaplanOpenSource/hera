import pandas
import dask.dataframe
from ..analytics.turbulencecalculator import TurbulenceCalculator
from .... import datalayer


def getTurbulenceCalculatorFromDB(projectName, samplingWindow, start, end, usePandas=False, isMissingData=False, **kwargs):
    """
    This method loads the raw data that corresponds to the requirements (projectName, station, instrument.. ) and
    creates a turbulence calculator with the desirable sampling window.

    :param projectName:     The name of the project.
    :param samplingWindow:  The desirable sampling window (as string).
    :param start:           Datetime of the begin.
    :param end:             Datetime of the end.
    :param usePandas:       A flag of whether or not to use pandas.
    :param isMissingData    A flag if there is a missing data to compute accordingly.

    :return: A turbulence calculator of the loaded raw data.
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
        stationData = projectData['stations'].query("station=='%s'" % identifier['station']).iloc[0]
        identifier['buildingHeight'] = stationData.get('buildingHeight', None)
        identifier['averagedHeight'] = stationData.get('averagedHeight', None)

    return TurbulenceCalculator(rawData = rawData, metadata=projectData, identifier=identifier, isMissingData=isMissingData)


def getTurbulenceCalculatorFromData(data, samplingWindow, isMissingData=False):
    """
    This method gets turbulence calculator

    :param data:           The raw data for the calculations.
    :param samplingWindow: The sampling window.
    :param isMissingData:  A flag if there is a missing data to compute accordingly.

    :return: A turbulence calculator of the loaded raw data.
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