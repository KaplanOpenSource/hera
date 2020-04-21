import pandas
from ...datalayer import Measurements


class DataLayer(object):
    _columnsDescDict = None

    @property
    def columnsDescDict(self):
        return self._columnsDescDict

    def __init__(self):
        self._columnsDescDict = dict(RH='Relative humidity',
                                     Latitude='Latitude',
                                     Longitude='Longitude'
                                     )

    def loadData(self, projectName, locationName, date, filePath, **kwargs):
        """
        Loads a radiosonde data to the database as JSON_pandas from a csv file.

        :param projectName: The project name
        :param locationName: The measurement location
        :param date: The measurement date
        :param filePath: csv data file path
        :param kwargs: Other description arguments
        :return:
        """
        data = pandas.read_csv(filePath)

        desc = dict(locationName=locationName,
                    date=date,
                    columns=list(data.columns),
                    columnsDesc=self.columnsDescDict,
                    DataSource='radiosonde'
                    )
        desc.update(kwargs)

        doc = dict(projectName=projectName,
                   resource=data.to_json(),
                   dataFormat='JSON_pandas',
                   type='meteorological',
                   desc=desc
                   )
        Measurements.addDocument(**doc)

    def getData(self, projectName, **kwargs):
        """
        Returns metadata documents list according to the input requirements.

        :param projectName:
        :param kwargs:
        :return:
        """
        docList = Measurements.getDocuments(projectName=projectName,
                                            dataFormat='JSON_pandas',
                                            type='meteorological',
                                            DataSource='radiosonde',
                                            **kwargs
                                            )
        return docList


class AnalysisLayer(object):
    def __init__(self):
        pass


class PresentationLayer(object):
    def __init__(self):
        pass