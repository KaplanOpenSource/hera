import pandas
from ... import datalayer


class DataLayer(datalayer.Project):
    _columnsDescDict = None

    @property
    def columnsDescDict(self):
        return self._columnsDescDict


    def __init__(self, projectName):
        """
            Initializes a datalayer for the Radiosonde data.

        Parameters
        ----------

        projectName: str
                The project name
        """
        super().__init__(projectName=projectName)

        self._columnsDescDict = dict(RH='Relative humidity',
                                     Latitude='Latitude',
                                     Longitude='Longitude')

    def loadData(self, locationName, date, filePath, **kwargs):
        """
        Loads a radiosonde data to the database as JSON_pandas from a csv file.

        Parameters
        ----------

        locationName : str
            The measurement location

        date : pandas.Timestamp / str
            The measurement date

        filePath : str
            csv data file path

        kwargs :
            Other description arguments

        Returns
        -------

        """
        data = pandas.read_csv(filePath)

        if type(date) is str:
            date = pandas.Timestamp(date)

        desc = dict(locationName=locationName,
                    date=date,
                    columns=list(data.columns),
                    columnsDesc=self.columnsDescDict,
                    DataSource='Radiosonde'
                    )
        desc.update(kwargs)

        doc = dict(resource=data.to_json(),
                   dataFormat='JSON_pandas',
                   type='meteorology',
                   desc=desc
                   )
        self.addMeasurementsDocument(**doc)

    def getDocFromDB(self, projectName, **query):
        """
        Returns list of metadata documents according to the input requirements.

        Parameters
        ----------
        projectName : str
            The project name

        query :
            Other query arguments

        Returns
        -------
        list
            List of the documents that fulfill the query
        """
        docList = self.getMeasurementsDocuments(dataFormat='JSON_pandas',
                                                type='meteorology',
                                                DataSource='Radiosonde',
                                                **query)
        return docList


class AnalysisLayer(object):
    def __init__(self):
        pass


class PresentationLayer(object):
    def __init__(self):
        pass