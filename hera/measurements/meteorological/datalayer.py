import hera.datalayer.project

import pydoc

class DataLayer(hera.datalayer.project.ProjectMultiDBPublic):
    _DataSource = None
    _parser = None
    _docType = 'meteorological'
    _np_size = "100Mb"


    def __init__(self, DataSource, projectName, publicProjectName, databaseNameList=None, useAll=False):
        """
            Initialization of the datalayer.

        Parameters
        ----------

        DataSource: str
                The name of the parser to load the data.

        projectName: str
                The project name at the local database(s).
        publicProjectName: str
                The project name at the public DB.
        databaseNameList: list
                A list of databases to search in.
        useAll: bool
                If true, return the first results.
        """
        super().__init__(projectName=projectName,
                         publicProjectName=publicProjectName,
                         databaseNameList=databaseNameList,
                         useAll=useAll)
        self.logger.info(f"Initializing a multiDBPublic project with name {projectName} and public name: {publicProjectName}")
        self._DataSource = DataSource
        classpath = f"hera.measurements.meteorological.parserClasses.Parser_{DataSource}"
        parserCls = pydoc.locate(classpath)
        self.logger.debug(f"Parser class is {classpath}")
        self._parser = parserCls()

    def getDocFromDB(self, resource=None, dataFormat=None, **desc):
        """
            Loads the data from the DB with the data source

        :param resource: str
                The resource to query on.
        :param dataFormat: str
                The dataFormat to query on.
        :param desc: dict
                MongoDB querying
        :return: list
                The doc list.
        """
        desc['DataSource'] = self._DataSource
        docList = self.getMeasurementsDocuments(resource=resource,
                                                dataFormat=dataFormat,
                                                type=self._docType,
                                                **desc)
        return docList

    def getDocFromFile(self, **kwargs):
        pass

    def parse(self, **kwargs):
        return self._parser().parse(**kwargs)

    def loadData(self, **kwargs):
        pass

