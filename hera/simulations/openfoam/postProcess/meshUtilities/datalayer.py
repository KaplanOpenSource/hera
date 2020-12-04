from .....datalayer import project
from .analysis import analysis
from .presentation import presentation
import os

class datalayer(project.ProjectMultiDBPublic):

    _projectName = None
    _analysis = None
    _presentation = None

    @property
    def presentation(self):
        return self._presentation

    @property
    def analysis(self):
        return self._analysis

    def __init__(self, projectName, databaseNameList=None, useAll=False,publicProjectName="postProcess"):

        self._projectName = projectName
        super().__init__(projectName=projectName, publicProjectName=publicProjectName,databaseNameList=databaseNameList,useAll=useAll)
        self._analysis = analysis(dataLayer=self)
        self._presentation = presentation(dataLayer=self,Analysis=self._analysis)

    def getDocuments(self,**kwargs):
        return self.getSimulationsDocuments(**kwargs)

    def getDocumentsAsDict(self,**kwargs):
        return self.getSimulationsDocumentsAsDict(**kwargs)

    def saveAndAddtoDB(self, save, addToDB, data, path, key, filter, user=None, **kwargs):

        path = path if path is not None else os.path.join("%s.hdf" % filter)
        if save:
            data.to_hdf(path, key=key, format="table")
            if addToDB:
                if user is None:
                    if self._databaseNameList[0] == "public" or self._databaseNameList[0] == "Public" and len(self._databaseNameList) > 1:
                        userName = self._databaseNameList[1]
                    else:
                        userName = self._databaseNameList[0]
                else:
                    userName = user
                self.addSimulationsDocument(desc=(dict(filter=filter, **kwargs)),resource=dict(path=path, key=key), type="OFsimulation", dataFormat="HDF",users=[userName])