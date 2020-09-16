from hera import GIS
import pathlib

class GISResources():

    _GISdatalayer = None
    _FilesDirectory = None
    _projectName = None

    @property
    def projectName(self):
        return self._projectName

    @projectName.setter
    def projectName(self, newName):
        self._projectName = newName
        self._GISdatalayer = GIS.GIS_datalayer(projectName=self._projectName, FilesDirectory=self._FilesDirectory)

    @property
    def FilesDirectory(self):
        return self._FilesDirectory

    @FilesDirectory.setter
    def FilesDirectory(self, newDirectory):
        self._FilesDirectory = newDirectory
        self._GISdatalayer = GIS.GIS_datalayer(projectName=self._projectName, FilesDirectory=self._FilesDirectory)

    def __init__(self, projectName="GISdata", FilesDirectory=None):

        self._FilesDirectory = pathlib.Path().absolute() if FilesDirectory is None else FilesDirectory
        self._projectName = projectName
        self._GISdatalayer = GIS.GIS_datalayer(projectName=self._projectName, FilesDirectory=self._FilesDirectory)