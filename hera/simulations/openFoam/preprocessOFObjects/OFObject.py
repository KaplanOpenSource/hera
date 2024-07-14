import pandas
import os
import glob
from ....utils.logging import get_classMethod_logger
from .. import FIELDTYPE_VECTOR, FIELDTYPE_TENSOR, FIELDTYPE_SCALAR
# from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile,WriteParameterFile
# from PyFoam.Basics.DataStructures import Field,Vector,Tensor,DictProxy,Dimension
# from .utils import extractFieldFile,ParsedParameterFileToDataFrame

class OFObject:
    """
        Represents an openfoam field object for the preprocessing phase.
        Hence, this phase is used to set the boundary conditions and the initial conditions of the field

        Reading for posprocess is perfmed with the readFieldAsDataFrame of the toolkit, or with the VTK pipeline.
    """
    data = None
    name = None  # The name of the field
    fileName = None  # The name of the file ont the disk
    dimensions = None

    REGION_INTERNSALFIELD = 'internalField'
    REGION_BOUNDARYFIELD = 'boundaryField'

    @staticmethod
    def getDimensions(kg=0, m=0, s=0, K=0, mol=0, A=0, cd=0):
        """
            Returns the openfaom dimensions vector.
        Parameters
        ----------
        kg : int

        m  : int
        s  : int
        K  : int
        mol: int
        A  : int
        cd : int

        Returns
        -------
            The openfoam unit vector  [kg m s K mol A cd]
        """
        return f"[{kg} {m} {s} {K} {mol} {A} {cd}]"

    @property
    def componentNames(self):
        if self.fieldType == FIELDTYPE_SCALAR:
            ret = [self.fileName]
        elif self.fieldType == FIELDTYPE_VECTOR:
            ret = [f"{self.fileName}{cn}" for cn in ['x', 'y', 'z']]
        else:
            ret = [f"{self.fileName}{cn}" for cn in ['xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz']]

        return ret


    def internalField(self,processorName='singleProcessor'):
        """
        Return the interinal field data
        Returns
        -------

        """
        return self.data[processorName]['internalField']

    @property
    def processors(self):
        return self.data.keys()

    @property
    def processorItems(self):
        return self.data.items()


    @property
    def dimensionsStr(self):
        return self.getDimensions(**self.dimensions)

    @property
    def dimensionsList(self):
        return [ self.dimensions.get('kg',0),
                 self.dimensions.get('m',0),
                 self.dimensions.get('s',0),
                 self.dimensions.get('K',0),
                 self.dimensions.get('mol',0),
                 self.dimensions.get('A',0),
                 self.dimensions.get('cd',0)]

    def boundaryField(self,processorName='singleProcessor'):
        return self.data[processorName]['boundaryField']

    def __init__(self, name, fileName, fieldType, dimensions=None):
        """
            Initializes the OF field.

            Use the internalField as data is not None.

        Parameters
        ----------
        name: str
            The name of the object
        fileName: str
            The name of the file to write.
        fieldType: str
            The type of the field. scalar, vector or tensor for scalar, vector or tensor object.
        dimensions: dict
            The dict of the dimenstins.
        """
        logger = get_classMethod_logger(self, "init")
        logger.info(f"---------Start : {logger.name}")
        self.name = name
        self.fileName = fileName
        self.data = dict()

        if dimensions is not None:
            self.dimensions = dimensions

        if fieldType not in [FIELDTYPE_TENSOR, FIELDTYPE_SCALAR, FIELDTYPE_VECTOR]:
            err = f"Field type must be {','.join([FIELDTYPE_TENSOR, FIELDTYPE_SCALAR, FIELDTYPE_VECTOR])}. Got {fieldType}"
            logger.error(err)
            raise ValueError(err)

        self.fieldType = fieldType


    def writeToCase(self,caseDirectory,timeOrLocation):
        """
            Writes the current field to a case directory.

            The file will be written in parallel if it is parallel, and in
            single format if it is single.

        Parameters
        ----------
        caseDirectory
        timeOrLocation : float / str
            The name of the subdirectory to write in.

        Returns
        -------

        """
        if 'singleProcessor' in self.data:
            outputdir = os.path.join(caseDirectory,str(timeOrLocation),self.fileName)
            with open(outputdir,'w') as outputdir:
                outputdir.writelines(str(self.data['singleProcessor']))
        else:
            for procName,procData in self.data.items():
                outputdir = os.path.join(caseDirectory,procName,str(timeOrLocation),self.fileName)
                with open(outputdir,'w') as outputdir:
                    outputdir.writelines(str(procData).replace("proc.*",'"proc.*"'))

