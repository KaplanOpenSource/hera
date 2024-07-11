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
        "fieldType": str
            The type of the field. scalar, vector or tensor for scalar, vector or tensor object.

        boundaryPatchList : list
            The list of the patches

        data : pandas.DataFrame
            The internalField structure is :

            - Column 1 : the data of columns [1]
                ...
            - Column N : The data of column [N]
            - time     : The time step of the data
            - type      : internalField/boundaryField
            - boundaryName : The name of the boundary. None for the internalField
            - processor : Decomponsed case: processor[i] the name of the processor that holds that data.
                          Reconstructed case: will be None
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

    def writeEmptyField(self, caseDirectory, timeOrLocation, parallel=False, parallelBoundary=False):
        """
            Overwrites an empty new field (without the data)

            Currently we have two version of write.
            (write new and write that writes the contents of the field).

            TODO In the future we should merge the two write procedures.

        Paramaters
        ----------
        caseDirectory: str
                The path of the case
        location: str
                The location in the caseDirectory to update.
                Can be system, constant or the time

        fileName: str
            Optional field name.
                If None, use the name of the field.

        parallel: bool
                If true, check if there is a decomposed case

        parallelBoundary: bool
            If true adds the
                "proc.*" {
                    type processor;
                }
            to the boundary.

            It is necessary as a template file for the 0.parallel. Hence we need to write as single processor
            but with the parallel boundary


        Returns
        -------
            String with the new data. if parallel, a dict with the processor number as key and the
            value as the data.

        """
        logger = get_classMethod_logger(self, "write")
        logger.execution("---- Start ----")
        logger.debug(
            f"caseDirectory {caseDirectory}, location {timeOrLocation}, parallel {parallel}, boundary parallel {parallelBoundary}")
        fileName = self.fileName
        data = '0' if self.fieldType == FIELDTYPE_SCALAR else f"({' '.join(['0' for x in self.componentNames])})"

        if parallel:
            logger.execution("Saving fields to parallel case")

            procPaths = [proc for proc in glob.glob(os.path.join(caseDirectory, "processor*"))]
            logger.debug(f"The processor found in the case : {','.join(procPaths)}")
            for processorName in procPaths:
                fullFileName = os.path.join(caseDirectory, processorName, str(timeOrLocation), fileName)
                self._writeNew(filename=fullFileName, data=data, parallel=parallel,
                               parallelBoundary=parallelBoundary)

        else:
            fullFileName = os.path.join(caseDirectory, str(timeOrLocation), fileName)
            logger.debug(f"Writinge {fullFileName} as composed file")

            self._writeNew(filename=fullFileName, data=data, parallel=parallel, parallelBoundary=parallelBoundary)

        logger.execution("---- End ----")

    def _getHeader(self):

        fieldType = self.fieldType

        return """
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  10
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       vol{fieldType}Field;
    location    "0";
    object      {name};
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    """.replace("{fieldType}", fieldType.title()).replace("{name}", self.name)

    def _writeNew(self, filename, data, parallel, parallelBoundary=False):
        raise NotImplementedError("Implemented specifically for field or list")

