import glob
import os

import numpy
import pandas
import tqdm
from foamlib import DimensionSet, FoamFieldFile, FoamFile

from ....utils.logging import get_classMethod_logger
from .. import FIELDTYPE_SCALAR, FIELDTYPE_VECTOR
from .OFObject import OFObject
from .utils import ParsedParameterFileToDataFrame, readFieldFile


def toFoamValue(value):
    """
        Converts a user value to the form foamlib writes.

        Accepts a number, a sequence ([5, 0, 0]) or OpenFOAM text
        ("uniform (5 0 0)", "fixedValue"), which is parsed the way OpenFOAM reads it.
    """
    if isinstance(value, str):
        return FoamFieldFile.loads(f"value {value};".encode())["value"]
    if isinstance(value, (list, tuple)):
        return numpy.asarray(value, dtype=float)
    return value


def readBoundaryFile(caseDirectory):
    """Returns the patches of constant/polyMesh/boundary as a dict name -> patch dict."""
    return dict(FoamFile(os.path.join(caseDirectory, "constant", "polyMesh", "boundary"))[None])


class OFField(OFObject):
    """OpenFOAM field with initialization, I/O, and boundary manipulation."""

    fieldComputation = None  # eulerian/ lagrangian.

    def __init__(self, name,fileName, fieldType, fieldComputation, dimensions=None,noOfProc=None,addParallelProc=False,initialize=True):
        """
            Initializing the OpenFOAM Field.

        Parameters
        ----------
        name : : str
            The name of the field
        fileName : str
            The name of the file to use.
        fieldType : str
            The type : scalar, vector or tensor.
            Use the constants FIELDTYPE_
        fieldComputation : str
            Eulerian/Lagrangian. Important for the loading of the field.
        dimensions : dicts
            The dict of the dimensions.
        noOfProc : None, number
            The number of processors.
            If the number is None, use single processor,
            otherwise initialize with the number of processors.
        addParallelProc : bool
            If true, add the parallel fields.
        initialize : bool
            If true, initializes the field, else leaves it empty.
            (the user can initialize it in the future).
        """
        super().__init__(name=name,fileName=fileName, fieldType=fieldType, dimensions=dimensions)
        self.fieldComputation = fieldComputation
        if initialize:
            self.initialize(noOfProc=noOfProc)
            if addParallelProc:
                self.addProcBoundary()


    def __getitem__(self, item):
        """Return the field data for the given processor key."""
        return self.data[item]

    def __iter__(self):
        """Iterate over processor keys in the field data."""
        return self.data.__iter__()

    def getZeroValue(self):
        """
            Returns a zero value of a scalar, vector or tensor
        Returns
        -------

        """
        if self.fieldType == FIELDTYPE_SCALAR:
            return 0
        elif self.fieldType == FIELDTYPE_VECTOR:
            return numpy.zeros(3)
        else: # tensor
            return numpy.zeros(9)

    def _newProcessorData(self, internalValue):
        """The dict of one processor: dimensions, a uniform internalField and no boundaries."""
        return {"dimensions": DimensionSet(*self.dimensionsList),
                "internalField": toFoamValue(internalValue),
                "boundaryField": {}}

    def initialize(self,noOfProc=None):
        """
            Initializes a field that is either single or multiple processor.
            If None, initialize as a single processor.

            Initializes the object as a general boundary:
            that is :
            the boundary is ".*" with type zeroGradient for all the
            processors.

        Parameters
        ----------
        noOfProc : int or None
            If None, initialize as a single processor

        Returns
        -------

        """
        logger = get_classMethod_logger(self,"initialize")
        if noOfProc is None:
            logger.debug(f"Initializing field {self.name} as single processor")
            self.parallel = False

            retDict = self._newProcessorData(self.getZeroValue())
            retDict['boundaryField']['".*"'] = {'type': 'zeroGradient'}

            self.data = dict(singleProcessor=retDict)

        else:
            try:
                noOfProc = int(noOfProc)
                if noOfProc < 1:
                    raise ValueError("must be positive")
            except Exception:
                err = f"noOfProc must be a positive integer (greater than 0). Got {noOfProc}"
                logger.error(err)
                raise ValueError(err)

            logger.debug(f"Initializing field {self.name} as a multi-processor")
            self.parallel = True
            procDict = dict()
            for procID in range(noOfProc):
                procName = f"processor{procID}"

                retDict = self._newProcessorData(self.getZeroValue())
                retDict['boundaryField']['".*"'] = {'type': 'zeroGradient'}

                procDict[procName] = retDict

            self.data = procDict

    def setInternalUniformFieldValue(self,value):
        """
            Sets the internal field values.
        Returns
        -------
            value : float, list
                Sets the value of the internal field.
                must be a flaot for scalar, 3-list and 9-list for vector and tensor.
        """
        value = toFoamValue(value)
        if 'singleProcessor' in self.data:
            self.data['singleProcessor']['internalField'] = value
        else:
            for procName,dictData in self.data.items():
                self.data[procName]['internalField'] = value

    def addBoundaryField(self,boundaryName,**kwargs):
        """
            Sets the boundary name with the key/value
        Parameters
        ----------
        boundaryName : The name of the boundary to set
        kwargs :
            key/value of the names.

        Returns
        -------

        """
        logger = get_classMethod_logger(self,"addBoundaryField")
        boundaryField = {key: toFoamValue(value) for key, value in kwargs.items()}

        if 'singleProcessor' in self.data:
            self.data['singleProcessor']['boundaryField'][boundaryName] = boundaryField
        else:
            for procName,dictData in self.data.items():
                self.data[procName]['boundaryField'][boundaryName] = boundaryField

    def addProcBoundary(self):
        """
            Adds the 'proc.*' boundary of typ processor to all the boundary fields.
        Returns
        -------

        """
        if 'singleProcessor' in self.data:
            self.data['singleProcessor']['boundaryField']['"proc.*"'] = {'type': 'processor'}
        else:
            for procName,dictData in self.data.items():
                self.data[procName]['boundaryField']['"proc.*"'] = {'type': 'processor'}


    def readBoundariesFromCase(self, caseDirectory,internalValue=None, readParallel=True):
        """
            Reads the boundaries from the case.
            If the mesh is parallel, the code sets zeroGradient for all the external boundaries,
            and a processor type for all the internal boundaries.

        Parameters
        ----------

        caseDirectory : str
            The directory
        internalValue : float or list
            The value of the internalField.

        timeStep : float
            The timestep to use for getting the directory.

        readParallel : bool
            Force read parallel if exists.

        Returns
        -------

        """
        logger = get_classMethod_logger(self, "readBoundariesFromCase")
        logger.debug(f"---- Start {logger.name}")

        readSingle = False
        if readParallel:
            logger.debug("Trying to read as a parallel case")

            procPaths = [proc for proc in glob.glob(os.path.join(caseDirectory, "processor*"))]
            if len(procPaths) == 0:
                logger.debug("Case is not parallel. Try to read it as single processor")
                readSingle = True
        else:
            readSingle = True

        if internalValue is None:
            internalValue = self.getZeroValue()
            logger.debug(f"Internal value is None, getting the zero of the current field: {internalValue}")

        boundaryTypeDict = dict(processor='processor',
                                symmetry='symmetry')

        try:
            if readSingle:
                logger.debug("Setting as a single processor")
                bndry = readBoundaryFile(caseDirectory)
                retDict = self._newProcessorData(internalValue)

                logger.debug("Setting boundary field")
                boundaryDict = retDict['boundaryField']

                if 'singleProcessor' in self.data:
                    logger.debug("The field was already initialized. Checking if the boundaries already exist. ")
                    for boundaryName, boundaryData in self.data['singleProcessor']['boundaryField'].items():
                        logger.debug(f"Setting the boundary {boundaryName} to existing value")
                        boundaryDict[boundaryName] = boundaryData

                for boundaryName,boundaryData in bndry.items():
                    boundaryType = boundaryData['type']
                    interiorBoundary = True if boundaryType =='processor' else False
                    logger.debug(f"... boundary {boundaryName}. An interior boundary ? {interiorBoundary}")
                    boundaryDict[boundaryName] = {'type': boundaryTypeDict.get(boundaryType,'zeroGradient')}

                self.data = dict(singleProcessor=retDict)
                self.parallel = False
            else:
                logger.debug("Setting as a parallel case")
                self.parallel = True
                for proc in procPaths:
                    bndry = readBoundaryFile(proc)
                    retDict = self._newProcessorData(internalValue)

                    for boundaryName,boundaryData in bndry.items():
                        logger.debug(f"... Setting {boundaryName}")
                        boundaryType = boundaryData['type']
                        interiorBoundary = True if boundaryType == 'processor' else False
                        logger.debug(f"... boundary {boundaryName}. An interior boundary ? {interiorBoundary}")
                        retDict['boundaryField'][boundaryName] = {'type': boundaryTypeDict.get(boundaryType,'zeroGradient')}

                    self.data[os.path.basename(proc)] = retDict
        except FileNotFoundError as e:
            err = f"File not found: {e}"
            logger.warning(err)
            # If reading in parallel mode and the case directory is empty, fall back to a single processor empty field.
            if readParallel:
                self.initialize(noOfProc=None)
                self.parallel = False
                return
            # Otherwise, propagate as a missing field error.
            raise ValueError("Field not found") from e

    def readFromCase(self, caseDirectory, timeStep=0, readParallel=True):
        """
            Reads one time step of the field from the case directory.

            If the file does not exist,

            if readParallel is true, tries to load as parallel if directories exist.
            else, reads only as single.

        Parameters
        ----------
        caseDirectory : string
            The case directory

        time : list, float
            Reads the field from the time

        readParallel : bool
            If True, tries to read as parallel if case exists.

        Returns
        -------

        """
        logger = get_classMethod_logger(self, "readFromCase")
        logger.debug(f"---- Start {logger.name}")

        readSingle = False
        if readParallel:
            logger.debug("Trying to read as a parallel case")

            procPaths = [proc for proc in glob.glob(os.path.join(caseDirectory, "processor*"))]
            if len(procPaths) == 0:
                logger.debug("Case is not parallel. Try to read it as single processor")
                readSingle = True
        else:
            readSingle = True

        try:
            if readSingle:
                self.data = dict(
                    singleProcessor=readFieldFile(os.path.join(caseDirectory, str(timeStep), self.fileName)))
                self.parallel = False
            else:
                self.parallel = True

                effective_level_int = logger.getEffectiveLevel()
                if effective_level_int == logger.DEBUG:
                    itrObj = tqdm.tqdm(procPaths)
                else:
                    itrObj = procPaths
                for proc in procPaths:
                    self.data.update({os.path.basename(proc): readFieldFile(
                        os.path.join(proc, str(timeStep), self.fileName))})

        except FileNotFoundError as e:
            err = f"File not found: {e}"
            logger.warning(err)
            # If reading in parallel mode and the case directory is empty, fall back to a single processor empty field.
            if readParallel:
                self.initialize(noOfProc=None)
                self.parallel = False
                return
            # Otherwise, propagate as a missing field error.
            raise ValueError("Field not found") from e

    def setFieldFromDataFrame(self,fieldDataFrame):
        """
            Gets a pandas/dask dataframe and sets the OFField to hold this data.

            The field holds it as a dict per processor, the way foamlib represents an OpenFOAM file.

        Returns
        -------

        """
        logger = get_classMethod_logger(self,"setFieldFromDataFrame")
        if 'processor' in fieldDataFrame.columns:
            logger.info("updating field as parallel")
            # Reset data to avoid retaining singleProcessor entry from previous state
            self.data = {}
            self.parallel = True
            for procName,procData in fieldDataFrame.groupby("processor"):
                internalData = fieldDataFrame.query(f"processor=={procName}")
                self.data[f"processor{procName}"] = self._processorToFoam(internalData,processorName=f"processor{procName}")
        else:
            logger.info("updating field as singleProcessor")
            self.data = dict(singleProcessor = self._processorToFoam(fieldDataFrame,processorName='singleProcessor'))

    def _processorToFoam(self,processordataFrame,processorName):
        """
            Converts the data of a single processor to the dict of an OpenFOAM field file.

        Parameters
        ----------
        processordataFrame

        Returns
        -------

        """
        def dataframeToField(regionData):
            """Convert a region DataFrame to a nonuniform field: an (N,) array for scalars, (N, 3) or (N, 9) otherwise."""
            values = regionData.sort_values("processorIndex")[self.componentNames].to_numpy(dtype=float)
            return values[:, 0] if self.fieldType == FIELDTYPE_SCALAR else values

        logger = get_classMethod_logger(self,"_processorToFoam")
        retDict = {'dimensions': DimensionSet(*self.dimensionsList)}

        logger.info("Setting internal field")
        internalField = processordataFrame.query("region=='internalField'")
        if len(internalField) > 0 :
            logger.debug("Setting the internal field to a new value")
            retDict['internalField'] = dataframeToField(internalField)
        else:
            logger.debug("Getting the existing internal field")
            retDict['internalField'] = self.data[processorName]['internalField']

        logger.info("Setting boundary field")
        boundaryDict = {}
        existing_boundary = self.data.get(processorName, {}).get('boundaryField', {})
        for boundaryName, boundaryData in existing_boundary.items():
            logger.debug(f"Setting the boundary {boundaryName} to existing value")
            boundaryDict[boundaryName] = boundaryData

        # Process boundary fields if present
        if "boundary" in processordataFrame.columns:
            for boundaryName, boundaryData in processordataFrame.query("region=='boundaryField'").groupby("boundary"):
                logger.debug(f"Setting the boundary {boundaryName} to new  value. data of length {len(boundaryData)}")
                boundaryDict[boundaryName] = {'type': 'fixedValue', 'value': dataframeToField(boundaryData)}
        else:
            logger.debug("No boundary column present; skipping boundary processing.")

        if self.parallel:
            boundaryDict['"proc.*"'] = {'type': 'processor'}

        retDict['boundaryField'] = boundaryDict

        return retDict

    def getDataFrame(self,filterInternalPatches=True):
        """
            Returns the data in the field as dataframe.
        Returns
        -------

        """
        ret = []
        if 'singleProcessor' in self.data:
            ret.append(ParsedParameterFileToDataFrame(self.data['singleProcessor'],columnNames=self.componentNames,filterInternalPatches=filterInternalPatches))
        else:
            for processorName,procParsedFile in self.data.items():
                ret.append(
                    ParsedParameterFileToDataFrame(self.data[processorName], columnNames=self.componentNames,filterInternalPatches=filterInternalPatches,processor=int(processorName[9:])))

        return pandas.concat(ret)
