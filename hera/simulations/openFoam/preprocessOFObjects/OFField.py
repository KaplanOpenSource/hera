import pandas
import os
import glob
from ....utils.logging import get_classMethod_logger
from .. import FIELDTYPE_VECTOR,  FIELDTYPE_SCALAR
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile,WriteParameterFile
from PyFoam.RunDictionary.BoundaryDict import BoundaryDict
from PyFoam.Basics.DataStructures import Field,Vector,Tensor,DictProxy,Dimension
from .OFObject import OFObject
from .utils import ParsedParameterFileToDataFrame


class OFField(OFObject):

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
        return self.data[item]

    def __iter__(self):
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
            return Vector(0,0,0)
        else: # tensor
            return Tensor(0,0,0,0,0,0,0,0,0)

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

            retDict = WriteParameterFile(name=self.fileName, className=f"vol{self.fieldType.title()}Field")
            retDict['dimensions'] = Dimension(*self.dimensionsList)
            newField = Field(self.getZeroValue())
            newField.uniform = True
            retDict['internalField'] = newField

            boundaryValue = DictProxy()
            boundaryValue['type'] = 'zeroGradient'

            boundaryDict = DictProxy()
            boundaryDict['".*"'] =boundaryValue
            retDict['boundaryField'] = boundaryDict

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

                retDict = WriteParameterFile(name=self.fileName, className=f"vol{self.fieldType.title()}Field")
                retDict['dimensions'] = Dimension(*self.dimensionsList)
                newField = Field(self.getZeroValue())
                newField.uniform = True
                retDict['internalField'] = newField

                boundaryValue = DictProxy()
                boundaryValue['type'] = 'zeroGradient'

                boundaryDict = DictProxy()
                boundaryDict['".*"'] = boundaryValue
                retDict['boundaryField'] = boundaryDict

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
        boundaryField = DictProxy()
        for key, value in kwargs.items():
            boundaryField[key] = value

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
        boundaryField = DictProxy()
        boundaryField['type'] = "processor"

        if 'singleProcessor' in self.data:
            self.data['singleProcessor']['boundaryField']["proc.*"] = boundaryField
        else:
            for procName,dictData in self.data.items():
                self.data[procName]['boundaryField']["proc.*"] = boundaryField




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
                bndry = BoundaryDict(caseDirectory)
                retDict = WriteParameterFile(name=self.fileName, className=f"vol{self.fieldType.title()}Field")
                retDict['dimensions'] = Dimension(*self.dimensionsList)
                newField = Field(internalValue)
                newField.uniform = True
                retDict['internalField'] = newField

                logger.debug("Setting boundary field")
                boundaryDict = DictProxy()

                if 'singleProcessor' in self.data:
                    logger.debug("The field was already initialized. Checking if the boundaries already exist. ")
                    for boundaryName, boundaryData in self.data['singleProcessor']['boundaryField'].items():
                        logger.debug(f"Setting the boundary {boundaryName} to existing value")
                        boundaryDict[boundaryName] = boundaryData

                for boundaryName,boundaryData in bndry.getValueDict().items():
                    boundaryType = boundaryData['type']
                    interiorBoundary = True if boundaryType =='processor' else False
                    logger.debug(f"... boundary {boundaryName}. An interior boundary ? {interiorBoundary}")
                    boundaryValue = DictProxy()
                    boundaryValue['type'] = boundaryTypeDict.get(boundaryType,'zeroGradient')
                    boundaryDict[boundaryName] = boundaryValue

                retDict['boundaryField'] = boundaryDict
                self.data = dict(singleProcessor=retDict)
                self.parallel = False
            else:
                logger.debug("Setting as a parallel case")
                self.parallel = True
                for proc in procPaths:
                    bndry = BoundaryDict(proc)
                    retDict = WriteParameterFile(name=self.fileName, className=f"vol{self.fieldType.title()}Field")
                    retDict['dimensions'] = Dimension(*self.dimensionsList)
                    newField =  Field(internalValue)
                    newField.uniform = True
                    retDict['internalField'] = newField

                    boundaryDict = DictProxy()
                    for boundaryName,boundaryData in bndry.getValueDict().items():
                        logger.debug(f"... Setting {boundaryName}")
                        boundaryType = boundaryData['type']
                        interiorBoundary = True if boundaryType == 'processor' else False
                        logger.debug(f"... boundary {boundaryName}. An interior boundary ? {interiorBoundary}")
                        boundaryValue = DictProxy()
                        boundaryValue['type'] = boundaryTypeDict.get(boundaryType,'zeroGradient')
                        boundaryDict[boundaryName] = boundaryValue

                    retDict['boundaryField'] = boundaryDict

                    self.data[os.path.basename(proc)] = retDict
        except FileNotFoundError as e:
            err = f"File not found, initializing as an empty field: {e} "
            logger.warning(err)
            raise ValueError("Field not found")

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
                    singleProcessor=ParsedParameterFile(os.path.join(caseDirectory, str(timeStep), self.fileName)))
                self.parallel = False
            else:
                self.parallel = True
                for proc in procPaths:
                    self.data.update({os.path.basename(proc): ParsedParameterFile(os.path.join( proc, str(timeStep), self.fileName))})
        except FileNotFoundError as e:
            err = f"File not found, initializing as an empty field: {e} "
            logger.warning(err)
            raise ValueError("Field not found")

    def setFieldFromDataFrame(self,fieldDataFrame):
        """
            Gets a pandas/dask dataframe and sets the OFField to hold this data.

            The of field will hold it as a ParsedParameterFile of PyFoam

        Returns
        -------

        """
        logger = get_classMethod_logger(self,"setFieldFromDataFrame")
        if 'processor' in fieldDataFrame.columns:
            logger.info("updating field as parallel")
            for procName,procData in fieldDataFrame.groupby("processor"):
                internalData = fieldDataFrame.query(f"processor=={procName}")
                self.data[f"processor{procName}"] = self._processorToPyFoam(internalData,processorName=f"processor{procName}")
        else:
            logger.info("updating field as singleProcessor")
            self.data = dict(singleProcessor = self._processorToPyFoam(fieldDataFrame,processorName='singleProcessor'))

    def _processorToPyFoam(self,processordataFrame,processorName):
        """
            Converts the data of a single processor to a ParsedParameterFile of openFOAM.

        Parameters
        ----------
        processordataFrame

        Returns
        -------

        """
        def dataframeToField(regionData):
            name = ""
            if self.fieldType == FIELDTYPE_SCALAR:
                fieldList = [l[0] for l in regionData.sort_values("processorIndex")[self.componentNames].values]
                name = 'List<scalar>'
            elif self.fieldType == FIELDTYPE_VECTOR:
                fieldList = [Vector(*l) for l in regionData.sort_values("processorIndex")[self.componentNames].values]
                name = 'List<vector>'
            else:
                fieldList = [Tensor(*l) for l in regionData.sort_values("processorIndex")[self.componentNames].values]
                name = 'List<tensor>'

            newField =  Field(fieldList)
            newField.name = name
            return newField

        logger = get_classMethod_logger(self,"_processorToPyFoam")
        retDict = WriteParameterFile(name=self.fileName,className=f"vol{self.fieldType.title()}Field")
        retDict['dimensions'] = Dimension(*self.dimensionsList)

        logger.info("Setting internal field")
        internalField = processordataFrame.query("region=='internalField'")
        if len(internalField) > 0 :
            logger.debug("Setting the internal field to a new value")
            retDict['internalField'] = dataframeToField(internalField)
        else:
            logger.debug("Getting the existing internal field")
            retDict['internalField'] = self.data[processorName]['internalField']

        logger.info("Setting boundary field")
        boundaryDict = DictProxy()
        for boundaryName, boundaryData in self.data[processorName]['boundaryField'].items():
            logger.debug(f"Setting the boundary {boundaryName} to existing value")
            boundaryDict[boundaryName] = boundaryData

        for boundaryName, boundaryData in processordataFrame.query("region=='boundaryField'").groupby("boundary"):
            logger.debug(f"Setting the boundary {boundaryName} to new  value. data of length {len(boundaryData)}")
            boundaryValue = DictProxy()
            boundaryValue['type'] = 'fixedValue'
            boundaryValue['value'] = dataframeToField(boundaryData)
            boundaryDict[boundaryName] = boundaryValue

        if self.parallel:
            processorValue = DictProxy()
            processorValue['type'] = 'processor'
            boundaryDict['proc.*'] = processorValue

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
