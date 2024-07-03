import pandas
import numpy
import os
import glob
import dask
from itertools import product
from ...utils import loadJSON
from ...utils.logging import get_classMethod_logger
from . import FIELDTYPE_VECTOR, FIELDTYPE_TENSOR, FIELDTYPE_SCALAR, FIELDCOMPUTATION_EULERIAN, \
    FIELDCOMPUTATION_LAGRANGIAN,FLOWTYPE_INCOMPRESSIBLE,FLOWTYPE_COMPRESSIBLE
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile,WriteParameterFile

from PyFoam.Basics.DataStructures import Field,Vector,Tensor,DictProxy,Dimension

#########################################################################
#               Fields
#########################################################################
class OFObjectHome:
    """
        A home for the openfoam preprocess field objects.
        Represents an openfoam field object for the preprocessing phase.
        Hence, this phase is used to set the boundary conditions and the initial conditions of the field

        Reading for posprocess is perfmed with the readFieldAsDataFrame of the toolkit, or with the VTK pipeline.
    """
    @property
    def fieldDefinitions(self):
        return self._fieldDefinitions

    def __init__(self):

        fieldJSON = """
        {
            "U" : { 
                "dimensions" : {
                    "default" : {
                        "m" :1, 
                        "s" :-1
                    }
                }, 
                "fieldType":"vector",
                "fieldComputation":"eulerian"
            },
            "p" : { 
                "dimensions" : {
                    "incompressible" : {
                        "m" :2, 
                        "s" :-2
                    },
                    "compressible" : {
                        "kg" : 1,
                        "m" : -1,
                        "s"  : -2
                    }
                }, 
                "fieldType":"scalar",
                "fieldComputation":"eulerian"
            },
            "p_rgh" : { 
                "dimensions" : {
                    "incompressible" : {
                        "m" :2, 
                        "s" :-2
                    },
                    "compressible": {
                        "kg" : 1,
                        "m" : -1,
                        "s"  : -2
                    }
                }, 
                "fieldType":"scalar",
                "fieldComputation":"eulerian"
            },
            "epsilon" : { 
                "dimensions" : {
                    "default" : {
                        "m" :2, 
                        "s" :-3
                    }
                }, 
                "fieldType":"scalar",
                "fieldComputation":"eulerian"
            },
            "nut" : { 
                "dimensions" : {
                    "default" : {
                        "m" :2, 
                        "s" :-1
                    }
                }, 
                "fieldType":"scalar",
                "fieldComputation":"eulerian"
            },
            "k" : { 
                "dimensions" : {
                    "default" : {
                        "m" :2, 
                        "s" :-2
                    }
                }, 
                "fieldType":"scalar",
                "fieldComputation":"eulerian"
            },
            "T" : { 
                "dimensions" : {
                    "default" : {
                        "K" :1
                    }
                },
                "fieldType" : "scalar",
                "fieldComputation":"eulerian"
            },
            "Tbackground" : { 
                "dimensions" : {
                    "default" : {
                        "K" :1
                    }
                }, 
                "fieldType":"scalar",
                "fieldComputation":"eulerian"
            },  
            "cellCenters" : {
                "fileName" : "C", 
                "dimensions" : {
                    "default" : {
                        "m" : 1
                    }
                }, 
                "fieldType":"vector",
                "fieldComputation":"eulerian"
            },
            "Hmix" : { 
                "dimensions" : {
                    "default" : {
                        "m" : 1
                    }
                }, 
                "fieldType":"scalar",
                "fieldComputation":"eulerian"
            },
            "ustar" : { 
                "dimensions" : {
                    "default" : {
                        "m" : 1,
                        "s"  : -1
                    }
                }, 
                "fieldType":"scalar",
                "fieldComputation":"eulerian"
            },            
            "CellHeights" : { 
                "dimensions" : {
                    "default" : {
                        "m" : 1
                    }
                }, 
                "fieldType":"scalar",
                "fieldComputation":"eulerian"
            }
        }"""
        self._fieldDefinitions = loadJSON(fieldJSON)

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

    @staticmethod
    def pandasToFoamFormat(data):
        """
            Converts pandas to a list of values in OF style.
            i.e

            [the number of records]
            (
                value1
                .
                .
            )

            where value can be a number (scalar) or a '(x0 x1 x2)' for a vector

        Parameters
        ----------
        data: pandas.DataFrame
            The data to convert

        Returns
        -------
            str.
        """
        D = data if self.componentNames is None else data[self.componentNames]

        newStr = f"{str(data.shape[0])}\n"
        newStr += "(\n"
        newStr += "\n".join([f"({x})" for x in D.to_csv(sep=' ', header=False, index=False).split("\n")[:-1]])
        newStr += "\n);\n"
        return newStr

    def addFieldDefinitions(self, fieldName, dimensions, fieldType, fieldComputation=FIELDCOMPUTATION_EULERIAN,
                            compressible_dimensions=None, overwrite=False):
        """
            Adds the field to the field home.
            TODO: Move the entire definition to be a config/datasource and so fields are appended for each project.

        Parameters
        ----------
        name : string
            The name of the field.
        dimensions : dict
            The dictionary of the dimensions.

        fieldType :  string
            scalar, vector or tensor. Use the FIELDTYPE constants.

        fieldComputation : string
            Can be (FIELDCOMPUTATION_EULERIAN or FIELDCOMPUTATION_LAGRANGIAN)
            Determine if the field is a solver field (eulerian) or a cloud field (lagrangian).


        compressible_dimensions : dict
            If the field has different units in the compressible case, pass them here.
            The dimensions will be the incompressible units.

        overwrite : bool
            If true, overwite the exisitng filed. Otherwise throw exception ValueError if the field exists.

        Returns
        -------
        OFField object.
        """
        if fieldName not in self.fieldDefinitions or overwrite:
            fielddimensions = dict()
            if compressible_dimensions is None:
                fielddimensions['default'] = dimensions
            else:
                fielddimensions['compressible'] = compressible_dimensions
                fielddimensions['incompressible'] = dimensions

            if fieldComputation not in [FIELDCOMPUTATION_LAGRANGIAN, FIELDCOMPUTATION_EULERIAN]:
                raise ValueError(
                    f"{fieldComputation} must be {FIELDCOMPUTATION_LAGRANGIAN} or {FIELDCOMPUTATION_EULERIAN}")

            if fieldType not in [FIELDTYPE_SCALAR, FIELDTYPE_VECTOR, FIELDTYPE_TENSOR]:
                raise ValueError(
                    f"The {fieldType} must be one of {','.join([FIELDTYPE_SCALAR, FIELDTYPE_VECTOR, FIELDTYPE_TENSOR])}")

            self.fieldDefinitions[fieldName] = dict(dimensions=fielddimensions, fieldType=fieldType,
                                                    fieldComputation=fieldComputation)
        else:
            raise ValueError(f"{fieldName} already exists!. Use overwrite=True to overwrite its definition")

    def getEmptyField(self, fieldName, flowType):
        """
            Return the field object with its dimensions.
            Since the dimensions of pressure change for compressible/incompressible
            solution, we haveto supply the group of parameters to which the group belongs.

            Currently holds a list of fields. In the future might read itfrom a list.

        Parameters
        ----------
        fieldName: str
            The field name

        flowType: str
            Compressible/incompressible.

        Returns
        -------
            OFField.

        """
        logger = get_classMethod_logger(self, "getEmptyField")
        logger.info(f"----- Start : {logger.name}")
        if fieldName not in self.fieldDefinitions.keys():
            err = f"Field {fieldName} not found. Must supply {','.join(self.fieldDefinitions.keys())}"
            logger.critical(err)
            raise ValueError(err)

        fieldData = self.fieldDefinitions[fieldName]
        dimensions = fieldData['dimensions'].get(flowType, fieldData['dimensions'].get('default',None))
        fileName = self.fieldDefinitions[fieldName].get("fileName", fieldName)

        ret = OFField(name=fieldName, fileName=fileName, dimensions=dimensions, fieldType=fieldData['fieldType'],
                      fieldComputation=fieldData['fieldComputation'])
        return ret

    def getFieldFromCase(self, fieldName, flowType,caseDirectory,timeStep=0, readParallel=True ):
        """
            Returns a field object, load the data from the case.

        Parameters
        ----------
        caseDirectory : string
            The directory  of the case

        flowType: str
            Compressible/incompressible.

        fieldName : string
            The name of the field

        Returns
        -------

        """
        logger = get_classMethod_logger(self, "getEmptyFieldFromCase")
        logger.info(f"----- Start : {logger.name}")
        if fieldName not in self.fieldDefinitions.keys():
            err = f"Field {fieldName} not found. Must supply {','.join(self.fieldDefinitions.keys())}"
            logger.critical(err)
            raise ValueError(err)

        fieldData = self.fieldDefinitions[fieldName]
        dimensions = fieldData['dimensions'].get(flowType, fieldData['dimensions']['default'])
        fileName = self.fieldDefinitions[fieldName].get("fileName", fieldName)

        ret = OFField(name=fieldName, fileName=fileName, dimensions=dimensions, fieldType=fieldData['fieldType'],
                      fieldComputation=fieldData['fieldComputation'])

        ret.readFromCase(caseDirectory,timeStep=timeStep, readParallel=readParallel)
        return ret



    def readFieldAsDataFrame(self, fieldName, caseDirectory, times=0, readParallel=True,filterInternalPatches=False):
        """
            Extracts a field to the disk from the requested times.
            If None, then reads from all the time steps.
            Reads only the internal field (and not the boundaries).

            If read in parallel, then add the processor to the output columns.

            Determine if the field is lagrnagian or eulerian and reads respectiveluy

        Parameters
        ----------
        caseDirectory
        times
        readParallel
        kwargs :
            Specialized fields

            Eulerian:
                - filterInternalPatches : bool [default True]
                    remove the proc* boundaries.
            Lagrangian:
                -  cloudName : string [default : kinematicCloud]
                    The default cloud name

        Returns
        -------

        """
        finalCasePath = os.path.abspath(caseDirectory)

        field = self.getEmptyField(fieldName=fieldName, flowType=FLOWTYPE_INCOMPRESSIBLE) # the type is important for the dimensions that are not considered here

        if readParallel:
            processorList = [os.path.basename(proc) for proc in glob.glob(os.path.join(finalCasePath, "processor*"))]
            if len(processorList) == 0:
                raise ValueError(f"There are no processor* directories in the case {finalCasePath}. Is it parallel?")

            if times is None:
                timeList = sorted([x for x in os.listdir(os.path.join(finalCasePath, processorList[0])) if (
                        os.path.isdir(os.path.join(finalCasePath, processorList[0], x)) and
                        x.isdigit() and
                        (not x.startswith("processor") and x not in ["constant", "system", "rootCase", 'VTK']))],
                                  key=lambda x: int(x))
            else:
                timeList = numpy.atleast_1d(times)

            data = pandas.concat([extractFieldFile(os.path.join(finalCasePath, processorName, str(timeName), field.fileName),
                                           columnNames=field.componentNames,
                                           time=timeName,filterInternalPatches=filterInternalPatches,
                                           processor=int(processorName[9:])) for processorName, timeName in
                 product(processorList, timeList)])
        else:

            if times is None:
                timeList = sorted([x for x in os.listdir(finalCasePath) if (
                        os.path.isdir(x) and
                        x.isdigit() and
                        (not x.startswith("processor") and x not in ["constant", "system", "rootCase", 'VTK']))],
                                  key=lambda x: int(x))
            else:
                timeList = numpy.atleast_1d(times)

            data = pandas.concat([extractFieldFile(os.path.join(finalCasePath, str(timeName), field.fileName),
                                           columnNames=field.componentNames,filterInternalPatches=filterInternalPatches,
                                           time=timeName) for timeName in timeList])

        return data


class OFObject:
    """
        Represents an openfoam field object for the preprocessing phase.
        Hence, this phase is used to set the boundary conditions and the initial conditions of the field

        Reading for posprocess is perfmed with the readFieldAsDataFrame of the toolkit, or with the VTK pipeline.
    """
    data = None
    name = None  # The name of the field
    fileName = None  # The name of the file ont the disk
    parallel = None

    dimensions = None

    REGION_INTERNSALFIELD = 'internalField'
    REGION_BOUNDARYFIELD = 'boundaryField'


    @property
    def componentNames(self):
        if self.fieldType == FIELDTYPE_SCALAR:
            ret = [self.fileName]
        elif self.fieldType == FIELDTYPE_VECTOR:
            ret = [f"{self.fileName}{cn}" for cn in ['x', 'y', 'z']]
        else:
            ret = [f"{self.fileName}{cn}" for cn in ['xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz']]

        return ret

    @property
    def internalField(self):
        """
        Return the interinal field data
        Returns
        -------

        """
        return self.data['internalField']

    @property
    def dimensionsStr(self):
        return OFObjectHome.getDimensions(**self.dimensions)

    @property
    def dimensionsList(self):
        return [ self.dimensions.get('kg',0),
                 self.dimensions.get('m',0),
                 self.dimensions.get('s',0),
                 self.dimensions.get('K',0),
                 self.dimensions.get('mol',0),
                 self.dimensions.get('A',0),
                 self.dimensions.get('cd',0)]


    @property
    def boundaryField(self):
        return self.data.get('boundaryField',None)

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
    """.replace("{fieldType}", fieldType).replace("{name}", self.name)

    def _writeNew(self, filename, data, parallel, parallelBoundary=False):
        raise NotImplementedError("Implemented specifically for field or list")


class OFField(OFObject):
    fieldComputation = None  # eulerian/ lagrangian.

    def __init__(self, name,fileName, fieldType, fieldComputation, dimensions=None):
        """
            Initializing the OpenFOAM Field.

            The internalField can be supplied either as a:
                - list (with the values)
                - a pandas .

            The boundary field as a:
                - dict of patch->[list/pandas],

        Parameters
        ----------
        name : str
            The name of the field
        "fieldType" :
            The type : scalar, vector or tensor.
        dimensions : dict
            The dict of parameters.
        boundaryPatch : list
            The name of the patches

        data : pandas.DataFrame
            The data of the internal field.

        boundaryField : pandas.DataFrame
            The
        """
        super().__init__(name=name,fileName=fileName, fieldType=fieldType, dimensions=dimensions)
        self.fieldComputation = fieldComputation

    def _writeNew(self, filename, data, parallel=False, parallelBoundary=False):
        """
            Inject the data to a new file.

        Parameters
        ----------
        filename : str
                The full file name
        data: pandas.Dataframe or pandas.Series or number.
                The data to inject
        dimensions: str
                The dimensions of the variable.
                Here just accept OF full format. i.e. [0 1 0 0 0 0 1]  ..

        Returns
        -------

        """
        logger = get_classMethod_logger(self, "_writeNew")
        logger.execution("----- Start -----")

        fileStrContent = self._getHeader()
        fileStrContent += "\n\n" + f"dimensions {self.dimensionsStr};\n\n"

        if type(data) in [float, int, list, tuple, str]:  # isinstance(data,float) or isinstance(data,int):
            if type(data) in [float, int, str]:
                fileStrContent += f"internalField  uniform {data};\n"
            else:
                fileStrContent += f"internalField  uniform ({' '.join([str(x) for x in data])});\n"
        else:
            fileStrContent += "internalField   nonuniform List<vector>\n"

            if isinstance(data, pandas.Series):
                componentNames = ['demo']
            else:
                componentNames = [x for x in data.columns if
                                  (
                                              x != 'processor' and x != 'time')] if self.componentNames is None else self.componentNames

            if len(componentNames) > 1:
                # vector
                fileStrContent += self.pandasToFoamFormat(data)

            else:
                # scalar
                if isinstance(data, pandas.Series):
                    data = data.values
                fileStrContent += ""
                fileStrContent = f"{str(data.shape[0])}\n"
                fileStrContent += "(\n"
                fileStrContent += "\n".join(data)
                fileStrContent += ");\n"

        fileStrContent += self._getBoundaryConditionsStr(parallel=parallelBoundary)
        with open(filename, 'w') as outfile:
            outfile.write(fileStrContent)

        return fileStrContent

    def _getBoundaryConditionsStr(self, parallel=False):
        boundaryConditions = """
boundaryField
{
"""
        if self.boundaryField is not None:
            for boundaryPatchName, boundaryData in self.boundaryField.items():
                bstr = f"\n{boundaryPatchName}\n"
                bstr += "{\n"
                for bcondProp, bcondData in boundaryData.items():
                    bstr += f"\t{bcondProp} {bcondData};\n"
                bstr += "}\n"
                boundaryConditions += bstr

        if parallel:
            boundaryConditions += """
"proc.*"
{
    type            processor;
}
                    """
        else:
            boundaryConditions += """
    ".*"
    {
        type            zeroGradient;
    }
"""
        boundaryConditions += """
}        
 """
        return boundaryConditions

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
        logger = get_classMethod_logger(self, "write")
        logger.execution(f"---- Start {logger.name}")

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

            Note that it overwrite the object, so in parallel case, the internal boundariesmust not be changed.
            In the future, we can take the old field values and just overwrite the existing bounadaries.
            rather than creating a whole new field.

        Returns
        -------

        """
        logger = get_classMethod_logger(self,"setFieldFromDataFrame")
        if 'processor' in fieldDataFrame.columns:
            logger.info("updating field as parallel")
            for procName,procData in fieldDataFrame.groupby("processor"):
                self.data[f"processor{procName}"] = self._processorToPyFoam(fieldDataFrame.query(f"processor=={procName}"))

        else:
            logger.info("updating field as singleProcessor")
            self.data = dict(singleProcessor = self._processorToPyFoam(fieldDataFrame))

    def _processorToPyFoam(self,processordataFrame):
        """
            Converts the data of a single processor to a ParsedParameterFile of openFOAM.

        Parameters
        ----------
        processordataFrame

        Returns
        -------

        """
        def dataframeToField(regionData):
            if self.fieldType == FIELDTYPE_SCALAR:
                fieldList = [l for l in regionData.sort_values("processorIndex")[self.componentNames].values]
            elif self.fieldType == FIELDTYPE_VECTOR:
                fieldList = [Vector(*l) for l in regionData.sort_values("processorIndex")[self.componentNames].values]
            else:
                fieldList = [Tensor(*l) for l in regionData.sort_values("processorIndex")[self.componentNames].values]
            newField =  Field(fieldList)
            return newField


        retDict = WriteParameterFile(name=self.fileName)
        retDict['dimensions'] = Dimension(*self.dimensionsList)
        retDict['internalField'] = dataframeToField(processordataFrame.query("region=='internalField'"))
        boundaryDict = DictProxy()
        for boundaryName, boundaryData in processordataFrame.query("region=='boundaryField'").groupby("boundary"):
             boundaryDict[boundaryName] = dataframeToField(boundaryData)

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

class OFList(OFObject):
    """
        Just data.
    """

    def _updateExisting(self, filename, data, parallel=False):
        """
            Just rewrite the field.

            This function exists to complete the interface similarly to field.


        Parameters
        ----------
        filename: str
            The file name
        data : str

        Returns
        -------

        """
        return self._writeNew(filename, data, parallel=parallel)

    def _writeNew(self, filename, data, parallel=False):
        """
            Writes an OF list file.

        Parameters
        ----------
        filename : str
                The name of the file

        data: pandas.DataFrame or pandas.Series
                Holds the data

        columnNames: list [optional]
                The list of names to use. If None, use all.

        Returns
        -------
            str,
        """
        if isinstance(data, pandas.Series):
            columnNames = ['demo']
        else:
            columnNames = [x for x in data.columns if
                           (x != 'processor' and x != 'time')] if self.columnNames is None else self.columnNames

        fileStrContent = self.getHeader()
        if len(columnNames) > 1:
            # vector
            fileStrContent += self.pandasToFoamFormat(data, columnNames)

        else:
            # scalar
            if isinstance(data, pandas.Series):
                fileStrContent += "\n".join(data)
            else:
                fileStrContent += "\n".join(data[columnNames])

        with open(filename, 'w') as outfile:
            outfile.write(fileStrContent)

        return fileStrContent


#########################################################################
#
#               Mesh handling
#########################################################################
class OFMeshBoundary:
    _case = None
    _checkIfParallel = None

    _boundaryNames = None

    @property
    def case(self):
        return self._case

    def __init__(self, directory: str, checkParallel: bool = True):
        """
            Reads the boundary from the boundary file of the mesh.

            If checkParallel is true, then checks if the parallel case exists. If it does,
            read the boundary from the parallel case. Else, read it from the constant.



        Parameters
        ----------
        baseFile
        """
        self._case = directory
        self._checkIfParallel = checkParallel

        self._boundaryNames = []
        if self._checkIfParallel and os.path.exists(os.path.join(self.case, "processor0")):
            for proc in glob.glob(os.path.join(self.case, "processor*")):
                self._boundaryNames += self._readBoundary(os.path.join(proc, "constant", "polyMesh", "boundary"))
        else:
            self._boundaryNames = self._readBoundary(os.path.join(directory, "constant", "poly", "boundary"))

    def getBoundary(self, filterProcessor: bool = True):
        """
            Return a list of all the boundaries. If fileterprocessor is true,
            remove all the processor*.
        Parameters
        ----------
        filterProcessor : bool
            If true remove all the processor* faces from the list.

        Returns
        -------

        """
        return list(set([x for x in self._boundaryNames if
                         'procBoundary' not in x] if filterProcessor else self._boundaryNames))

    def _readBoundary(self, boundaryFile):
        """
                Reads the boundary file and extracts the boundary names.
        Parameters
        ----------
        boundaryFile

        Returns
        -------

        """

        def isInt(line):
            try:
                return int(line)
            except:
                return None

        with open(boundaryFile, "r") as inFile:
            data = inFile.readlines()

        firstLine = [i for i, x in enumerate(data) if isInt(x) is not None][0]
        braceList = [i for i, x in enumerate(data[firstLine:]) if x.strip() == '{']
        return [data[firstLine:][x - 1].strip() for x in braceList]


#########################################################################
#
#               Utils
#########################################################################
def extractFieldFile(casePath, columnNames, patchNameList=None,filterInternalPatches=True, **kwargs):
    try:
        dataParsedFile = ParsedParameterFile(casePath)
    except Exception as e:
        print(casePath)
        raise ValueError(e)
    return ParsedParameterFileToDataFrame(columnNames=columnNames, patchNameList=patchNameList,filterInternalPatches=filterInternalPatches, **kwargs)

def ParsedParameterFileToDataFrame(dataParsedFile,columnNames, patchNameList=None,filterInternalPatches=True, **kwargs):
    ret = []
    pndsData = pandas.DataFrame([[x for x in item] for item in dataParsedFile['internalField'].val],
                                columns=columnNames).assign(**kwargs, region='internalField')

    ret.append(pndsData)
    for patchName in dataParsedFile['boundaryField']:

        if patchNameList is not None:
            addPatch = True if patchName in patchNameList else False
        else:
            addPatch = True

        if filterInternalPatches and 'proc' in patchName:
            addPatch = False

        if addPatch:
            pndsData = pandas.DataFrame(
                [[x for x in item] for item in dataParsedFile['boundaryField'][patchName]['value'].val],
                columns=columnNames).assign(**kwargs, region='boundaryField', boundary=patchName)
        ret.append(pndsData)
    return pandas.concat(ret).reset_index().rename(columns=dict(index="processorIndex"))

#     def emptyParallelField(self, caseDirectory,timeName:"0.parallel",processor:"", data=None,boundaryField=None):
#         """
#             Writes a null file with definitions of the processor boundries
#             because sometimes, the snappyHex mesh of decomposed pars breaks the parallel
#             fields. This file is used to correct this problem.
#
#         Parameters
#         ----------
#         filename: str
#             The name of the file to write to.
#
#         data: float, str, pandas.DataFrame, dask.DataFrame
#             The data to write to the field.
#
#         Returns
#         -------
#
#         """
#         if data is None:
#             data = '0' if self.componentNames is None else f"({' '.join(['0' for x in self.componentNames])})"
#
#         fileStrContent = self._getHeader()
#         fileStrContent += "\n\n" + f"dimensions {self.dimensions};\n\n"
#
#         if type(data) in [float,int,list,tuple,str]: #isinstance(data,float) or isinstance(data,int):
#             if type(data) in [float,int,str]:
#                 fileStrContent += f"internalField  uniform {data};\n"
#             else:
#                 fileStrContent += f"internalField  uniform ({' '.join([str(x) for x in data])});\n"
#         else:
#             fileStrContent += "internalField   nonuniform List<vector>\n"
#
#             if isinstance(data,pandas.Series):
#                 componentNames = ['demo']
#             else:
#                 componentNames = [x for x in data.columns if (x != 'processor' and x != 'time')] if self.componentNames is None else self.componentNames
#
#             if len(componentNames) > 1:
#                 # vector/tensor
#                 fileStrContent += self.pandasToFoamFormat(data,componentNames)
#             else:
#                 # scalar
#                 if isinstance(data, pandas.Series):
#                     data = data.values
#                 fileStrContent += ""
#                 fileStrContent = f"{str(data.shape[0])}\n"
#                 fileStrContent += "(\n"
#                 fileStrContent += "\n".join(data)
#                 fileStrContent += ");\n"
#
#         boundaryConditions = ""
#         if boundaryField is not None:
#             for boundaryPatchName,boundaryData in boundaryField.items():
#                 bstr = f"\n{boundaryPatchName}\n"
#                 bstr += "{\n"
#                 for bcondProp,bcondData in boundaryData.items():
#                     bstr += f"\t{bcondProp} {bcondData};\n"
#                 bstr += "}\n"
#                 boundaryConditions += bstr
#
#         fileStrContent += """
#    boundaryField
# {
#     "proc.*"
#     {
#         type            processor;
#     }
#     """  + boundaryConditions + """
# }
# """
#         filename = os.path.join(caseDirectory,processor, timeName, self.name)
#         logger.debug(f"Saving the file {self.name} to {filename} ")
#         with open(filename,'w') as outfile:
#             outfile.write(fileStrContent)
#
# def extractFile(path, columnNames, vector=True, skiphead = 20,skipend = 4):
#     """
#         Extracts data from a openFOAM list file.
#
#         list files has no boundary and so, we can just skip head and end.
#
#     Parameters
#     ----------
#     path: str
#         The path of the file
#     time: str
#         The files' time step
#     columnNames: list of str
#         The names of the columns
#     skiphead: int
#         Number of lines to skip from the beginning of the file
#     skipend: int
#         Number of lines to skip from the ending of the file
#
#     Returns
#     -------
#         Pandas with the data.
#     """
#
#     cnvrt = lambda x: float(x.replace("(", "").replace(")", ""))
#     cnvrtDict = dict([(x, cnvrt) for x in columnNames])
#
#     try:
#         newData = pandas.read_csv(path,
#                                   skiprows=skiphead,
#                                   skipfooter=skipend,
#                                   engine='python',
#                                   header=None,
#                                   delim_whitespace=True,
#                                   converters=cnvrtDict,
#                                   names=columnNames)
#     except ValueError:
#         newData = []
#
#     if len(newData) == 0:
#         with open(path, "r") as thefile:
#             lines = thefile.readlines()
#
#         vals = lines[17]
#         data = []
#
#         if vector:
#             if "{" in vals:
#                 inputs = vals.split("{")
#                 repeat = int(inputs[0])
#                 valuesList = inputs[1][inputs[1].find("(") + 1:inputs[1].find(")")]
#                 data = dict(
#                     [(colname, [float(x)] * repeat) for colname, x in zip(columnNames, valuesList.split(" "))])
#             else:
#                 for rcrdListTuple in vals.split("(")[2:]:
#                     record = dict(
#                         [(name, float(y)) for name, y in zip(columnNames, rcrdListTuple.split(")")[0].split(" "))])
#                     data.append(record)
#         else:
#             if "{" in vals:
#                 inputs = vals.split("{")
#                 repeat = int(inputs[0])
#                 value = float(inputs[1].split("}")[0])
#                 data = [{columnNames[0]: value} for x in range(repeat)]
#
#             else:
#                 valuesList = vals.split("(")[1]
#                 for rcrdListItem in valuesList.split(" "):
#                     record = {columnNames[0]: float(rcrdListItem.split(")")[0])}
#                     data.append(record)
#
#         newData = pandas.DataFrame(data)
#
#     return newData.astype(float)
#
# def extractBoundaryField(path,boundaryName, columnNames, **kwargs):
#     """
#         Return the boundry condition as a dict:
#         - typenames
#         - value : pandas.
#             pandas with the column names.
#         - [other values]
#
#     Parameters
#     ----------
#     path : string
#         The dir
#     boundaryName : string
#         the name of the boundary
#
#     columnNames :
#         The type of the data (None for scalar).
#     kwargs :
#         Additional attributes that will be appended to the values pandas.
#
#     Returns
#     -------
#         dict
#     """
#     # 1. find the number of points in the file.
#     logger = logging.getLogger("extractBoundaryField")
#     L = []
#     resDict = dict()
#     with open(path, "r") as thefile:
#         lines = thefile.readline()
#         try:
#             if lines.strip() == boundaryName:
#                 lineCount = True
#             else:
#                 lineCount = False
#         except ValueError:
#             lineCount = None
#
#         while lineCount:
#             if lines.strip() == boundaryName:
#                 lineCount = True
#             else:
#                 lineCount = False
#
#         logger.debug(f"Found the boundary {boundaryName}")
#         alldata = []
#         while '}' not in lines:
#             prsed =[x.replace("\n"," ").strip() for x in lines.expandtabs().split(" ") if len(x.strip())>0]
#             logger.debug(f"Parsing the next line: {prsed}")
#
#             if 'value' in prsed[0]:
#
#
#                 if 'nonuniform' in prsed[1]:
#
#                 elif 'uniform' in prsed[1]:
#
#                 else:
#                     resDict[prsed[0]] = " ".join(prsed)
#
#             else:
#                 resDict[prsed[0]] =prsed[-1]
#
#             lines = thefile.readline()
#
#         for line in islice(thefile, 1, lineCount+1):
#           L.append(line.replace("(","").replace(")",""))
#
#     return pandas.read_csv(StringIO("".join(L)),names=columnNames,sep:" ").assign(**kwargs)
