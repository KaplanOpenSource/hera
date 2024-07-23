import pandas
import numpy
import os
import glob
from itertools import product
from ....utils import loadJSON
from ....utils.logging import get_classMethod_logger
from .. import FIELDTYPE_VECTOR, FIELDTYPE_TENSOR, FIELDTYPE_SCALAR, FIELDCOMPUTATION_EULERIAN, \
    FIELDCOMPUTATION_LAGRANGIAN,FLOWTYPE_INCOMPRESSIBLE,FLOWTYPE_COMPRESSIBLE
from .OFField import OFField

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

    FLOWTYPE_INCOMPRESSIBLE = FLOWTYPE_INCOMPRESSIBLE
    FLOWTYPE_COMPRESSIBLE   = FLOWTYPE_COMPRESSIBLE


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
    def pandasToFoamFormat(self,data):
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


    def getEmptyField(self, fieldName, flowType,noOfProc = None, addParallelProc = False):
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
                      fieldComputation=fieldData['fieldComputation'],noOfProc = noOfProc, addParallelProc = addParallelProc)
        return ret


    def getEmptyFieldFromCase(self,fieldName, flowType,caseDirectory,internalValue=0, readParallel=True ):
        """
            Reads the field structure (processors, boundary fields) from case, but not the data.


        Parameters
        ----------

        fieldName : string
            The name of the field

        flowType: str
            Compressible/incompressible.

        caseDirectory : string
            The directory  of the case

        internalValue : float
            Initialize the field. should be a list if the field is a vector

        readParallel : bool
            If false, read as single processor even if parallel case exists.

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
        dimensions = fieldData['dimensions'].get(flowType, fieldData['dimensions'].get('default', None))
        fileName = self.fieldDefinitions[fieldName].get("fileName", fieldName)

        ret = OFField(name=fieldName, fileName=fileName, dimensions=dimensions, fieldType=fieldData['fieldType'],
                      fieldComputation=fieldData['fieldComputation'],initialize=False)

        ret.readBoundariesFromCase(caseDirectory, readParallel=readParallel,internalValue=internalValue)
        return ret

    def readFieldFromCase(self, fieldName, flowType,caseDirectory,timeStep=0, readParallel=True ):
        """
            Returns a field object, load the data from the case.

        Parameters
        ----------

        fieldName : string
            The name of the field

        flowType: str
            Compressible/incompressible.

        caseDirectory : string
            The directory  of the case

        timeStep : float
            The time step to read

        readParallel : bool
            If false, read as single processor even if parallel case exists.
        Returns
        -------

        """
        logger = get_classMethod_logger(self, "readFieldFromCase")
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
