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
from .utils import extractFieldFile

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
        """Return the dictionary of field definitions."""
        return self._fieldDefinitions

    def __init__(self):
        """Initialize the field definitions from the built-in JSON catalog."""
        # The field catalog is an inline JSON string that maps each OpenFOAM field name
        # (e.g. "U", "p", "k") to its metadata: dimensions, fieldType, and fieldComputation.
        #
        # Each field entry has the following structure:
        #   - "dimensions": a dictionary keyed by flow regime. Fields whose units do not
        #     change between compressible and incompressible solvers use a single "default"
        #     key. Fields whose physical dimensions differ (notably pressure "p" and "p_rgh")
        #     carry separate "incompressible" and "compressible" keys so the correct SI
        #     exponents are selected at runtime via the flowType argument.
        #     The dimension dictionaries themselves map SI unit symbols ("kg", "m", "s", "K",
        #     "mol", "A", "cd") to their integer exponents; only non-zero exponents are listed.
        #   - "fieldType": one of "scalar", "vector", or "tensor" — controls how values are
        #     written and read (single number vs. 3-component tuple vs. 9-component tuple).
        #   - "fieldComputation": "eulerian" (mesh-based solver field) or "lagrangian"
        #     (particle cloud field). All fields in this built-in catalog are eulerian.
        #   - "fileName" (optional): the actual file name on disk when it differs from the
        #     logical field name (e.g. "cellCenters" is stored as "C").
        #
        # Solver field definitions below cover the standard CFD quantities:
        #   U         — velocity vector [m/s]
        #   p, p_rgh  — pressure fields with flow-type-dependent dimensions
        #   epsilon   — turbulent dissipation rate [m^2/s^3]
        #   omega     — specific dissipation rate [1/s]
        #   alphat    — turbulent thermal diffusivity [kg/(m*s)]
        #   nut       — turbulent kinematic viscosity [m^2/s]
        #   k         — turbulent kinetic energy [m^2/s^2]
        #   T, Tbackground — temperature fields [K]
        #   cellCenters    — mesh cell centre coordinates (file "C") [m]
        #   Hmix, ustar, distanceFromWalls — auxiliary fields used by specific solvers
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
            "omega" : { 
                "dimensions" : {
                    "default" : {
                        "s" :-1
                    }
                }, 
                "fieldType":"scalar",
                "fieldComputation":"eulerian"
            },
            "alphat" : { 
                "dimensions" : {
                    "default" : {
                        "kg" : 1,
                        "m" :-1, 
                        "s" :-1
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
            "distanceFromWalls" : { 
                "dimensions" : {
                    "default" : {
                        "m" : 1
                    }
                }, 
                "fieldType":"vector",
                "fieldComputation":"eulerian"
            }
        }"""
        # Parse the JSON catalog string into a Python dict and store it as the
        # authoritative registry of field definitions for this instance.
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
        # If componentNames is set (i.e. vector/tensor field), select only those
        # columns from the DataFrame; for scalar fields (componentNames is None)
        # use the full DataFrame as-is.
        D = data if self.componentNames is None else data[self.componentNames]

        # Build the OpenFOAM list format:
        #   <numberOfEntries>
        #   (
        #       (x0 x1 x2)   <-- vector values wrapped in parentheses
        #       ...
        #   );
        # First line: the record count (number of rows).
        newStr = f"{str(data.shape[0])}\n"
        # Opening parenthesis of the OpenFOAM list block.
        newStr += "(\n"
        # Convert the DataFrame to space-separated CSV (no header, no index),
        # split into individual lines (dropping the trailing empty line produced
        # by the CSV serialiser), and wrap each line in parentheses to form
        # OpenFOAM vector notation "(v0 v1 v2)".
        newStr += "\n".join([f"({x})" for x in D.to_csv(sep=' ', header=False, index=False).split("\n")[:-1]])
        # Closing parenthesis with a semicolon terminator as required by OF syntax.
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
        # Validate that the requested field exists in the catalog.
        if fieldName not in self.fieldDefinitions.keys():
            err = f"Field {fieldName} not found. Must supply {','.join(self.fieldDefinitions.keys())}"
            logger.critical(err)
            raise ValueError(err)

        fieldData = self.fieldDefinitions[fieldName]
        # Resolve the correct dimension dictionary based on flowType.
        # The catalog stores dimensions under flow-type-specific keys
        # ("incompressible" / "compressible") for fields like pressure whose
        # units depend on the solver formulation, or under a single "default"
        # key for fields with flow-type-independent units (e.g. velocity "U").
        # We first try an exact match on flowType; if that key is absent we
        # fall back to "default". If neither exists, dimensions will be None.
        dimensions = fieldData['dimensions'].get(flowType, fieldData['dimensions'].get('default',None))
        # The on-disk file name may differ from the logical field name (e.g.
        # "cellCenters" is stored as file "C"). Fall back to fieldName itself.
        fileName = self.fieldDefinitions[fieldName].get("fileName", fieldName)

        ret = OFField(name=fieldName, fileName=fileName, dimensions=dimensions, fieldType=fieldData['fieldType'],
                      fieldComputation=fieldData['fieldComputation'],noOfProc = noOfProc, addParallelProc = addParallelProc)
        return ret


    def getEmptyFieldFromCase(self,fieldName, flowType,caseDirectory,internalValue=None, readParallel=True ):
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
                      fieldComputation=fieldData['fieldComputation'],initialize=False)

        ret.readFromCase(caseDirectory,timeStep=timeStep, readParallel=readParallel)
        return ret

    def readFieldAsDataFrame(self, fieldName, caseDirectory, times=0, readParallel=True,filterInternalPatches=False):
        """
            DEPRACATED. read the field and use getDataFrame.

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

        # Construct a lightweight field object to obtain the file name and component
        # names. The flowType is set to incompressible because only the metadata
        # (not the dimensions) matters for reading raw data.
        field = self.getEmptyField(fieldName=fieldName, flowType=FLOWTYPE_INCOMPRESSIBLE) # the type is important for the dimensions that are not considered here

        # --- Parallel case detection ---
        # A parallel OpenFOAM case stores results in processor0/, processor1/, ...
        # directories. If readParallel is True we glob for these directories to
        # build the processor list and verify that at least one exists.
        if readParallel:
            processorList = [os.path.basename(proc) for proc in glob.glob(os.path.join(finalCasePath, "processor*"))]
            if len(processorList) == 0:
                raise ValueError(f"There are no processor* directories in the case {finalCasePath}. Is it parallel?")

            # --- Time list discovery (parallel) ---
            # When times=None, discover available time directories inside the first
            # processor folder. Only numeric directory names are kept (filtering out
            # "constant", "system", etc.), then sorted numerically.
            # When times is provided, wrap it with atleast_1d so a single scalar
            # value becomes iterable.
            if times is None:
                timeList = sorted([x for x in os.listdir(os.path.join(finalCasePath, processorList[0])) if (
                        os.path.isdir(os.path.join(finalCasePath, processorList[0], x)) and
                        x.isdigit() and
                        (not x.startswith("processor") and x not in ["constant", "system", "rootCase", 'VTK']))],
                                  key=lambda x: int(x))
            else:
                timeList = numpy.atleast_1d(times)

            # itertools.product generates the Cartesian product of processorList and
            # timeList, yielding every (processorName, timeName) combination. Each
            # pair maps to a unique field file path:
            #   <case>/processorN/<time>/<fieldFile>
            # All extracted DataFrames are concatenated into a single frame; the
            # processor index is recorded as an integer column stripped from the
            # "processorN" directory name (processorName[9:] drops the "processor" prefix).
            data = pandas.concat([extractFieldFile(os.path.join(finalCasePath, processorName, str(timeName), field.fileName),
                                           columnNames=field.componentNames,
                                           time=timeName,filterInternalPatches=filterInternalPatches,
                                           processor=int(processorName[9:])) for processorName, timeName in
                 product(processorList, timeList)])
        else:
            # --- Serial (single-processor) case ---
            # Time directories live directly under the case root instead of under
            # processor sub-directories.

            # Discover time directories the same way as in the parallel branch,
            # but scanning the case root directly.
            if times is None:
                timeList = sorted([x for x in os.listdir(finalCasePath) if (
                        os.path.isdir(x) and
                        x.isdigit() and
                        (not x.startswith("processor") and x not in ["constant", "system", "rootCase", 'VTK']))],
                                  key=lambda x: int(x))
            else:
                timeList = numpy.atleast_1d(times)

            # In serial mode the file path is simply <case>/<time>/<fieldFile>.
            # No processor column is added to the resulting DataFrame.
            data = pandas.concat([extractFieldFile(os.path.join(finalCasePath, str(timeName), field.fileName),
                                           columnNames=field.componentNames,filterInternalPatches=filterInternalPatches,
                                           time=timeName) for timeName in timeList])

        return data
