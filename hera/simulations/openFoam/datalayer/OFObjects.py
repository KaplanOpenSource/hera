import pandas
import numpy
import os
import glob
from dask.delayed import delayed
import dask
from itertools import product,islice
from io import StringIO
from ....utils import loggedObject

class ofObjectHome(loggedObject):
    """
        A factory for field classes.

        This class creates the field and populates the data.
    """

    GROUP_COMPRESSIBLE = "compressible"
    GROUP_INCOMPRESSIBLE = "incompressible"
    GROUP_DISPERSION = "dispersion"

    @property
    def predifinedFields(self):
        return self._predefinedfields

    def __init__(self):
        super().__init__(loggerName=None)

        incompressibleDict = dict(U=dict(dimensions=self.getDimensions(m=1, s=-1), componentNames=['Ux', 'Uy', 'Uz']),
                                  p=dict(dimensions=self.getDimensions(m=2, s=-2), componentNames=None),
                                  epsilon=dict(dimensions="[ 0 2 -3 0 0 0 0 ]", componentNames=None),
                                  nut=dict(dimensions="[ 0 2 -1 0 0 0 0 ]", componentNames=None),
                                  k=dict(dimensions="[ 0 2 -2 0 0 0 0 ]", componentNames=None),
                                  )

        compressibleDict = dict()

        dispersionDict = dict(Hmix=dict(dimensions=self.getDimensions(m=1), componentNames=None),
                              ustar=dict(dimensions=self.getDimensions(m=1,s=-1), componentNames=None),
                              CellHeights=dict(dimensions=self.getDimensions(m=1), componentNames=None))


        self._predefinedfields=dict(compressible=compressibleDict,
                      incompressible=incompressibleDict,
                      dispersion = dispersionDict)


    def getDimensions(self,kg=0,m=0,s=0,K=0,mol=0,A=0,cd=0):
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

    def getPredefinedFields(self,fieldGroup):
        """
            Return a list of all the fields in a gropu.
        Parameters
        ----------
        fieldGroup: str
            The group to list.

        Returns
        -------
            A list of names
        """
        return [x for x in self.predifinedFields[fieldGroup].keys()]

    def getField(self,fieldName,fieldGroup="incompressible",componentNames=None,dimensions=None):
        """
            Return the field with its dimensions.
            Since the dimensions of pressure change for compressible/incompressible
            solution, we haveto supply the group of parameters to which the group belongs.

            Currently holds a list of fields. In the future might read itfrom a list.

        Parameters
        ----------
        fieldName: str
            The field name

        fieldGroup: str
            The group to which the parameter belongs.
            Currently:
                - compressible
                - incompressible
                - dispersion

        componentNames: list of str
            The names of the components.
            The length of the list determines whether the field s scalar, vector or tensor.
            (for scalars, the componenetName is None).

        dimensions : dict
            The dimensions of the object.
            the key name is the unit name and the value is its exponent:
            dict(m=1,s=-1) is the units of velocity (m/s).

            The units are:
                kg,m,s,K,mol,A,cd

        Returns
        -------
            OFField.

        """
        self.logger.info("----- Start ----")

        if fieldName not in self.predifinedFields[fieldGroup]:
            self.logger.warning(f"{fieldName} not found in pre-existing fields. Using inputs")
            if dimensions is  None:
                raise ValueError(f"Must supply dimensions to the new field {fieldName}")
            objData = dict(dimensions= self.getDimensions(**dimensions),componentNames = componentNames)
        else:
            objData = dict(self.predifinedFields[fieldGroup][fieldName])

        objData['name'] = fieldName

        self.logger.info(f"Getting the object {objData}")
        self.logger.info("----- End -----")
        return OFField(**objData)


    def loadLagrangianDataParallel(self,
                                   casePath,
                                   times=None,
                                   parallelCase=True,
                                   withVelocity=False,
                                   withReleaseTimes=False,
                                   withMass=False,
                                   cloudName="kinematicCloud"):
        """
            Extracts results of an LSM run.

        Parameters
        -----------

            times: list
                 a list of time steps to extract.
                 If None, it extracts all time steps in the casePath.
            file: str
                The name for a file, in which the data is saved. Default is "<cloud name>_data.parquet", in the current working directory.

            withVelocity: bool
                    True: extract the particles' velocities in addition to their positions.
                    Default is False.

            **kwargs:
                Any additional parameters to add to the description in the DB.
        Returns
        --------
            document of the data.
        """
        finalCasePath = os.path.abspath(casePath)
        loader = lambda timeName: readLagrangianRecord(timeName,
                                                       casePath=finalCasePath,
                                                       withVelocity=withVelocity,
                                                       withReleaseTimes=withReleaseTimes,
                                                       withMass=withMass,
                                                       cloudName=cloudName)

        if parallelCase:

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

            data = dask.dataframe.from_delayed(
                [delayed(loader)(os.path.join(processorName, timeName)) for processorName, timeName in
                 product(processorList, timeList[1:])])
        else:

            if times is None:
                timeList = sorted([x for x in os.listdir(finalCasePath) if (
                        os.path.isdir(x) and
                        x.isdigit() and
                        (not x.startswith("processor") and x not in ["constant", "system", "rootCase", 'VTK']))],
                                  key=lambda x: int(x))
            else:
                timeList = numpy.atleast_1d(times)

            data = dask.dataframe.from_delayed([delayed(loader)(timeName) for timeName in timeList])

        return data





class OFObject(loggedObject):
    """
        Represents an OF object. i.e a field or a list.

        This will be the interface to handle fields in OF.

        It will provide the interface to load, store and create empty fields.

        Currently, the support for the boundary conditions is limited.

    """

    _name = None        # The name of the field
    _data = None        # The data (after loading).
    _componentNames = None # If none, the field is scalar, else, this is the list of fields.

    @property
    def name(self):
        return self._name

    @property
    def componentNames(self):
        return self._componentNames

    @property
    def data(self):
        return self.data

    @property
    def fieldType(self):
        fieldType = "Scalar"
        if self.componentNames is not None:
            if len(self.componentNames) == 3:
                fieldType = "Vector"
            elif len(self.componentNames) == 9:
                fieldType = "Tensor"
            else:
                raise ValueError(f"Component names for field {self.name} is invalid. ")
        return fieldType

    def __init__(self,name,componentNames=None):
        """
            Initializes the OF field.
        Parameters
        ----------
        name: str
            The name of the object

        componentNames: list, None
            If None, then it is a scalar.
            If list, then it is a list.
        """
        super().__init__(loggerName=None)
        self._name = name
        self._componentNames =componentNames

    def pandasToFoamFormat(self, data):
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

    def _getHeader(self):

        fieldType = self.fieldType

        return """
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  7
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       vol{fieldType}Field;
    location    "0";
    object      C;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    """.replace("{fieldType}",fieldType)

    def write(self, caseDirectory, location, data=None, fileName=None, parallel=False):
        """
            Updates the field internal data to the input.
            Should **not** update the boundaries data.


            If parallel is true, and there is a decomposed case, write the data to the different processors.
            In that case, the series must have processorNumber and index columns. Otherwise, only index columns is required.

            The new data is written to the disk

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

        data: pandas
                The data to update. Must include the column index (for composed cases) and processorName/index for decomposed cases.

        parallel: bool
                If true, check if there is a decomposed case

        Returns
        -------
            String with the new data. if parallel, a dict with the processor number as key and the
            value as the data.

        """
        self.logger.info("---- Start----")
        fileName = self.name if fileName is None else fileName

        if data is None:
            self.logger.debug("Data was not supplied, using 0. ")
            data = '0' if self.componentNames is None else f"({' '.join(['0' for x in self.componentNames])})"

        if parallel:
            self.logger.info("Saving fields to parallel case")

            procPaths = [proc for proc in glob.glob(os.path.join(caseDirectory, "processor*"))]

            if isinstance(data,pandas.DataFrame) or isinstance(data,dask.dataframe.DataFrame):

                if 'processor' not in data:
                    raise ValueError("data was read from composed case. Does not have multiprocessor structure.")

                dataProcessorList = data.processor.unique()
                if isinstance(data,dask.dataframe):
                    dataProcessorList = dataProcessorList.compute()

                if len(dataProcessorList) != len(procPaths):
                    errStr = f"Processor on disk and the data provided have mismatch in processor count. Disk: {len(procPaths)} ; Data {len(dataProcessorList)}"
                    raise ValueError(errStr)

            for processorName in procPaths:

                procID = int(processorName[9:])

                if isinstance(data, pandas.DataFrame) or isinstance(data, dask.dataframe):
                    procData = data.query("processor == @procID")
                    if isinstance(data, dask.dataframe):
                        procData = procData.compute()
                    else:
                        procData = data

                fullFileName = os.path.join(caseDirectory, processorName, str(location), fileName)

                if os.path.exists(fullFileName):
                    self._updateExisting(filename=fullFileName,data=procData)
                else:
                    self._writeNew(filename=fullFileName,data=procData)

        else:
            fullFileName = os.path.join(caseDirectory, str(location), fileName)
            self.logger.debug(f"Writinge {fullFileName} as composed file")

            if os.path.exists(fullFileName):
                self.logger.debug("File exists, updating")
                self._updateExisting(filename=fullFileName, data=data)
            else:
                self.logger.debug("File does not exist, write new")
                self._writeNew(filename=fullFileName, data=data)

    def _updateExisting(self, filename, data):
        raise NotImplementedError("Implemented specifically for field or list")

    def _writeNew(self, filename, data):
        raise NotImplementedError("Implemented specifically for field or list")


class OFField(OFObject):
    """
        Holds a field.

        This object has:
            - dimensions
            - internalField
            - boundaryField


    There is a bug in updating an OF field -> removes the boundary conditions.

    """
    _dimensions = None  # The dimensions str (currently just a string without the parsing.

    @property
    def dimensions(self):
        return self._dimensions


    def __init__(self,name,dimensions,componentNames=None):
        """
            Initializes a field.

        Parameters
        ----------
        name :  str
            The name of the field.

        dimensions : str
            The dimensions of the field.

        cloudName : str, or None
            If None, then the field is eulerian. Otherwise,
            use this property as the field name.

        componentNames: list
            The name of the components. If None, then the field is a scalar.

        """
        super().__init__(name=name,componentNames=componentNames)
        self._dimensions = dimensions

    def load(self, caseDirectory, times=None, parallelCase=True):
        """
            Extracts a field to the disk from the requested times.
            If None, then reads from all the time steps.
            Reads only the internal field (and not the boundaries).

            If read in parallel, then add the processorto the output columns.


        Parameters
        -----------

            times: list
                 a list of time steps to extract.
                 If None, it extracts all time steps in the casePath.
            file: str
                The name for a file, in which the data is saved. Default is "<cloud name>_data.parquet", in the current working directory.

            withVelocity: bool
                    True: extract the particles' velocities in addition to their positions.
                    Default is False.

            **kwargs:
                Any additional parameters to add to the description in the DB.
        Returns
        --------
            document of the data.
        """
        finalCasePath = os.path.abspath(caseDirectory)

        if parallelCase:
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

            data = dask.dataframe.from_delayed(
                [delayed(extractFieldFile)(os.path.join(finalCasePath,processorName, str(timeName),self.name),
                                           columnNames=self.componentNames,
                                           time=timeName,
                                           processor=int(processorName[9:])) for processorName, timeName in product(processorList, timeList)])
        else:

            if times is None:
                timeList = sorted([x for x in os.listdir(finalCasePath) if (
                        os.path.isdir(x) and
                        x.isdigit() and
                        (not x.startswith("processor") and x not in ["constant", "system", "rootCase", 'VTK']))],
                                  key=lambda x: int(x))
            else:
                timeList = numpy.atleast_1d(times)

            data = dask.dataframe.from_delayed([delayed(extractFieldFile)(os.path.join(finalCasePath,str(timeName),self.name),
                                                                          columnNames = self.componentNames,
                                                                          time=timeName)
                                                for timeName in timeList])

        return data


    def _updateExisting(self,filename,data):
        """
            Injects the data into existing file.

        Parameters
        ----------
        filename: str
            The full path to the file to write.
        data: pandas
            The data to write (we assume it is in the right order.

        Returns
        -------
            None
        """
        if isinstance(data,(pandas.DataFrame,dask.dataframe.DataFrame)):
            columnNames = [x for x in data.columns if (x!='processor' and x!='time') ] if self.componentNames is None else self.componentNames
        else:
            columnNames = None

        filecontentStr = ""

        with open(filename,'r') as inFile:
            for inline in inFile:
                firstToken = inline.split(" ")[0].strip()

                if firstToken != "internalField":
                    filecontentStr += inline
                else:

                    if columnNames is None:
                        filecontentStr += f"internalField   uniform {data};\n"
                    else:
                        filecontentStr += "internalField   nonuniform List<vector>\n"

                        if len(columnNames) > 1:
                            # vector
                            filecontentStr += self.pandasToFoamFormat(data)
                        else:
                            # scalar
                            if isinstance(data, pandas.Series):
                                data = data.values
                            fileStrContent += ""
                            fileStrContent = f"{str(data.shape[0])}\n"
                            fileStrContent += "(\n"
                            fileStrContent += "\n".join(data)
                            fileStrContent += ");\n"

                        filecontentStr += "\n"
                        # Skip internal field until the boundaryField.
                        firstTokenInner = ""
                        while firstTokenInner != 'boundaryField':
                            inline = inFile.readline()
                            firstTokenInner = inline.split(" ")[0].strip()
                            if not inline:
                                break

            with open(filename, 'w') as outFile:
                outFile.write(filecontentStr)

            return filecontentStr

    def _writeNew(self,filename,data):
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
        self.logger.info("----- Start -----")

        fileStrContent = self._getHeader()
        fileStrContent += "\n\n" + f"dimensions {self.dimensions};\n\n"

        if type(data) in [float,int,list,tuple,str]: #isinstance(data,float) or isinstance(data,int):
            if type(data) in [float,int,str]:
                fileStrContent += f"internalField  uniform {data};\n"
            else:
                fileStrContent += f"internalField  uniform ({' '.join([str(x) for x in data])});\n"
        else:
            fileStrContent += "internalField   nonuniform List<vector>\n"

            if isinstance(data,pandas.Series):
                componentNames = ['demo']
            else:
                componentNames = [x for x in data.columns if
                               (x != 'processor' and x != 'time')] if self.componentNames is None else self.componentNames

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

        fileStrContent += """
    boundaryField
    {
        ".*"
        {
            type            zeroGradient;
        }
    }"""

        with open(filename,'w') as outfile:
            outfile.write(fileStrContent)

        return fileStrContent


    def emptyParallelField(self, caseDirectory,timeName="0.parallel",processor="", data=None,boundaryField=None):
        """
            Writes a null file with definitions of the processor boundries
            because sometimes, the snappyHex mesh of decomposed pars breaks the parallel
            fields. This file is used to correct this problem.

        Parameters
        ----------
        filename: str
            The name of the file to write to.

        data: float, str, pandas.DataFrame, dask.DataFrame
            The data to write to the field.

        Returns
        -------

        """
        if data is None:
            data = '0' if self.componentNames is None else f"({' '.join(['0' for x in self.componentNames])})"

        fileStrContent = self._getHeader()
        fileStrContent += "\n\n" + f"dimensions {self.dimensions};\n\n"

        if type(data) in [float,int,list,tuple,str]: #isinstance(data,float) or isinstance(data,int):
            if type(data) in [float,int,str]:
                fileStrContent += f"internalField  uniform {data};\n"
            else:
                fileStrContent += f"internalField  uniform ({' '.join([str(x) for x in data])});\n"
        else:
            fileStrContent += "internalField   nonuniform List<vector>\n"

            if isinstance(data,pandas.Series):
                componentNames = ['demo']
            else:
                componentNames = [x for x in data.columns if
                               (x != 'processor' and x != 'time')] if self.componentNames is None else self.componentNames

            if len(componentNames) > 1:
                # vector/tensor
                fileStrContent += self.pandasToFoamFormat(data,componentNames)
            else:
                # scalar
                if isinstance(data, pandas.Series):
                    data = data.values
                fileStrContent += ""
                fileStrContent = f"{str(data.shape[0])}\n"
                fileStrContent += "(\n"
                fileStrContent += "\n".join(data)
                fileStrContent += ");\n"

        boundaryConditions = ""
        if boundaryField is not None:
            for boundaryPatchName,boundaryData in boundaryField.items():
                bstr = f"\n{boundaryPatchName}\n"
                bstr += "{\n"
                for bcondProp,bcondData in boundaryData.items():
                    bstr += f"\t{bcondProp} {bcondData};\n"
                bstr += "}\n"
                boundaryConditions += bstr

        fileStrContent += """
   boundaryField
{
    "proc.*"
    {
        type            processor;
    }
    """  + boundaryConditions + """
}
"""
        filename = os.path.join(caseDirectory,processor, timeName, self.name)
        self.logger.debug(f"Saving the file {self.name} to {filename} ")
        with open(filename,'w') as outfile:
            outfile.write(fileStrContent)


class OFList(OFObject):
    """
        Just data.
    """

    def _updateExisting(self,filename,data):
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
        return self._writeNew(filename,data)

    def _writeNew(self,filename,data):
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
        if isinstance(data,pandas.Series):
            columnNames = ['demo']
        else:
            columnNames = [x for x in data.columns if
                           (x != 'processor' and x != 'time')] if self.columnNames is None else self.columnNames

        fileStrContent = self.getHeader()
        if len(columnNames) > 1:
            # vector
            fileStrContent += self.pandasToFoamFormat(data,columnNames)

        else:
            # scalar
            if isinstance(data, pandas.Series):
                fileStrContent +=  "\n".join(data)
            else:
                fileStrContent += "\n".join(data[columnNames])

        with open(filename,'w') as outfile:
            outfile.write(fileStrContent)

        return fileStrContent


#########################################################################
#
#
#########################################################################

def extractFieldFile(path, columnNames, **kwargs):
    """
        Extract data from openFOAM field file.

        The difference brom list file is that field files has boundary conditions as
        well, and so we must read from the file their length.


    :param path:
    :param columnNames:
    **kwargs:
        Adds this result to the pandas as column.
    :return:
        pandas.
    """

    # 1. find the number of points in the file.
    L = []
    with open(path, "r") as thefile:
        lines = thefile.readline()
        try:
            lineCount = int(lines)
        except ValueError:
            lineCount = None

        while lineCount is None:
            lines = thefile.readline()
            try:
                lineCount = int(lines)
            except ValueError:
                lineCount = None

        for line in islice(thefile, 1, lineCount+1):
          L.append(line.replace("(","").replace(")",""))

    return pandas.read_csv(StringIO("".join(L)),names=columnNames,sep=" ").assign(**kwargs)

def extractFile(path, columnNames, vector=True, skiphead = 20,skipend = 4):
    """
        Extracts data from a openFOAM list file.

        list files has no boundary and so, we can just skip head and end.

    Parameters
    ----------
    path: str
        The path of the file
    time: str
        The files' time step
    columnNames: list of str
        The names of the columns
    skiphead: int
        Number of lines to skip from the beginning of the file
    skipend: int
        Number of lines to skip from the ending of the file

    Returns
    -------
        Pandas with the data.
    """

    cnvrt = lambda x: float(x.replace("(", "").replace(")", ""))
    cnvrtDict = dict([(x, cnvrt) for x in columnNames])

    try:
        newData = pandas.read_csv(path,
                                  skiprows=skiphead,
                                  skipfooter=skipend,
                                  engine='python',
                                  header=None,
                                  delim_whitespace=True,
                                  converters=cnvrtDict,
                                  names=columnNames)
    except ValueError:
        newData = []

    if len(newData) == 0:
        with open(path, "r") as thefile:
            lines = thefile.readlines()

        vals = lines[17]
        data = []

        if vector:
            if "{" in vals:
                inputs = vals.split("{")
                repeat = int(inputs[0])
                valuesList = inputs[1][inputs[1].find("(") + 1:inputs[1].find(")")]
                data = dict(
                    [(colname, [float(x)] * repeat) for colname, x in zip(columnNames, valuesList.split(" "))])
            else:
                for rcrdListTuple in vals.split("(")[2:]:
                    record = dict(
                        [(name, float(y)) for name, y in zip(columnNames, rcrdListTuple.split(")")[0].split(" "))])
                    data.append(record)
        else:

            if "{" in vals:
                inputs = vals.split("{")
                repeat = int(inputs[0])
                value = float(inputs[1].split("}")[0])
                data = [{columnNames[0]: value} for x in range(repeat)]

            else:
                valuesList = vals.split("(")[1]
                for rcrdListItem in valuesList.split(" "):
                    record = {columnNames[0]: float(rcrdListItem.split(")")[0])}
                    data.append(record)

        newData = pandas.DataFrame(data)

    return newData.astype(float)



def readLagrangianRecord(timeName, casePath, withVelocity=False, withReleaseTimes=False, withMass=False,
                         cloudName="kinematicCloud"):

    print(f"Processing {timeName}")

    columnsDict = dict(x=[], y=[], z=[], id=[], procId=[], globalID=[], globalX=[], globalY=[], globalZ=[])
    if withMass:
        columnsDict['mass'] = []
    if withReleaseTimes:
        columnsDict['age'] = []
    if withVelocity:
        columnsDict['U_x'] = []
        columnsDict['U_y'] = []
        columnsDict['U_z'] = []

    newData = pandas.DataFrame(columnsDict, dtype=numpy.float64)

    try:
        newData = extractFile(
            os.path.join(casePath, timeName, "lagrangian", cloudName, "globalSigmaPositions"),
            ['x', 'y', 'z'])
        for fld in ['x', 'y', 'z']:
            newData[fld] = newData[fld].astype(numpy.float64)

        dataID = extractFile(os.path.join(casePath, timeName, "lagrangian", cloudName, "origId"),
                             ['id'], vector=False)
        newData['id'] = dataID['id'].astype(numpy.float64)

        dataprocID = extractFile(
            os.path.join(casePath, timeName, "lagrangian", cloudName, "origProcId"), ['procId'],
            vector=False)
        newData['procId'] = dataprocID['procId'].astype(numpy.float64)

        newData = newData.ffill().assign(globalID=1000000000 * newData.procId + newData.id)

        dataGlobal = extractFile(
            os.path.join(casePath, timeName, "lagrangian", cloudName, "globalPositions"),
            ['globalX', 'globalY', 'globalZ'])

        for col in ['globalX', 'globalY', 'globalZ']:
            newData[col] = dataGlobal[col].astype(numpy.float64)

        if withVelocity:
            dataU = extractFile(os.path.join(casePath, timeName, "lagrangian", cloudName, "U"),
                                ['U_x', 'U_y', 'U_z'])
            for col in ['U_x', 'U_y', 'U_z']:
                newData[col] = dataU[col]

        if withReleaseTimes:
            dataM = extractFile(os.path.join(casePath, timeName, "lagrangian", cloudName, "age"),
                                ['age'], vector=False)
            # newData["releaseTime"] = dataM["time"] - dataM["age"] + releaseTime

        if withMass:
            dataM = extractFile(os.path.join(casePath, timeName, "lagrangian", cloudName, "mass"),
                                ['mass'], vector=False)
            try:
                newData["mass"] = dataM["mass"]
            except:
                newData = newData.compute()
                newData["mass"] = dataM["mass"]

    except:
        pass

    theTime = os.path.split(timeName)[-1]
    newData['time'] = float(theTime)
    return newData



