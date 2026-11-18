import numpy
import os
import tempfile
from foamlib import FoamFieldFile
from ....utils.logging import get_classMethod_logger
from .. import FIELDTYPE_VECTOR, FIELDTYPE_TENSOR, FIELDTYPE_SCALAR

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
        """Return the list of component column names based on field type."""
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
        """Return the processor names in this field."""
        return self.data.keys()

    @property
    def processorItems(self):
        """Return the processor name-data pairs."""
        return self.data.items()

    @property
    def dimensionsStr(self):
        """Return the dimensions as an OpenFOAM-formatted string."""
        return self.getDimensions(**self.dimensions)

    @property
    def dimensionsList(self):
        """Return the dimensions as an ordered list of exponents."""
        return [ self.dimensions.get('kg',0),
                 self.dimensions.get('m',0),
                 self.dimensions.get('s',0),
                 self.dimensions.get('K',0),
                 self.dimensions.get('mol',0),
                 self.dimensions.get('A',0),
                 self.dimensions.get('cd',0)]

    def boundaryField(self,processorName='singleProcessor'):
        """Return the boundary field data for a given processor."""
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
            directory = os.path.join(caseDirectory,str(timeOrLocation))
            os.makedirs(directory,exist_ok=True)
            with open(os.path.join(directory,self.fileName),'wb') as outputFile:
                outputFile.write(self.toFoamText(self.data['singleProcessor']))
        else:
            for procName,procData in self.data.items():
                directory = os.path.join(caseDirectory,procName,str(timeOrLocation))
                os.makedirs(directory, exist_ok=True)
                with open(os.path.join(directory,self.fileName),'wb') as outputFile:
                    outputFile.write(self.toFoamText(procData))

    def toFoamText(self, processorData):
        """
            Returns the OpenFOAM file content (bytes) of one processor's data.

        Parameters
        ----------
        processorData : dict
            dimensions, internalField and boundaryField, as foamlib represents them.

        Returns
        -------
            bytes
        """
        header = {"version": 2.0, "format": "ascii", "class": f"vol{self.fieldType.title()}Field", "object": self.fileName}
        content = {"FoamFile": header, **processorData}

        # ponytail: foamlib < 1.10 writes any 3, 6 or 9 long array as a uniform
        # vector/symmTensor/tensor, whatever the file class.  A scalar list of that
        # length (a patch with 3 faces) is swapped for a placeholder word and written
        # by hand.  Drop this once foamlib >= 1.10 (it needs numpy >= 2).
        lists = {}

        def protect(value):
            if isinstance(value, numpy.ndarray) and value.ndim == 1 and len(value) in (3, 6, 9):
                placeholder = f"heraScalarList{len(lists)}"
                lists[placeholder] = value
                return placeholder
            return value

        if self.fieldType == FIELDTYPE_SCALAR:
            content["internalField"] = protect(content.get("internalField"))
            content["boundaryField"] = {patch: {key: protect(value) for key, value in patchData.items()}
                                        for patch, patchData in content.get("boundaryField", {}).items()}

        # Written through a file, not FoamFieldFile.dumps: dumps puts the whole file on
        # one line, and line-based readers (the Hmix template of the LSM toolkit) need
        # OpenFOAM's own layout, one entry per line.
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, self.fileName)
            with FoamFieldFile(path) as foamFile:
                for key, value in content.items():
                    foamFile[key] = value
            with open(path, "rb") as written:
                text = written.read()
        for placeholder, value in lists.items():
            listText = f"nonuniform List<scalar> {len(value)}({' '.join(repr(float(x)) for x in value)})"
            text = text.replace(placeholder.encode(), listText.encode())
        return text

