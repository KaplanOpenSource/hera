import pandas
import numpy
import os
import glob
from itertools import product
from ....utils import loadJSON
from ....utils.logging import get_classMethod_logger
from .. import FIELDTYPE_VECTOR, FIELDTYPE_TENSOR, FIELDTYPE_SCALAR, FIELDCOMPUTATION_EULERIAN, \
    FIELDCOMPUTATION_LAGRANGIAN,FLOWTYPE_INCOMPRESSIBLE
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile,WriteParameterFile
from PyFoam.Basics.DataStructures import Field,Vector,Tensor,DictProxy,Dimension
from .OFObject import OFObject



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

