import pandas
import numpy

def centersToPandas(filepath='C', skiphead = 22,skipend = 8, saveToTxt=False, fileName="CellCenters.txt"):
    """
        Extract pandas from openfoam cell centers file.

        It is also possible to save it to  a new txt file

    Parameters
    -----------

    filepath: str
        The path to the cell centers file.
        Default: 'C'

    saveToTxt:  boolean
        Weather to save the centers coords to txt file
        Default: False

    fileName: str
        The file name to save.
        default: None

    return
    ------

    cellData: pandas DF

    """

    cellData = pandas.read_csv(filepath, skiprows=skiphead,
                               skipfooter=skipend,
                               engine='python',
                               header=None,
                               delim_whitespace=True, names=['x', 'y', 'z'])

    cellData['x'] = cellData['x'].str[1:]
    cellData['z'] = cellData['z'].str[:-1]
    cellData = cellData.astype(float)

    if saveToTxt:

        L = numpy.array(cellData[['x', 'y', 'z']])
        numberOfCells = L.shape[0]
        out=[]
        for i in range(numberOfCells):
            out.append(f"({L[i, 0]} {L[i, 1]} {L[i, 2]})\n")

        with open(fileName, "w") as outfile:
            outfile.write("".join(out))

    return cellData