import pandas
import numpy

def centersToPandas(filepath='C', saveToTxt=False):

    skiphead = 22
    skipend = 45043


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

        with open("CellCenters.txt", "w") as outfile:
            outfile.write("".join(out))

    return cellData