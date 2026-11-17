import numpy
import pandas
from foamlib import FoamFieldFile


#########################################################################
#
#               Utils
#########################################################################
def readFieldFile(fieldPath):
    """Parse an OpenFOAM field file into a dict, without its FoamFile header."""
    return FoamFieldFile(fieldPath).as_dict()


def extractFieldFile(casePath, columnNames, patchNameList=None,filterInternalPatches=True, **kwargs):
    """Parse an OpenFOAM field file and return its data as a DataFrame."""
    try:
        dataParsedFile = readFieldFile(casePath)
    except Exception as e:
        raise ValueError(f"Cannot read the OpenFOAM field file {casePath}: {e}") from e
    return ParsedParameterFileToDataFrame(dataParsedFile=dataParsedFile,columnNames=columnNames, patchNameList=patchNameList,filterInternalPatches=filterInternalPatches, **kwargs)

def _rows(value, columnNames):
    """One row per cell/face: a uniform value is one row, a list is one row per item."""
    return numpy.reshape(value, (-1, len(columnNames)))


def ParsedParameterFileToDataFrame(dataParsedFile,columnNames, patchNameList=None,filterInternalPatches=True, **kwargs):
    """Convert a parsed OpenFOAM field (the dict that readFieldFile returns) to a pandas DataFrame."""
    ret = []
    pndsData = pandas.DataFrame(_rows(dataParsedFile['internalField'], columnNames),
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
            if 'value' in dataParsedFile['boundaryField'][patchName]:
                if numpy.size(dataParsedFile['boundaryField'][patchName]['value']) >0:
                    pndsData = pandas.DataFrame(
                        _rows(dataParsedFile['boundaryField'][patchName]['value'], columnNames),
                        columns=columnNames).assign(**kwargs, region='boundaryField', boundary=patchName,type=dataParsedFile['boundaryField'][patchName]['type'])
                else:
                    continue # skip that boundary
            else:
                pndsData = pandas.DataFrame([[numpy.nan]*len(columnNames)],columns=columnNames).assign(**kwargs, region='boundaryField', boundary=patchName,type=dataParsedFile['boundaryField'][patchName]['type'])

            ret.append(pndsData)

    return pandas.concat(ret).reset_index().rename(columns=dict(index="processorIndex"))
