import pandas
import numpy
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile,WriteParameterFile
from PyFoam.Basics.DataStructures import Field,Vector,Tensor,DictProxy,Dimension

#########################################################################
#
#               Utils
#########################################################################
def extractFieldFile(casePath, columnNames, patchNameList=None,filterInternalPatches=True, **kwargs):
    """Parse an OpenFOAM field file and return its data as a DataFrame."""
    try:
        dataParsedFile = ParsedParameterFile(casePath)
    except Exception as e:
        print(casePath)
        raise ValueError(e)
    return ParsedParameterFileToDataFrame(dataParsedFile=dataParsedFile,columnNames=columnNames, patchNameList=patchNameList,filterInternalPatches=filterInternalPatches, **kwargs)

def ParsedParameterFileToDataFrame(dataParsedFile,columnNames, patchNameList=None,filterInternalPatches=True, **kwargs):
    """Convert a PyFoam ParsedParameterFile to a pandas DataFrame."""
    ret = []
    pndsData = pandas.DataFrame([[x for x in item] for item in numpy.atleast_2d(dataParsedFile['internalField'].val)],
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
                if len(dataParsedFile['boundaryField'][patchName]['value'].val) >0:
                    pndsData = pandas.DataFrame(
                        [[x for x in item] for item in numpy.atleast_2d(dataParsedFile['boundaryField'][patchName]['value'].val)],
                        columns=columnNames).assign(**kwargs, region='boundaryField', boundary=patchName,type=dataParsedFile['boundaryField'][patchName]['type'])
                else:
                    continue # skip that boundary 
            else:
                pndsData = pandas.DataFrame([[numpy.nan]*len(columnNames)],columns=columnNames).assign(**kwargs, region='boundaryField', boundary=patchName,type=dataParsedFile['boundaryField'][patchName]['type'])

            ret.append(pndsData)

    return pandas.concat(ret).reset_index().rename(columns=dict(index="processorIndex"))
