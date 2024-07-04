import pandas

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
