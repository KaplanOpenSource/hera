from unum.units import *
import os
import json
import pandas
from json.decoder import JSONDecodeError

def convertJSONtoConf(JSON):
    """
        Traverse the JSON and replace all the unum values with objects.

    :param JSON:
    :return:
    """
    ret ={}
    for key,value in JSON.items():
        print(key,value)
        if isinstance(value,dict):
            ret[key] = convertJSONtoConf(JSON[key])
        elif isinstance(value,list):
            ret[key] = value
        else:

            if "'" in str(value):
                ret[key] = value
            else:
                try:
                    ret[key] = eval(str(value))
                except:
                    ret[key] = value

    return ret

def loadJSON(jsonData):
    """
        Reads the json object to the memory.

        Could be:

            * file object: any file-like object with the property 'read'.
            * str: either the JSON or a path to the directory.
            * dict: the JSON object.

    Parameters
    ----------
    jsonData : str, object file, path to disk, dict
        The object that contains the dict.

    Returns
    -------
        dict
        The loaded JSON.

    """

    if hasattr(jsonData, 'read'):
        loadedjson = json.load(jsonData)
    elif isinstance(jsonData, str):
        if os.path.exists(jsonData):
            with open(jsonData) as jsonFile:
                loadedjson = json.load(jsonFile)
        else:
            try:
                loadedjson = json.loads(jsonData)
            except JSONDecodeError as e:
                raise ValueError(f" {str(e)}: Got {jsonData}, either bad format of file does not exist")

    elif isinstance(jsonData, dict):
        loadedjson = jsonData
    else:
        err = f"workflow type: {type(jsonData)} is unknonw. Must be str, file-like or dict. "
        raise ValueError(err)


    return  loadedjson


def processJSONToPandas(jsonData, nameColumn="parameterName", valueColumn="value"):
    """
        Trnasforms a JSON to pandas, flattens a list items and names them according to their order

        The example the JSON :
        {
          "nodes": {
              "a" : {
                "x" : 1,
                "y" : 2,
                "z" : 3
              },
               "b" : {
                "r" : 1,
                "t" : 2,
                "y" : [3,2,4,5,6]
              }
          }
        }

        will be converted to

        parameterName  value
        --------------------
    0     nodes.a.x     1
    1     nodes.a.y     2
    2     nodes.a.z     3
    3     nodes.b.r     1
    4     nodes.b.t     2
    5   nodes.b.y_0     3
    6   nodes.b.y_1     2
    7   nodes.b.y_2     4
    8   nodes.b.y_3     5
    9   nodes.b.y_4     6

        Notes:

            - Currently does not support JSON whose root is a list.
              [
                {  "a" : 1 },
                {  "b" : 2}

              ]
                It will be supported if needed in the future.

        Parameters
        ----------

        jsonData : dict
            the JSON data (a dict)

        nameColumn: str
            The name of the parameter column name

        valueColumn : str
            The name of the value

    """
    pnds = pandas.json_normalize(jsonData).T.reset_index().rename(columns={'index': nameColumn, 0: valueColumn})\
        .explode(valueColumn)\
        .reset_index()

    listParameters = pnds.groupby(nameColumn).count().query(f"{valueColumn}>1").index

    for pname in listParameters:
        counter = 0
        for I, dta in pnds.iterrows():
            if dta.loc[nameColumn] == pname:
                pnds.loc[I, nameColumn] = f"{pname}_{counter}"
                counter += 1

    return pnds[[nameColumn,valueColumn]]


def convertJSONtoPandas(jsonData, nameColumn="parameterName", valueColumn="value"):
    """
        converts a JSON (either in file or loaded, or json str) to pandas.
        The pandas flattens the JSON using the json path convection.
        e.g
        {
            "a" : {
                "b" : 1,
                "c" : [1,2,3]
            }
        }

        will be converted to
         a.b  1
         a.c_0 1
         a.c_1 2
         a.c_2 3


        Does not support (currently) JSON whose root is a list but only supports dict

    Parameters
    ----------
    jsonData : str,dict
        A json data either a file name, a json dict string, or a dict.


        nameColumn: str
            The name of the parameter column name

        valueColumn : str
            The name of the value

    Returns
    -------
            pandas.DataFrame

            with the fields nameColumn (the path of the json) and valueColumn
    """
    param1 =  loadJSON(jsonData)
    pnds1 = processJSONToPandas(param1,nameColumn=nameColumn,valueColumn=valueColumn)
    dictIndex = pnds1.apply(lambda x: isinstance(x.value,dict),axis=1)
    while dictIndex.sum()>0:
        base = [pnds1[~dictIndex]]

        toProcessList = pnds1[dictIndex].set_index("parameterName")[['value']]
        for pname,data in toProcessList.iterrows():
            newdata = processJSONToPandas(data.value,nameColumn=nameColumn,valueColumn=valueColumn)
            newdata = newdata.assign(parameterName=newdata.parameterName.apply(lambda x: f"{pname}.{x}"))
            base.append(newdata)

        pnds1 = pandas.concat(base,ignore_index=True)
        dictIndex = pnds1.apply(lambda x: isinstance(x.value, dict), axis=1)

    return pnds1
