from unum.units import *
import os
import json
import pandas
from json.decoder import JSONDecodeError
from .unum import *
from unum import Unum


def unumToStr(obj):
    if isinstance(obj, Unum):
        ret = str(obj).replace(" [", "*").replace("]", "")
    else:
        ret = str(obj)
    return ret

def strToUnum(value):
    try:
        ret = eval(str(value))
    except:
        ret = value
    return ret

def ConfigurationToJSON(conf):
    """
        Converts a configuration dict (that might include unum objects) to
        JSON dict (where all the values are strings).

        The unum objects are converted to Str in a way that allows for their
        retrieval. (see the JSONToConfiguration function)

    Parameters
    ----------
    conf : dict
        A key-value dict where the values are string.
        Converts unum objects to string like representations:
            1*m/s --> '1*m/s'.
    Returns
    -------
        dict with all the values as string

    """

    ret = {}
    for key,value in conf.items():
        if isinstance(value,dict):
            ret[key] = ConfigurationToJSON(value)
        elif isinstance(value,list):
            ret[key] = [unumToStr(x) for x in value]
        elif isinstance(value,Unum):
            ret[key] = unumToStr(value)
        else:
            ret[key] = value


    return ret

def JSONToConfiguration(JSON):
    """
        Converts a dictionary (all the values are string) to
        a JSON where all the values are string.

        Convert the JSON back to configuration object using  the ConfigurationToJSON function.

    Parameters
    ----------
        JSON : dict
        A key-value where all the values are strings.
        The unum objects has the format '1*<unit>' (for exampe '1*m/s')

    Returns
    -------
        A dict with the values convected to unum.

    """
    ret ={}
    for key,value in JSON.items():
        if isinstance(value,dict):
            ret[key] = JSONToConfiguration(JSON[key])
        elif isinstance(value,list):
            ret[key] = [strToUnum(x) for x in value]
        elif isinstance(value, Unum):
            ret[key] = unumToStr(value)
        else:
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
        .explode(valueColumn,ignore_index=True)\
        .reset_index()

    # Handles nested lists. keep on exploding until all is flat!.
    while True:
        listParameters = pnds.groupby(nameColumn).count().query(f"{valueColumn}>1").index
        for pname in listParameters:
            counter = 0
            for I, dta in pnds.iterrows():
                if dta.loc[nameColumn] == pname:
                    pnds.loc[I, nameColumn] = f"{pname}_{counter}"
                    counter += 1

        # Handles lists with 1 item.
        for I, dta in pnds.iterrows():
            if isinstance(dta[valueColumn],list):
                if len(dta[valueColumn]) ==1:
                    pnds.loc[I, nameColumn] = f"{pnds.loc[I][nameColumn]}_{0}"
                    pnds.at[I, valueColumn] = pnds.loc[I][valueColumn][0]

        # Handling nested lists.
        tmp = pnds.explode(valueColumn,ignore_index=True)
        if len(tmp) == len(pnds):
            break
        else:
            pnds = tmp



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



