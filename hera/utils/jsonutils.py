from unum.units import *
import os
import json
from json.decoder import JSONDecodeError

def ConvertJSONtoConf(JSON):
    """
        Traverse the JSON and replace all the unum values with objects.

    :param JSON:
    :return:
    """
    ret ={}
    for key,value in JSON.items():
        print(key,value)
        if isinstance(value,dict):
            ret[key] = ConvertJSONtoConf(JSON[key])
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

