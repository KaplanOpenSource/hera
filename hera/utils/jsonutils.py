from unum.units import *
import os
import json

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

def loadJSON(json):
    """
        Reads the json object to the memory.

        Could be:

            * file object: any file-like object with the property 'read'.
            * str: either the JSON or a path to the directory.
            * dict: the JSON object.

    Parameters
    ----------
    json : str, object file, path to disk, dict
        The object that contains the dict.

    Returns
    -------
        dict
        The loaded JSON.

    """

    if hasattr(json, 'read'):
        loadedjson = json.load(json)
    elif isinstance(json, str):

        if os.path.exists(json):
            with open(json) as jsonFile:
                loadedjson = json.load(jsonFile)
        else:
            loadedjson = json.loads(json)
    elif isinstance(json, dict):
        loadedjson = json
    else:
        err = f"workflow type: {type(json)} is unknonw. Must be str, file-like or dict. "
        raise ValueError(err)


    return  loadedjson