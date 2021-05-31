from unum.units import *
def ConvertJSONtoConf(JSON):
    """
        Traverse the JSON and replace all the unum values with objects.

    :param JSON:
    :return:
    """
    ret ={}
    for key,value in JSON.items():
        if isinstance(value,dict):
            ret[key] = ConvertJSONtoConf(JSON[key])
        elif isinstance(value,list):
            ret[key] = value
        else:
            try:
                ret[key] = eval(str(value))
            except:
                ret[key] = value

    return ret
