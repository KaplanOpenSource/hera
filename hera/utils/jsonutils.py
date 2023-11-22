from unum.units import *
import os
import json
import pandas
from json.decoder import JSONDecodeError
from .unum import *
from unum import Unum
import logging
from .logging import get_classMethod_logger
from jsonpath_ng import parse
from itertools import product
from .dataframeutils import compareDataframeConfigurations

def compareJSONS(longFormat=False,**kwargs):
    """
        Recieves a group of name->JSONs (as file, string or dict)
        and returns the pandas that compares them.
    Parameters
    ----------
    kwargs

    longFormat : bool
        Return the value as long or wide

    Returns
    -------

    """
    fulldata = pandas.concat([convertJSONtoPandas(data).assign(datasetName=name) for name,data in kwargs.items()])
    return compareDataframeConfigurations(fulldata,datasetName="datasetName",parameterName="parameterNameFullPath",longFormat=longFormat)


def ConfigurationToJSON(valueToProcess):
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
    logger = logging.getLogger("hera.bin.ConfigurationToJSON")
    ret = {}
    logger.debug(f"Processing {valueToProcess}")
    if isinstance(valueToProcess,dict):
        for key, value in valueToProcess.items():
            logger.debug("\t dictionary, calling recursively")
            ret[key] = ConfigurationToJSON(value)
    elif isinstance(valueToProcess,list):
        logger.debug("\t list, transforming to str every item")
        ret = [ConfigurationToJSON(x) for x in valueToProcess]
    elif "'" in str(valueToProcess):
        logger.debug(f"\t {valueToProcess} is String, use as is")
        ret = valueToProcess
    else:
        logger.debug(f"\t Try to convert *{valueToProcess}* to string")
        ret = unumToStr(valueToProcess)

    return ret

def JSONToConfiguration(valueToProcess):
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
    logger = logging.getLogger("hera.util.JSONToConfiguration")
    ret ={}
    logger.debug(f"Processing {valueToProcess} of type {type(valueToProcess)}")
    if isinstance(valueToProcess,dict):
        for key, value in valueToProcess.items():
            logger.debug(f"Transforming key: {key}")
            ret[key] = JSONToConfiguration(value)
    elif isinstance(valueToProcess,list):
        logger.debug("\t list, transforming to unum every item")
        ret = [JSONToConfiguration(x) for x in valueToProcess]
    elif isinstance(valueToProcess,Unum):
        ret = valueToProcess
    elif isinstance(valueToProcess,int):
        ret = valueToProcess
    elif isinstance(valueToProcess,float):
        ret = valueToProcess
    elif "'" in str(valueToProcess):
        logger.debug(f"\t {valueToProcess} is a String, use as is")
        ret = valueToProcess
    else:
        logger.debug(f"\t Try to convert {valueToProcess} to unum")
        ret = strToUnum(valueToProcess)
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
    elif isinstance(jsonData, list):
        loadedjson = jsonData
    else:
        err = f"workflow type: {type(jsonData)} is unknown. Must be str, file-like or dict. "
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

def convertJSONtoPandas(jsonData, nameColumn="parameterNameFullPath", valueColumn="value"):
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
            newdata = newdata.assign(parameterNameFullPath=newdata.parameterName.apply(lambda x: f"{pname}.{x}"))
            base.append(newdata)

        pnds1 = pandas.concat(base,ignore_index=True)
        dictIndex = pnds1.apply(lambda x: isinstance(x.value, dict), axis=1)

    return pnds1

def JSONVariations(base,variationJSON):
    """
        Return a list of JSONs with all the variations.
    Parameters
    ----------
    base
    variationJSON

    Returns
    -------

    """
    variations = [[x for x in JSONvariationItem(dict(base),variation)] for variation in variationJSON]
    for var in product(*variations):
        params = {}
        local_base = dict(base)
        for item in var:
            params.update(item)

        for key, value in params.items():
            curkey = key if key.startswith("$.") else f"$.{key}"
            jsonpath_expr = parse(curkey)
            match = [match for match in jsonpath_expr.find(local_base)][0]
            jsonpath_expr.update(match, value)


        yield local_base

class JSONvariationItem:
    """
        An iterator that creates the variations of a single parameter block.

        Each group is a list of parameters.
        {
            p1 : list of values
            p2 :  list of values

        }

        where p1, p2 are json paths. (If $. is not specified in the path, add them).
    """

    base = None
    variationItem = None

    _curIter = None # the ID of the
    _itemCount = None

    def __init__(self,base,variationItem):
        """
            The base is the json that will be changes.

            variationItem is a map of json paths -> values. All the paths change together.
        Parameters
        ----------
        base
        variationItem
        """
        logger = get_classMethod_logger(self,name="init")

        self.base = base
        self._itemCount = None
        for key,value in variationItem.items():
            if self._itemCount is None:
                self._itemCount = len(value)
                logger.debug(f"Got {self._itemCount} items in key {key}. Now have to make sure that it equal to all keys ")
            else:
                if len(value) != self._itemCount:
                    err = f"The key {key} does not have the right number of values. Got {len(value)}, and should be {self._itemCount}"
                    logger.error(err)
                    raise ValueError(err)

        self.variationItem = variationItem
        self._curIter = 0


    def __iter__(self):
        return self

    def __next__(self):

        if self._curIter > self._itemCount-1:
            raise StopIteration

        result = {}
        for key, value in self.variationItem.items():
            result[key] = value[self._curIter]
        self._curIter += 1
        return result

