import logging
import numpy
import pandas

def compareDataframeConfigurations(data,datasetName="datasetName",parameterName="parameterName",valueName="value",indexList=None,longFormat=False,changeDotToUnderscore=False):
    """
        Compares datasets and outputs the keys that differ between them.

        The datasets can be given as:
            - list of tuples : (name,dataframe)
            - list of dict   : { <setName> : name, "data":dataframe}, where setName is given as input
            - dataframe, with a column setName as the distinguising between the sets (i.e long format)

        The structure of the dataframe is:

            - <parameterName> (the path of the parameter with '.' as a separator).
                parameterName is an input

            - <value>     - The value of the parameter
                value is an input

            - <index1,index2,...>  - The name of sub directories in the dataset.
                the list is an input.



    Parameters
    ----------
    data : pandas.DataFrame, list
            The data to compare.
            - list of tuples : (name,dataframe)
            - list of dict   : { <setName> : dataframe}, where setName is given as input
            - dataframe, with a column setName as the distinguising between the sets (i.e long format)

    datasetName: str
        the name of the datasetName key (if configuration is a list of dicts) or column
        (if the configuration is a pandas.Dataframe). Ignored if configuration is a list of tuples.

    parameterName: str
        The column that contains the json path.

    valueName: str
        The name of the value column.

    indexList: list, None
        A list of columns that create sub dictionaries in the data.

    longFormat : bool
        The format to return.

    changeDotToUnderscore : bool
        If true, change the columns names from x.y -> x_y. This will allow using pandas.query function

    Returns
    -------
        pandas.dataframe
        with the parameters that differ between the different datasets.

        The dataframe has the structure:

        If longFormat is True:

                    index                              parameter name value
        --------------------------------------------
        <index 1> .... <index 2> <parameter name> |


        If longFormat is False:

                            index                       dataset 1 | dataset 2 | dataset 3
        --------------------------------------------
        <index 1> .... <index 2> <parameter name> |

    """
    logger = logging.getLogger("hera.utils.compareConfigurationsPandas")
    logger.info("--- Start ---")
    diffList = []
    indexList = [] if indexList is None else [x for x in numpy.atleast_1d(indexList)]


    # Convert the configuration to a single dataframe in longformat.
    if isinstance(data,list):
        if len(data) == 0:
            logger.debug("Empty input, returning {}")
            return pandas.Dataframe({})
        else:
            tmpList = []
            for ditem in data:
                if isinstance(ditem,tuple):
                    item_data = ditem[1]
                    item_name = ditem[0]

                elif isinstance(ditem,dict):
                        try:
                            item_name = ditem[datasetName]
                        except KeyError:
                            logger.error(f"{ditem} does not have the key for the dataset name = {datasetName}")
                            raise KeyError(f"{ditem} does not have the key for the dataset name = {datasetName}")
                        try:
                            item_data = ditem['data']
                        except KeyError:
                            logger.error(f"{ditem} does not have key 'data' ")
                            raise KeyError(f"{ditem} does not have key 'data' ")
                else:
                    logger.error("The data is in incorrect format. List must consits a tuple (name,data) or dict {<namekey>:dat}")
                    raise ValueError("The data is in incorrect format. List must consits a tuple (name,data) or dict {<namekey>:dat}")

                tmpList.append(item_name.assign(**{datasetName: item_name}))

            configurations = pandas.concat(tmpList)
    elif isinstance(data,pandas.DataFrame):
        configurations = data
    else:
        logger.error("data must be a list of tuples, a list of dict or a pandas.Dataframe")
        raise KeyError("data must be a list of tuples, a list of dict or a pandas.Dataframe")

    datasetCount = configurations[datasetName].unique().shape[0]

    for grpid, grpdata in configurations.groupby([parameterName]+indexList):
        logger.debug(f"Testing {grpid} --> {grpdata[valueName]} ")
        if grpdata[valueName].unique().shape[0] > 1:
            logger.debug(f"{grpid}:: Normal Field. Different  ")
            diffList.append(grpdata.copy())

        if 0 < grpdata[valueName].count() < datasetCount:
            diffList.append(grpdata.copy())

    if len(diffList) > 0:
        ret = pandas.concat(diffList)
        if longFormat is False:
            ret = ret.pivot(index=indexList+[parameterName], columns=datasetName, values=valueName)

    else:
        ret = data[[datasetName]].drop_duplicates()

    if changeDotToUnderscore:
        try:
            newColNames = [(oldName,oldName.replace(".","_")) for oldName in ret.T.columns]
        except AttributeError:
            import pdb
            pdb.set_trace()
        ret_tmp = ret.T.rename(columns=dict(newColNames))
        ret = ret_tmp.T
    return ret

