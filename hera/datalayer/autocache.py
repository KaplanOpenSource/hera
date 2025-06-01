import os
import shutil
import numpy
from hera import Project
import inspect
from functools import wraps
from hera.datalayer.datahandler import getHandler,datatypes
from hera.utils import dictToMongoQuery,ConfigurationToJSON


def clearFunctionCache(functionName,projectName=None):
    """
        Removes all the cache documents of the function with the data from the disk.
    Parameters
    ----------
    functionName : str
        The name of the function
    projectName : str
        The name of the project that holds the cache. If None, load the name from the caseConfiguration.

    Returns
    -------

    """
    proj = Project(projectName=projectName)
    docList = proj.deleteCacheDocuments(type ="functionCacheData",functionName=functionName)
    for doc in docList:
        if os.path.exists(doc['resource']):
            if os.path.isdir(doc['resource']):
                shutil.rmtree(doc['resource'])
            else:
                os.remove(doc['resource'])

    return True

def cacheFunction(_func=None, *, returnFormat=None, projectName=None, postProcessFunction=None, getDataParams={}):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            return cacheDecorators(
                func=func,
                dataFormat=returnFormat,
                projectName=projectName,
                postProcessFunction=postProcessFunction,
                getDataParams=getDataParams
            )(*args, **kwargs)
        return wrapper

    if _func is None:
        # Decorator called with parentheses and arguments
        return decorator
    else:
        # Decorator used directly like @cacheFunction
        return decorator(_func)
##############################################################################################
#####                                                                                    #####
#####                       Main cacheDecorators                                         #####
#####                                                                                    #####
##############################################################################################
class cacheDecorators:

    prepareFunction = None
    postProcessFunction = None
    projectName = None
    dataFormat = None

    def __init__(self, func,dataFormat,projectName = None,postProcessFunction=None,getDataParams={}):
        self.func = func
        self.postProcessFunction = postProcessFunction
        self.projectName = projectName
        self.getDataParams = getDataParams
        self.dataFormat = dataFormat

    def __call__(self, *args, **kwargs):


        sig = inspect.signature(self.func)

        # Bind the passed args and kwargs to the signature
        bound = sig.bind(*args, **kwargs)

        # Apply defaults to fill in any missing optional arguments
        bound.apply_defaults()

        call_info = dict(bound.arguments)

        # Combine them
        call_info['functionName'] =  self._get_full_func_name(self.func)

        # convert any pint/unum to standardized MKS and dict with the magnitude and units seperated.
        # This will allow the query of the querys even if they are given in different units
        call_info_JSON = ConfigurationToJSON(call_info, standardize=True, splitUnits=True, keepOriginalUnits=True)

        data = self.checkIfFunctionIsCached(call_info_JSON)
        if data is None:
            data = self.func(*args, **kwargs)
            # query without the original units. This allows the user to query different units
            #call_info_query_JSON = ConfigurationToJSON(call_info, standardize=True, splitUnits=True, keepOriginalUnits=False)
            doc = self.saveFunctionCache(call_info,data)

        ret = data if self.postProcessFunction is None else self.postProcessFunction(data)
        return ret

    def _get_full_func_name(self,func):
        """Returns the full qualified path: module.[class.]function_name"""
        if not callable(func):
            raise TypeError("Provided object is not callable.")

        # Handle bound methods by unwrapping them
        if inspect.ismethod(func):
            # Get the original function and its class
            cls = func.__self__.__class__
            method_name = func.__name__
            class_qualname = cls.__qualname__
            module = func.__module__
            if module == "__main__":
                ret = f"{class_qualname}.{method_name}"
            else:
                ret = f"{module}.{class_qualname}.{method_name}"

        elif inspect.isfunction(func):
            ret = func.__qualname__
        else:

            # Handle unbound class or static methods and plain functions
            qualname = func.__qualname__
            module = func.__module__
            ret = f"{module}.{qualname}"
        return ret


    def checkIfFunctionIsCached(self,call_info):
        """
            Check if the function and the parameters are stored in the DB.
        Parameters
        ----------
        call_info : dict
            A dict with the info on the function that was called.
            functionName and functionParameters as parameters.

        Returns
        -------
            None if the data does not exist,
            the data otherwise.
        """

        proj = Project(self.projectName)
        docList = proj.getCacheDocuments(type="functionCacheData",**call_info)
        return None if len(docList)==0 else docList[0].getData()

    def saveFunctionCache(self,call_info,data):
        """
            Save the data to the disk.
        Parameters
        ----------
        data

        Returns
        -------

        """
        proj = Project(self.projectName)
        qry  = ConfigurationToJSON(call_info,standardize=True,splitUnits=True,keepOriginalUnits=True)

        guessedDataFormat = datatypes.getDataFormatName(data) if self.dataFormat is None else self.dataFormat
        handler = datatypes.getHandler(guessedDataFormat)
        file_extension = datatypes.getDataFormatExtension(data)

        cacheDirectory = os.path.join(proj.filesDirectory,"cache")
        os.makedirs(cacheDirectory,exist_ok=True)
        fileID = proj.getCounterAndAdd(call_info['functionName'])
        fileName = os.path.join(cacheDirectory,f"{call_info['functionName']}_{fileID}.{file_extension}")

        handler.saveData(data,fileName)

        doc = proj.addCacheDocument(type="functionCacheData",dataFormat = guessedDataFormat,resource=fileName  ,desc = qry)
        return doc







