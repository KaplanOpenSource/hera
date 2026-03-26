import os
import shutil
import numpy
from hera import Project
import inspect
from functools import wraps
from hera.datalayer.datahandler import getHandler,datatypes
from hera.utils import dictToMongoQuery,ConfigurationToJSON
import pickle
import base64
from bson import BSON
from bson.errors import InvalidDocument


def clearAllFunctionsCache(projectName=None):
    """
        Remove the cache of all functions.
    Parameters
    ----------
    projectName

    Returns
    -------

    """
    clearFunctionCache(functionName=None,projectName=projectName)

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
    paramDict = dict()
    if functionName is not None:
        paramDict['functionName'] = functionName

    docList = proj.deleteCacheDocuments(type ="functionCacheData",**paramDict)
    for doc in docList:
        if os.path.exists(doc['resource']):
            if os.path.isdir(doc['resource']):
                shutil.rmtree(doc['resource'])
            else:
                os.remove(doc['resource'])

    return True

def cacheFunction(_func=None, *, returnFormat=None, projectName=None, postProcessFunction=None, getDataParams={},storeDataParams={}):
    """
    Decorator that caches a function's return value in the project database.

    On first call, the function executes and its result is saved as a cache document.
    On subsequent calls with the same arguments, the cached result is returned
    instead of re-executing the function.

    Can be used with or without arguments::

        @cacheFunction
        def my_func(x):
            ...

        @cacheFunction(returnFormat=datatypes.PARQUET, projectName="myproject")
        def my_func(x):
            ...

    Parameters
    ----------
    returnFormat : str, optional
        The data format to use when storing the result. If None, auto-detected.
    projectName : str, optional
        The project to store the cache in. If None, loaded from caseConfiguration.
    postProcessFunction : callable, optional
        A function applied to the result before returning it.
    getDataParams : dict, optional
        Extra keyword arguments passed to ``getData`` when loading from cache.
    storeDataParams : dict, optional
        Extra keyword arguments passed when saving to cache.
    """
    def decorator(func):
        """Wrap the target function with caching logic."""
        @wraps(func)
        def wrapper(*args, **kwargs):
            """Invoke the cached version of the wrapped function."""
            return cacheDecorators(
                func=func,
                dataFormat=returnFormat,
                projectName=projectName,
                postProcessFunction=postProcessFunction,
                getDataParams=getDataParams,
                storeDataParams=storeDataParams
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
    """
    Internal implementation of the function caching mechanism.

    Wraps a function call, serializes its arguments, checks the database for
    a cached result, and stores the result if not already cached.
    """

    prepareFunction = None
    postProcessFunction = None
    projectName = None
    dataFormat = None

    @staticmethod
    def is_mongo_serializable(obj):
        """Check whether ``obj`` can be BSON-encoded for MongoDB storage."""
        try:
            # BSON expects a dict at the top level
            BSON.encode({'test': obj})
            return True
        except InvalidDocument:
            return False

    # Serialize an object into a plain text
    @staticmethod
    def obj_to_txt(obj):
        """Serialize ``obj`` to a base64-encoded text string via pickle."""
        message_bytes = pickle.dumps(obj)
        base64_bytes = base64.b64encode(message_bytes)
        txt = base64_bytes.decode('ascii')
        return txt

    # De-serialize an object from a plain text
    @staticmethod
    def txt_to_obj(txt):
        """Deserialize an object from a base64-encoded text string."""
        base64_bytes = txt.encode('ascii')
        message_bytes = base64.b64decode(base64_bytes)
        obj = pickle.loads(message_bytes)
        return obj

    def __init__(self, func,dataFormat,projectName = None,postProcessFunction=None,getDataParams={},storeDataParams={}):
        """
        Parameters
        ----------
        func : callable
            The function whose results are cached.
        dataFormat : str or None
            Storage format for the cached data.
        projectName : str or None
            Project that owns the cache collection.
        postProcessFunction : callable or None
            Optional transform applied to the result before returning.
        getDataParams : dict
            Extra keyword arguments forwarded to ``getData``.
        storeDataParams : dict
            Extra keyword arguments forwarded when saving data.
        """
        self.func = func
        self.postProcessFunction = postProcessFunction
        self.projectName = projectName
        self.getDataParams = getDataParams
        self.storeDataParams = storeDataParams
        self.dataFormat = dataFormat

    def __call__(self, *args, **kwargs):
        """Execute the function, returning a cached result when available."""
        sig = inspect.signature(self.func)

        # Bind the passed args and kwargs to the signature
        bound = sig.bind(*args, **kwargs)

        # Apply defaults to fill in any missing optional arguments
        bound.apply_defaults()

        call_info = dict(bound.arguments)

        if 'self' in call_info:
            call_info['context'] = call_info.pop('self')

        # convert any pint/unum to standardized MKS and dict with the magnitude and units seperated.
        # This will allow the query of the querys even if they are given in different units
        call_info_JSON = ConfigurationToJSON(call_info, standardize=True, splitUnits=True, keepOriginalUnits=True)

        call_info_serialized = dict()

        for key,value in call_info_JSON.items():
            serializable = cacheDecorators.is_mongo_serializable(value)
            serialized_value = value if serializable else cacheDecorators.obj_to_txt(value)
            call_info_serialized[key] = (serializable,serialized_value)


        # Add the function name

        call_info_serialized['functionName'] =  self._get_full_func_name(self.func)

        data = self.checkIfFunctionIsCached(call_info_serialized)
        if data is None:
            data = self.func(*args, **kwargs)
            # query without the original units. This allows the user to query different units
            #call_info_query_JSON = ConfigurationToJSON(call_info, standardize=True, splitUnits=True, keepOriginalUnits=False)
            if data is not None:
                doc = self.saveFunctionCache(call_info_serialized,data)

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
        return None if len(docList)==0 else docList[0].getData(**self.getDataParams)

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
        return proj.saveCacheData(name=call_info['functionName'], data=data, desc=call_info, type="functionCacheData",dataFormat=self.dataFormat)






