import os
import shutil
import numpy
from hera import Project
from hera.utils import dictToMongoQuery
import inspect
from functools import wraps
from hera.datalayer.datahandler import getHandler,datatypes

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
    docList = proj.get(type ="functionCacheData",functionName=functionName)
    for doc in docList:
        if os.path.exists(doc.resource):
            if os.path.isdir(doc.resource):
                shutil.rmtree(doc.resource)
            else:
                os.remove(doc.resource)

    [doc.delete() for doc in docList]
    return True

def cache(returnFormat,projectName=None,postProcessFunction=None,getDataParams={}):
    def wrapper(func):
        return cacheDecorators(func=func,
                               projectName=projectName,
                               returnFormat=returnFormat,
                               postProcessFunction=postProcessFunction,
                               getDataParams = getDataParams)
    return wrapper

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
        call_info = self._top_caller_info()
        data = self.checkIfFunctionIsCached(call_info)
        if data is None:
            data = self.func(*args, **kwargs)
            doc = self.saveFunctionCache(call_info,data)

        ret = data if self.postProcessFunction is None else self.postProcessFunction(data)
        return ret

    def _top_caller_info(self):
        """
            Since it is called from the __call__
            it gets all the data from the function 2 stacks up.
        Returns
        -------

        """
        # Get the frame of the caller function (1 level back)
        frame = inspect.currentframe()
        caller_frame = frame.f_back.f_back  # caller of this method

        # Extract the code object of caller
        code = caller_frame.f_code
        func_name = code.co_name

        # Try to get module name
        module = inspect.getmodule(caller_frame)
        module_name = module.__name__ if module else None

        # Try to get class name from 'self' or 'cls' in caller locals
        cls_name = None
        # Common names for instance and class methods
        for possible_self in ('self', 'cls'):
            if possible_self in caller_frame.f_locals:
                cls_name = caller_frame.f_locals[possible_self].__class__.__name__
                break

        # Build full function name string
        if cls_name:
            full_name = f"{cls_name}.{func_name}" if module_name is None else f"{module_name}.{cls_name}.{func_name}"
        else:
            full_name = func_name if module_name is None else f"{module_name}.{func_name}"

        # Get the parameters and their values of the caller function
        try:
            sig = inspect.signature(caller_frame.f_code)
            # Note: signature() works better with function objects,
            # here we use the frame to get arguments passed in the frame's locals
            # We'll try to map parameters to local variables by name
        except Exception:
            sig = None

        # Get argument names from code object
        argcount = code.co_argcount
        varnames = code.co_varnames[:argcount]

        # Get actual parameter values from caller's locals
        params = {name: caller_frame.f_locals.get(name, '<no value>') for name in varnames}

        if 'self' in params:
            del params['self']

        params['functionName'] = full_name
        return params

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
        qry = dictToMongoQuery(call_info)
        docList = proj.getCacheDocuments(type="functionCacheData",**qry)
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
        fileID = proj.getCounterAndAdd(call_info['functionName'])
        qry = dictToMongoQuery(call_info)

        if self.dataFormat is None:
             guessedDataFormat = datatypes.type_to_datatype_string[data]
        else:
            guessedDataFormat = self.dataFormat


        file_extension = datatypes.datatype_to_extension[guessedDataFormat]

        cacheDirectory = os.path.join(proj.filesDirectory,"cache")
        os.makedirs(cacheDirectory,exist_ok=True)
        fileName = os.path.join(cacheDirectory,f"{call_info['functionName']}_{fileID}.{file_extension}")

        handler = guessHandler(guessedDataFormat)
        handler.saveData(data,fileName)

        doc = proj.addCacheDocument(type="functionCacheData",format = guessedDataFormat,resource=fileName  ,desc = qry)
        return doc







