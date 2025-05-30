import os
import shutil
import numpy
from hera import Project
from hera.utils import dictToMongoQuery
import inspect
from functools import wraps

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

    def __init__(self, func,returnFormat,projectName = None,postProcessFunction=None,getDataParams={}):
        self.func = func
        self.postProcessFunction = postProcessFunction
        self.projectName = projectName
        self.getDataParams = getDataParams

    def __call__(self, *args, **kwargs):

        sig = inspect.signature(self.func)
        bound = sig.bind(*args, **kwargs)
        bound.apply_defaults()

        # Build the info dictionary
        call_info = {
            'functionName': self.func.__name__,
            'functionParameters': dict(bound.arguments)
        }

        self.checkIfFunctionIsCached(call_info)
        result = self.func(*args, **kwargs)
        self.saveFunctionCache(result)
        return result

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


    def saveFunctionCache(self,data):
        """
            Save the data to the disk.
        Parameters
        ----------
        data

        Returns
        -------

        """
        proj = Project(self.projectName)
        fileID = proj.
        qry = dictToMongoQuery(call_info)
        docList = proj.getCacheDocuments(type="functionCacheData", **qry)







