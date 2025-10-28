import numpy
import SALib
import json


from hera.utils.jsonutils import setJSONPath

class SALibUtils:
    
    @classmethod
    def buildSAProblem(cls, **kwargs):
        """
            Gets parameters and returns a dict with :
                - problem : a dict with the problem JSON
                - types   : a dict that correlates number->type and allows for the maping of the valus that were
                            samples to the real values (int, category, gaussian, log and ect.  
        Parameters
        ----------
        kwargs:
            A dict of variable name -> type.

            The types are:
                * float - a list with 2 number , lower and upper.
                * Distribution - dict(type='gaussian',params = dict of distribution parameters)
                                Not implemented yet.
                * int - dict(type="int",params=[range])
                * log - dict(type="log",params=[range])
                * List - dict(type="list", params=[list items])


        Returns
        -------
            dict

        """
        boundsList   = []
        typeDict = kwargs.copy()
        for key,value in kwargs.items():
            if isinstance(value,dict):
                boundsFunc = getattr(cls,f"_getBounds_{value['type']}")
                bounds = boundsFunc(value['parameters'])
                boundsList.append(bounds)
            else:
                # Then, it is a regular float number between lower...upper
                typeDict[key] = dict(type="float",parameters=value)
                boundsList.append([value[0],value[1]])

        problemJSON = dict(
            problem = dict(num_vars=len(kwargs.keys()),
                           names = [x for x in kwargs.keys()],
                           bounds = boundsList),
            type=typeDict
        )
        return problemJSON

    @classmethod
    def transformSample(cls, batchList, problemContainer):

        typeDict = problemContainer['type']
        newbatchList = []
        for batch in batchList:
            newValueList = []
            for paramName,paramValue in zip(problemContainer['problem']['names'],batch):
                transformerFunc = getattr(cls, f"_transform_{typeDict[paramName]['type']}")
                newValue = transformerFunc(paramValue,typeDict[paramName])
                newValueList.append(newValue)
            newbatchList.append(newValueList)

        return newbatchList

    @classmethod
    def typeInt(cls,lower,upper):
        return dict(type="int",parameters=[lower,upper+1])

    @classmethod
    def typeList(cls,items):
        return dict(type="list", parameters=items)

    @classmethod
    def typeLog(cls,lower,upper):
        return dict(type="log",parameters=[lower,upper])


    ##----------------------------------------------------------
    ##----------------------------------------------------------
    ##----------------------------------------------------------
    @classmethod
    def _getBounds_list(cls,parameters):
        return [0,len(parameters)]

    @classmethod
    def _getBounds_int(cls,parameters):
        return [parameters[0],parameters[1]+1]

    @classmethod
    def _getBounds_log(cls,parameters):
        return [parameters[0],parameters[1]+1]

    ##----------------------------------------------------------
    ##----------------------------------------------------------
    ##----------------------------------------------------------
    @classmethod
    def _transform_list(cls,value,meta):
        indx = int(numpy.floor(value))
        if indx==len(meta['parameters']):
            indx -=1
        return meta['parameters'][indx]

    @classmethod
    def _transform_int(cls,value,meta):
        return int(numpy.floor(value))

    @classmethod
    def _transform_log(cls,value,meta):
        return 10**float(value)

    @classmethod
    def _transform_float(cls,value,meta):
        return float(value)



