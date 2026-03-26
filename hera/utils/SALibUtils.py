import numpy
import SALib
import json


from hera.utils.jsonutils import setJSONPath

class SALibUtils:
    """Utilities for building and transforming SALib sensitivity analysis problems."""

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
        """Transform raw SALib samples to their actual typed values."""
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
        """Build a type descriptor for an integer parameter."""
        return dict(type="int",parameters=[lower,upper+1])

    @classmethod
    def typeList(cls,items):
        """Build a type descriptor for a categorical list parameter."""
        return dict(type="list", parameters=items)

    @classmethod
    def typeLog(cls,lower,upper):
        """Build a type descriptor for a log-scale parameter."""
        return dict(type="log",parameters=[lower,upper])


    ##----------------------------------------------------------
    ##----------------------------------------------------------
    ##----------------------------------------------------------
    @classmethod
    def _getBounds_list(cls,parameters):
        """Return sampling bounds for a list-type parameter."""
        return [0,len(parameters)]

    @classmethod
    def _getBounds_int(cls,parameters):
        """Return sampling bounds for an integer-type parameter."""
        return [parameters[0],parameters[1]+1]

    @classmethod
    def _getBounds_log(cls,parameters):
        """Return sampling bounds for a log-scale parameter."""
        return [parameters[0],parameters[1]+1]

    ##----------------------------------------------------------
    ##----------------------------------------------------------
    ##----------------------------------------------------------
    @classmethod
    def _transform_list(cls,value,meta):
        """Map a sampled float to the corresponding list item."""
        indx = int(numpy.floor(value))
        if indx==len(meta['parameters']):
            indx -=1
        return meta['parameters'][indx]

    @classmethod
    def _transform_int(cls,value,meta):
        """Transform a sampled float to an integer by flooring."""
        return int(numpy.floor(value))

    @classmethod
    def _transform_log(cls,value,meta):
        """Transform a sampled value to log-scale (10**value)."""
        return 10**float(value)

    @classmethod
    def _transform_float(cls,value,meta):
        """Transform a sampled value to a Python float."""
        return float(value)



