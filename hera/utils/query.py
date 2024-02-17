import pandas
import json

def andClause(excludeFields=[], **kwargs):
    """
        Builds a pandas query str
    Parameters
    ----------
    excludeFields
    kwargs

    Returns
    -------

    """

    L = []
    for key, value in kwargs.items():
        if key in excludeFields:
            continue

        if isinstance(value, list):
            conditionStr = "%s in %s"
        elif isinstance(value, str):
            conditionStr = "%s == '%s'"
        elif isinstance(value, dict):
            conditionStr = "%s " + value['operator'] + " %s"
            value = value['value']
        else:
            conditionStr = "%s == %s"

        L.append(conditionStr % (key, value))

    return " and ".join(L)


def dictToMongoQuery(dictObj,prefix=""):
    """
        Convets a dict object to a mongodb query.

        That is, if the JSON is:

    .. code-block:: JSON

        "fieldname1" : {
                "subfield" : 1,
                "subfield1" : "d"
        },
        "fieldname2" : {
                "subfield3" : "hello",
                "subfield4" : "goodbye",
        }


    translates to the dict:

    .. code-block:: python

        {
            "fieldname1__subfield"  : 1,
            "fieldname1__subfield1"  : "d",
            "fieldname2__subfield3"  : "hello",
            "fieldname2__subfield4"  : "goodbye",
        }

    if the

    This will allow to use the returned dict in collection.getDocumets.

    :param  prefix: prefix for the mongoDB query fields.
                    used to select only a part of the document description.

    :param dictObj:
            A dictionary with fields and values.
    :return:
        dict
    """
    ret = {}

    def determineType(value, prefix):
        if isinstance(value, dict):
            _dictTomongo(value, local_perfix=prefix)
        elif isinstance(value, list):
            for indx,listValue in enumerate(value):
                new_prefix = f"{prefix}__{indx}"
                determineType(listValue,new_prefix)
        else:
            ret[prefix] = value

    def _dictTomongo(dictObj,local_perfix):
        for key,value in dictObj.items():
            if key in ['type','in','ne','lt','lte','gt','gte','not','all','size','exists','nin','not']:
                key = f"{key}__"

            new_prefix = key if local_perfix=="" else "%s__%s" % (local_perfix, key)
            determineType(value,new_prefix)

    _dictTomongo(dictObj,local_perfix=prefix)
    return ret




