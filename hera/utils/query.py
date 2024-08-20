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


def dictToMongoQuery(dictObj,prefix="",prefixExclude="desc"):
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

    prefixExclude : str
            Prefix to exclude from the addition to the dict.
            This is used to query a document from the DB, and we to pass to the query
            the parameters without the 'desc'



    :param dictObj:
            A dictionary with fields and values.
    :return:
        dict
    """
    ret = {}

    def determineType(value, prefix,prefixExclude):
        if isinstance(value, dict):
            _dictTomongo(value, local_perfix=prefix,prefixExclude=prefixExclude)
        elif isinstance(value, list):
            for indx,listValue in enumerate(value):
                new_prefix = f"{prefix}__{indx}"
                determineType(listValue,new_prefix,prefixExclude)
        else:
            ret[prefix] = value

    def _dictTomongo(dictObj,local_perfix,prefixExclude):
        for key,value in dictObj.items():
            if key==prefixExclude:
                new_prefix = local_perfix
            else:
                new_prefix = key if local_perfix=="" else "%s__%s" % (local_perfix, key)

            determineType(value,new_prefix,prefixExclude)

    _dictTomongo(dictObj,local_perfix=prefix,prefixExclude=prefixExclude)

    # Now convert all the keys that end with __type to __type__.
    # This fixes the pymongo queyr where __type is used to determine the type of the format.
    keyList = [key for key in ret.keys()]
    for key in keyList:
        if key.endswith("__type"):
            ret[f"{key}__"] = ret[key]
            del ret[key]

    return ret




