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

    def _dictTomongo(dictObj,local_perfix):

        for key,value in dictObj.items():
            new_prefix = key if local_perfix=="" else "%s__%s" % (local_perfix, key)
            if isinstance(value,dict):
                _dictTomongo(dictObj[key],local_perfix=new_prefix)
            else:
                ret[new_prefix] = value

    _dictTomongo(dictObj,local_perfix=prefix)
    return ret




