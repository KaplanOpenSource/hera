def andClause(excludeFields=None, **kwargs):
    """
        Builds a pandas query str
    Parameters
    ----------
    excludeFields
    kwargs

    Returns
    -------

    """
    if excludeFields is None:
        excludeFields = []

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
        Converts a dict object to a mongodb query.

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

    This will allow to use the returned dict in collection.getDocuments.

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
    # The algorithm works by recursively walking the nested dict structure
    # and building flattened keys using MongoEngine's double-underscore (__)
    # convention for nested field access.
    #
    # Example: {"config": {"model": "LSM"}} → {"config__model": "LSM"}
    #
    # The prefixExclude parameter (default "desc") strips a known wrapper
    # key from the path. This is needed because Hera documents store user
    # metadata under "desc", but MongoEngine queries on desc sub-fields
    # use "desc__fieldname" — so when building queries from a dict that
    # already has "desc" as a top-level key, we skip adding "desc" to the
    # prefix to avoid generating "desc__desc__fieldname".

    ret = {}

    def determineType(value, prefix,prefixExclude):
        """Dispatch value by type and recurse or store in ret."""
        if isinstance(value, dict):
            # Nested dict: recurse deeper, extending the __ prefix chain.
            _dictToMongo(value, local_prefix=prefix,prefixExclude=prefixExclude)
        elif isinstance(value, list):
            # Lists: index each element (e.g. "field__0", "field__1").
            for indx,listValue in enumerate(value):
                new_prefix = f"{prefix}__{indx}"
                determineType(listValue,new_prefix,prefixExclude)
        else:
            # Leaf value: store with the accumulated prefix as key.
            ret[prefix] = value

    def _dictToMongo(dictObj,local_prefix,prefixExclude):
        """Recursively flatten a dict into mongo-style __ separated keys."""
        for key,value in dictObj.items():
            if key==prefixExclude:
                # Skip the excluded prefix (e.g. "desc") — don't add it to
                # the key chain. This lets us pass {"desc": {"field": 1}}
                # and get {"field": 1} instead of {"desc__field": 1}.
                new_prefix = local_prefix
            else:
                new_prefix = key if local_prefix=="" else "%s__%s" % (local_prefix, key)

            determineType(value,new_prefix,prefixExclude)

    _dictToMongo(dictObj,local_prefix=prefix,prefixExclude=prefixExclude)

    # MongoEngine workaround: keys ending with "__type" clash with
    # MongoEngine's internal type discrimination. Appending an extra "__"
    # escapes the collision (e.g. "desc__type" → "desc__type__").
    keyList = [key for key in ret.keys()]
    for key in keyList:
        if key.endswith("__type"):
            ret[f"{key}__"] = ret[key]
            del ret[key]

    return ret




