def andClause(excludeFields=[], **kwargs):

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
