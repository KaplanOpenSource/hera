export const IS_REPO_JSON_PYTHON = `
def loadRepoJson(path):
    SECTIONS = {'config', 'datasource', 'measurements', 'simulations', 'cache', 'function'}
    try:
        with open(path, 'r') as f:
            doc = json.load(f)
    except Exception:
        return None
    if not isinstance(doc, dict) or not doc:
        return None
    for toolkitValue in doc.values():
        if not isinstance(toolkitValue, dict):
            return None
        for sectionKey in toolkitValue.keys():
            if sectionKey.lower() not in SECTIONS:
                return None
    return doc
`;
