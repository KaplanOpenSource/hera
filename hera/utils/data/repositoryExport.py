"""
Pure helpers for exporting Hera Metadata documents into a repository JSON.

These functions operate only on plain dicts (the output of
``MetadataFrame.asDict(with_id=True)``) and plain repository-JSON dicts. They
perform no database access and no file IO, so they are fully unit-testable in
isolation. ``hera.utils.data.toolkit.dataToolkit.exportDocumentsToRepository``
is the thin facade that supplies the document dicts and writes the result.
"""

import copy
import hashlib
import json

# Section name (repository-JSON key) per document collection class.
_SECTION_BY_CLS = {
    "Measurements": "Measurements",
    "Simulations": "Simulations",
    "Cache": "Cache",
}


def _sourceObjectId(docDict):
    """Return the document ObjectId as a plain string, or None if absent.

    ``asDict(with_id=True)`` renders ``_id`` either as ``{"$oid": "..."}`` (the
    mongoengine/bson JSON form) or, when hand-built, as a plain string.
    """
    raw = docDict.get("_id")
    if raw is None:
        return None
    if isinstance(raw, dict):
        return raw.get("$oid")
    return str(raw)


def documentContentHash(docDict, idStrategy="contentHash"):
    """Return a stable identity string for a document dict.

    idStrategy="contentHash" (default): sha256 over the canonical JSON of
        {type, dataFormat, resource, desc}. The ``_id`` is intentionally
        excluded, so identical content in different projects collapses.
    idStrategy="objectId": the document ObjectId string is the identity.
    """
    if idStrategy == "objectId":
        oid = _sourceObjectId(docDict)
        if oid is None:
            raise ValueError("idStrategy='objectId' but the document has no _id")
        return oid

    if idStrategy != "contentHash":
        raise ValueError(f"Unknown idStrategy: {idStrategy!r}")

    payload = {
        "type": docDict.get("type"),
        "dataFormat": docDict.get("dataFormat"),
        "resource": docDict.get("resource"),
        "desc": docDict.get("desc", {}),
    }
    canonical = json.dumps(payload, sort_keys=True, default=str)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _sectionForCls(cls_value):
    """Map a document ``_cls`` ('Metadata.Measurements') to a repo section name."""
    if not cls_value or "." not in cls_value:
        raise ValueError(f"Document _cls is missing or malformed: {cls_value!r}")
    short = cls_value.split(".")[1]
    if short not in _SECTION_BY_CLS:
        raise ValueError(f"Unsupported document _cls: {cls_value!r}")
    return _SECTION_BY_CLS[short]


def documentToRepositoryItem(docDict, idStrategy="contentHash"):
    """Convert one document dict into a (section, itemName, entry) triple.

    The entry is a repository-JSON record:
        {"isRelativePath": "False", "contentHash": ..., "sourceId": ...,
         "item": {type, resource, dataFormat, desc}}
    Raises ValueError if ``_cls`` is missing or unrecognised.
    """
    section = _sectionForCls(docDict.get("_cls"))
    content_hash = documentContentHash(docDict, idStrategy=idStrategy)
    source_id = _sourceObjectId(docDict)
    item_name = source_id if source_id is not None else content_hash[:16]

    entry = {
        "isRelativePath": "False",
        "contentHash": content_hash,
        "sourceId": source_id,
        "item": {
            "type": docDict.get("type"),
            "resource": docDict.get("resource"),
            "dataFormat": docDict.get("dataFormat"),
            "desc": docDict.get("desc", {}),
        },
    }
    return section, item_name, entry


def _existingIdentities(toolkitDict):
    """Return {identity: (section, itemName)} for all entries under a toolkit.

    Identity prefers the stored ``contentHash``; falls back to ``sourceId`` for
    entries (e.g. hand-authored) that lack a hash.
    """
    index = {}
    for section, entries in toolkitDict.items():
        if not isinstance(entries, dict):
            continue
        for item_name, entry in entries.items():
            if not isinstance(entry, dict):
                continue
            identity = entry.get("contentHash") or entry.get("sourceId")
            if identity is not None:
                index[identity] = (section, item_name)
    return index


def mergeDocumentsIntoRepository(repoJSON, docDicts, toolkitName,
                                 idStrategy="contentHash", overwrite=False):
    """Merge document dicts under repoJSON[toolkitName].

    Returns (newRepoJSON, report). ``report`` has keys ``added``,
    ``skipped_existing`` and ``overwritten``, each a list of itemNames. A
    document already present (matching identity anywhere under ``toolkitName``)
    is skipped unless ``overwrite=True``. The input ``repoJSON`` is never
    mutated.
    """
    repo = copy.deepcopy(repoJSON)
    toolkitDict = repo.setdefault(toolkitName, {})
    index = _existingIdentities(toolkitDict)
    report = {"added": [], "skipped_existing": [], "overwritten": []}

    for docDict in docDicts:
        section, item_name, entry = documentToRepositoryItem(docDict, idStrategy=idStrategy)
        identity = entry.get("contentHash") or entry.get("sourceId")

        if identity in index:
            if not overwrite:
                report["skipped_existing"].append(item_name)
                continue
            # Remove the previously stored entry (possibly in another section).
            old_section, old_name = index[identity]
            toolkitDict.get(old_section, {}).pop(old_name, None)
            report["overwritten"].append(item_name)
        else:
            report["added"].append(item_name)

        sectionDict = toolkitDict.setdefault(section, {})
        item_name = _uniqueItemName(sectionDict, item_name, entry["contentHash"])
        sectionDict[item_name] = entry
        index[identity] = (section, item_name)

    return repo, report


def _uniqueItemName(sectionDict, item_name, content_hash):
    """Return an itemName not already used in ``sectionDict``.

    ObjectId-derived names are unique in practice, but if a collision occurs
    (e.g. distinct content sharing a name) the contentHash disambiguates so no
    entry is silently overwritten.
    """
    if item_name not in sectionDict:
        return item_name
    candidate = f"{item_name}_{content_hash[:8]}"
    suffix = 1
    while candidate in sectionDict:
        candidate = f"{item_name}_{content_hash[:8]}_{suffix}"
        suffix += 1
    return candidate


def deduplicateRepository(repoJSON):
    """Collapse entries sharing the same identity to a single entry.

    Identity is ``contentHash`` (or ``sourceId`` fallback). The first occurrence
    (iterating toolkits -> sections -> items) is kept; later duplicates are
    removed. Returns (newRepoJSON, report) where ``report['removed']`` lists
    (toolkit, section, itemName) tuples. Input is not mutated.
    """
    repo = copy.deepcopy(repoJSON)
    report = {"removed": []}

    for toolkitName, toolkitDict in repo.items():
        if not isinstance(toolkitDict, dict):
            continue
        seen = set()
        for section, entries in toolkitDict.items():
            if not isinstance(entries, dict):
                continue
            for item_name in list(entries.keys()):
                entry = entries[item_name]
                if not isinstance(entry, dict):
                    continue
                identity = entry.get("contentHash") or entry.get("sourceId")
                if identity is None:
                    continue
                if identity in seen:
                    del entries[item_name]
                    report["removed"].append((toolkitName, section, item_name))
                else:
                    seen.add(identity)

    return repo, report
