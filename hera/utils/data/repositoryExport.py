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
