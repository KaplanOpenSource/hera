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
