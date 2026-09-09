from __future__ import annotations

from typing import Any, Dict


class SentWorkflowDoc:
    """Doc-like stand-in for a saved workflow document, built from the dict the
    client sent instead of read from the DB.

    A real saved document is read two ways by ``prepareWorkflowRunFromDoc``:
    attribute access for ``doc.desc`` and item access for ``doc['resource']``. This
    wraps the sent dict to support both, so the same prep runs with no DB lookup.
    """

    def __init__(self, doc: Dict[str, Any]) -> None:
        self._doc = doc

    @property
    def desc(self) -> Dict[str, Any]:
        return self._doc["desc"]

    def __getitem__(self, key: str) -> Any:
        return self._doc[key]
