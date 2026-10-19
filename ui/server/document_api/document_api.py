"""One document, as the named document routes see it."""
from __future__ import annotations

from typing import Any, Dict, List

from api_models import Collection

# Fields the editor must never send back: they say which document this is, and
# which project and collection it belongs to.
FORBIDDEN_FIELDS = ("_id", "_cls", "projectName")


class Refused(Exception):
    """The request asks for something the document layer will not do."""


class NoSuchDocument(Exception):
    """No document in the database has this id."""


class DocumentApi:
    """What the UI can do to one document.

    Holds nothing but the document's id: it is built for one request and
    dropped, so every call reads the database again.
    """

    def __init__(self, docId: str):
        self.docId = docId

    @staticmethod
    def _collection(collection: Collection):
        from hera.datalayer import (Cache_Collection, Measurements_Collection,
                                    Simulations_Collection)
        classes = {
            Collection.Measurements: Measurements_Collection,
            Collection.Simulations: Simulations_Collection,
            Collection.Cache: Cache_Collection,
        }
        try:
            return classes[Collection(collection)]()
        except ValueError as error:
            raise Refused(f"no such collection: {collection}") from error

    @staticmethod
    def create(projectName: str,
               collection: Collection,
               type: str = "",
               resource: Any = "",
               dataFormat: str = "string",
               desc: Dict[str, Any] = None) -> Dict[str, Any]:
        """Add a document to the collection, and answer it with its id."""
        if desc is not None and not isinstance(desc, dict):
            raise Refused("desc must be a dict")
        added = DocumentApi._collection(collection).addDocument(
            projectName=projectName,
            resource=resource,
            dataFormat=dataFormat,
            type=type,
            desc=desc or {})
        return added.asDict(with_id=True)

    @staticmethod
    def one_project(docIds: List[str]) -> str:
        """The one project every given document belongs to.

        Read before anything is deleted, so the answer still has a project to carry.
        """
        names = {DocumentApi(docId).project_name() for docId in docIds}
        if len(names) > 1:
            raise Refused(
                f"the documents are from more than one project: {', '.join(sorted(names))}")
        return names.pop()

    def _document(self):
        from hera.datalayer import All
        try:
            return All.getDocumentByID(self.docId)
        except Exception as error:
            raise NoSuchDocument(f"no such document: {self.docId}") from error

    def project_name(self) -> str:
        """The project this document belongs to."""
        return self._document().projectName

    def update(self, fields: Dict[str, Any]) -> Dict[str, Any]:
        """Write the given top-level fields, and answer the saved document."""
        forbidden = [f for f in fields if f in FORBIDDEN_FIELDS]
        if forbidden:
            raise Refused(f"cannot change {', '.join(sorted(forbidden))}")
        document = self._document()
        for field, value in fields.items():
            setattr(document, field, value)
        document.save()
        return self._document().asDict(with_id=True)

    def delete(self) -> None:
        """Remove the document from the database."""
        from hera.datalayer import All
        # Reading it first is what turns an unknown id into a 404.
        self._document()
        All.deleteDocumentByID(self.docId)
