"""``AbstractCollection.addDocumentFromJSON``: importing a document from its
own JSON serialisation.

This is the read side of the export/import round trip.  Its counterparts
are ``MetadataFrame.asDict(with_id=True)`` and ``Document.to_json()``, both
of which emit an ``_id`` field, and ``Project.exportProject``, which is
built on ``asDict``.  So the shape this method is fed in practice always
carries an ``_id``.

An earlier batch recorded the symptom as "``.save()`` reports success and a
new id, but the document is not found immediately after by that id" and
flagged it as possibly a mongomock quirk (see the docstring of
``test_datalayer_project_export.py``).  It was not a mongomock quirk: it was
``mongoengine.Document.from_json``'s ``created=False`` default, whose own
docstring spells out the rule --

    If False and an ID is provided, assume that the object has already been
    persisted (this has an impact on the subsequent call to .save()).

-- so a body carrying an ``_id`` was marked already-persisted with an empty
change set and ``.save()`` wrote nothing.  ``addDocumentFromJSON`` now passes
``created=True``.  The tests below cover both halves: JSON without an ``_id``
and JSON with one.

Deliberately not covered here: the rest of ``AbstractCollection``
(``getDocuments``, ``addDocument``, ``deleteDocuments``,
``deleteDocumentByID``, ``getProjectList``, ``getDocumentsAsDict``), which
is already exercised through ``Project`` in
``test_datalayer_project*.py``.
"""
import json

import pytest

from hera.tests.unit.conftest import UNIT_PROJECT_NAME


@pytest.fixture()
def collection():
    from hera.datalayer import Measurements_Collection

    return Measurements_Collection()


def _document_body(**overrides):
    """The minimal JSON body a Measurements document validates against."""
    body = {
        "_cls": "Metadata.Measurements",
        "projectName": UNIT_PROJECT_NAME,
        "type": "ImportedType",
        "resource": "/data/imported.parquet",
        "dataFormat": "parquet",
        "desc": {"station": "YAVNEEL"},
    }
    body.update(overrides)
    return body


def _stored(unit_project, type="ImportedType"):
    return list(unit_project.getMeasurementsDocuments(type=type))


# ---------------------------------------------------------------------------
# The working half: JSON with no _id
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestImportingJSONWithoutAnID:
    def test_the_document_is_inserted(self, collection, unit_project):
        collection.addDocumentFromJSON(json.dumps(_document_body()))

        assert len(_stored(unit_project)) == 1

    def test_every_field_survives_the_round_trip(self, collection, unit_project):
        collection.addDocumentFromJSON(json.dumps(_document_body()))

        doc = _stored(unit_project)[0]
        assert doc["projectName"] == UNIT_PROJECT_NAME
        assert doc["type"] == "ImportedType"
        assert doc["resource"] == "/data/imported.parquet"
        assert doc["dataFormat"] == "parquet"
        assert doc.desc == {"station": "YAVNEEL"}

    def test_the_document_lands_in_the_collection_it_was_called_on(self, collection,
                                                                   unit_project):
        """CLAUDE.md: measurements and simulations must never mix.  The
        collection object, not the JSON, decides where the row goes."""
        collection.addDocumentFromJSON(json.dumps(_document_body()))

        assert len(_stored(unit_project)) == 1
        assert list(unit_project.getSimulationsDocuments(type="ImportedType")) == []

    def test_a_fresh_id_is_assigned_by_the_database(self, collection, unit_project):
        collection.addDocumentFromJSON(json.dumps(_document_body()))

        assert _stored(unit_project)[0].id is not None

    def test_importing_the_same_body_twice_inserts_two_documents(self, collection,
                                                                 unit_project):
        """No _id means "brand new" each time, so this is an insert, not an
        upsert."""
        body = json.dumps(_document_body())
        collection.addDocumentFromJSON(body)
        collection.addDocumentFromJSON(body)

        assert len(_stored(unit_project)) == 2

    def test_the_imported_document_can_be_queried_by_its_description(self, collection,
                                                                     unit_project):
        collection.addDocumentFromJSON(json.dumps(_document_body()))

        assert len(list(unit_project.getMeasurementsDocuments(
            type="ImportedType", station="YAVNEEL"))) == 1

    def test_an_omitted_description_is_allowed(self, collection, unit_project):
        body = _document_body()
        body.pop("desc")
        collection.addDocumentFromJSON(json.dumps(body))

        assert _stored(unit_project)[0].desc == {}

    def test_nothing_is_returned(self, collection):
        """Unlike addDocument, which hands back the mongoengine document, this
        method returns None -- so a caller cannot reach what it inserted."""
        assert collection.addDocumentFromJSON(json.dumps(_document_body())) is None

    @pytest.mark.parametrize("missing", ["projectName", "type", "resource", "dataFormat"])
    def test_a_body_missing_a_required_field_is_refused(self, collection, unit_project,
                                                        missing):
        from mongoengine import ValidationError

        body = _document_body()
        body.pop(missing)
        with pytest.raises(ValidationError, match="Field is required"):
            collection.addDocumentFromJSON(json.dumps(body))

        assert _stored(unit_project) == []

    def test_an_empty_body_is_refused(self, collection):
        from mongoengine import ValidationError

        with pytest.raises(ValidationError, match="Field is required"):
            collection.addDocumentFromJSON("{}")

    def test_malformed_json_is_refused(self, collection):
        with pytest.raises(Exception):
            collection.addDocumentFromJSON("{not json")


# ---------------------------------------------------------------------------
# JSON carrying an _id
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestImportingJSONWithAnID:
    """See the module docstring: the body is imported under its own _id."""

    @staticmethod
    def _exported(unit_project, resource="/data/original.parquet"):
        """A document as the datalayer itself serialises it, then removed --
        i.e. exactly the export/re-import case."""
        doc = unit_project.addMeasurementsDocument(
            resource=resource, dataFormat="parquet", type="ImportedType",
            desc={"station": "YAVNEEL"})
        body = doc.to_json()
        doc.delete()
        return body

    def test_an_exported_document_should_be_reimported(self, collection, unit_project):
        collection.addDocumentFromJSON(self._exported(unit_project))

        assert len(_stored(unit_project)) == 1

    def test_the_same_body_without_its_id_is_inserted_too(self, collection,
                                                          unit_project):
        """The ``_id`` field makes no difference to whether the body lands."""
        body = json.loads(self._exported(unit_project))
        assert "_id" in body
        body.pop("_id")

        collection.addDocumentFromJSON(json.dumps(body))

        assert len(_stored(unit_project)) == 1

    def test_an_id_that_does_exist_is_replaced(self, collection, unit_project):
        """Re-importing a body whose id is still live overwrites that
        document rather than silently doing nothing."""
        doc = unit_project.addMeasurementsDocument(
            resource="/data/original.parquet", dataFormat="parquet",
            type="ImportedType", desc={})
        body = json.loads(doc.to_json())
        body["resource"] = "/data/replacement.parquet"

        collection.addDocumentFromJSON(json.dumps(body))

        remaining = _stored(unit_project)
        assert len(remaining) == 1
        assert remaining[0]["resource"] == "/data/replacement.parquet"

    def test_the_rebuilt_document_is_marked_as_new(self, collection):
        """The mechanism the fix rests on: without created=True mongoengine
        marks a body carrying an _id as already persisted, with no pending
        changes for .save() to write."""
        body = json.dumps(_document_body(_id={"$oid": "0" * 24}))

        assert collection._metadataCol.from_json(body)._created is False
        assert collection._metadataCol.from_json(body, created=True)._created is True
