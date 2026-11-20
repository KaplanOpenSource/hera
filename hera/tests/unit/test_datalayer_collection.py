"""``AbstractCollection.addDocumentFromJSON``: importing a document from its
own JSON serialisation.

This is the read side of the export/import round trip.  Its counterparts
are ``MetadataFrame.asDict(with_id=True)`` and ``Document.to_json()``, both
of which emit an ``_id`` field, and ``Project.exportProject``, which is
built on ``asDict``.  So the shape this method is fed in practice always
carries an ``_id``.

An earlier batch left this method uncovered, recording the symptom as
"``.save()`` reports success and a new id, but the document is not found
immediately after by that id" and flagging it as possibly a mongomock
quirk (see the docstring of ``test_datalayer_project_export.py``).  It is
not a mongomock quirk: it is ``mongoengine.Document.from_json``'s
documented ``created=False`` default, and the tests below pin it as B261.
``from_json``'s own docstring spells out the rule --

    If False and an ID is provided, assume that the object has already been
    persisted (this has an impact on the subsequent call to .save()).

-- so a document rebuilt from JSON that contains an ``_id`` is marked as
already-persisted with an empty change set, and ``.save()`` degenerates to
an update that matches nothing (when the id is absent from the collection)
or that changes nothing (when it is present).  Either way no error is
raised.  The tests separate the two halves: JSON *without* an ``_id`` is
inserted correctly and is asserted normally; JSON *with* an ``_id`` is
xfailed strict plus characterised.

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
# The broken half: JSON carrying an _id
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestImportingJSONWithAnIDIsANoOp:
    """B261: see the module docstring."""

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

    @pytest.mark.xfail(
        strict=True,
        reason="B261: addDocumentFromJSON calls "
               "self._metadataCol.from_json(json_data).save(), and "
               "mongoengine's from_json defaults to created=False, which its "
               "own docstring defines as 'if an ID is provided, assume that "
               "the object has already been persisted'.  The rebuilt document "
               "therefore has _created=False and an empty _changed_fields, so "
               ".save() issues an update instead of an insert: nothing is "
               "written and no error is raised.  Since asDict(with_id=True), "
               "Document.to_json() and Project.exportProject all emit _id, "
               "re-importing exported documents silently drops every one of "
               "them.  The fix is from_json(json_data, created=True).  This "
               "supersedes the 'inconsistent mongomock behavior' note in "
               "test_datalayer_project_export.py -- the cause is mongoengine, "
               "not mongomock.  See the consolidated findings issue.",
    )
    def test_an_exported_document_should_be_reimported(self, collection, unit_project):
        collection.addDocumentFromJSON(self._exported(unit_project))

        assert len(_stored(unit_project)) == 1

    def test_nothing_is_inserted_and_nothing_is_raised(self, collection, unit_project):
        """Characterisation of B261."""
        collection.addDocumentFromJSON(self._exported(unit_project))

        assert _stored(unit_project) == []

    def test_the_same_body_without_its_id_does_get_inserted(self, collection,
                                                            unit_project):
        """Characterisation of B261: this isolates the ``_id`` field as the
        whole cause -- the body is otherwise byte-identical."""
        body = json.loads(self._exported(unit_project))
        assert "_id" in body
        body.pop("_id")

        collection.addDocumentFromJSON(json.dumps(body))

        assert len(_stored(unit_project)) == 1

    def test_an_id_that_does_exist_is_not_updated_either(self, collection, unit_project):
        """Characterisation of B261: the change set is empty as well as the
        insert being suppressed, so the method cannot even be used to
        overwrite a live document."""
        doc = unit_project.addMeasurementsDocument(
            resource="/data/original.parquet", dataFormat="parquet",
            type="ImportedType", desc={})
        body = json.loads(doc.to_json())
        body["resource"] = "/data/replacement.parquet"

        collection.addDocumentFromJSON(json.dumps(body))

        remaining = _stored(unit_project)
        assert len(remaining) == 1
        assert remaining[0]["resource"] == "/data/original.parquet"

    def test_mongoengine_confirms_the_created_default(self, collection):
        """Characterisation of B261's mechanism, so the diagnosis does not
        rest on the symptom alone: the document comes back already marked as
        persisted, with no pending changes for .save() to write."""
        body = json.dumps(_document_body(_id={"$oid": "0" * 24}))
        rebuilt = collection._metadataCol.from_json(body)

        assert rebuilt._created is False
        assert rebuilt._get_changed_fields() == []
