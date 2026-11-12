"""The document entrypoints: create, update, delete."""
import json
import sys
import types

import pytest
from fastapi import HTTPException

import dependencies
from api_models import Collection
from document_api.document_api import DocumentApi, NoSuchDocument, Refused
from document_api.notebook_file import NotebookFile
from routers import documents as routes


class FakeDocument:
    """A stored document: fields the route writes, and a save that keeps them."""

    def __init__(self, store, docId, projectName, **fields):
        self._store = store
        self._id = docId
        self.projectName = projectName
        self.resource = fields.get("resource", "")
        self.dataFormat = fields.get("dataFormat", "string")
        self.type = fields.get("type", "")
        self.desc = fields.get("desc", {})
        self.saves = 0

    def save(self):
        self.saves += 1

    def asDict(self, with_id=False):
        asDict = {
            "_cls": "Metadata.Cache",
            "projectName": self.projectName,
            "resource": self.resource,
            "dataFormat": self.dataFormat,
            "type": self.type,
            "desc": self.desc,
        }
        if with_id:
            asDict["_id"] = {"$oid": self._id}
        return asDict


class FakeDatabase:
    """The little bit of hera the document entrypoints touch, in memory."""

    def __init__(self):
        self.documents = {}
        self.addedTo = []
        self.nextId = 1

    def add(self, docId, projectName, **fields):
        self.documents[docId] = FakeDocument(self, docId, projectName, **fields)
        return self.documents[docId]

    # --- hera.datalayer.All ---

    def getDocumentByID(self, docId):
        return self.documents[docId]

    def deleteDocumentByID(self, docId):
        del self.documents[docId]

    def getDocumentsAsDict(self, projectName, with_id=False):
        documents = [d.asDict(with_id=with_id) for d in self.documents.values()
                     if d.projectName == projectName]
        return {"documents": documents}

    # --- the three collection classes ---

    def collection(self, name):
        def make():
            return types.SimpleNamespace(
                addDocument=lambda projectName, resource, dataFormat, type, desc:
                    self._addDocument(name, projectName, resource, dataFormat, type, desc))
        return make

    def _addDocument(self, collectionName, projectName, resource, dataFormat, type, desc):
        docId = f"new{self.nextId}"
        self.nextId += 1
        self.addedTo.append(collectionName)
        return self.add(docId, projectName, resource=resource,
                        dataFormat=dataFormat, type=type, desc=desc)


@pytest.fixture
def fake_hera(monkeypatch):
    """Install the fake hera and return the database behind it."""
    db = FakeDatabase()

    datalayer = types.ModuleType("hera.datalayer")
    datalayer.All = types.SimpleNamespace(
        getDocumentByID=db.getDocumentByID,
        deleteDocumentByID=db.deleteDocumentByID,
        getDocumentsAsDict=db.getDocumentsAsDict)
    datalayer.Measurements_Collection = db.collection("Measurements")
    datalayer.Simulations_Collection = db.collection("Simulations")
    datalayer.Cache_Collection = db.collection("Cache")

    hera = types.ModuleType("hera")
    hera.datalayer = datalayer
    monkeypatch.setitem(sys.modules, "hera", hera)
    monkeypatch.setitem(sys.modules, "hera.datalayer", datalayer)
    return db


@pytest.fixture
def warmed(monkeypatch):
    monkeypatch.setattr(dependencies.warmup, "_ready", True)


def _create(**fields):
    payload = {
        "projectName": "Alpha",
        "collection": Collection.Measurements,
        "type": "",
        "resource": "",
        "dataFormat": "string",
        "desc": {},
        "notebook": None,
    }
    payload.update(fields)
    return types.SimpleNamespace(**payload)


class TestCreateDocument:
    def test_adds_to_the_collection_it_is_given(self, fake_hera):
        document = DocumentApi.create("Alpha", Collection.Simulations,
                                      type="hermesWorkflow", resource="/w.json",
                                      dataFormat="JSON_dict", desc={"datasourceName": "W"})

        assert fake_hera.addedTo == ["Simulations"]
        assert document["type"] == "hermesWorkflow"
        assert document["_id"]["$oid"] == "new1"

    def test_a_resource_can_be_data_and_not_a_path(self, fake_hera):
        document = DocumentApi.create("Alpha", Collection.Measurements,
                                      resource={"effects": {}})

        assert document["resource"] == {"effects": {}}

    def test_refuses_an_unknown_collection(self, fake_hera):
        with pytest.raises(Refused):
            DocumentApi.create("Alpha", "Nowhere")

    def test_refuses_a_desc_that_is_not_a_dict(self, fake_hera):
        with pytest.raises(Refused):
            DocumentApi.create("Alpha", Collection.Measurements, desc=["not", "a", "dict"])

    def test_route_answers_the_document_its_id_and_the_project(self, fake_hera, warmed):
        answer = routes.create(_create(desc={"datasourceName": "MyDoc"}))

        assert answer["addedDocId"] == "new1"
        assert answer["document"]["desc"] == {"datasourceName": "MyDoc"}
        assert answer["project"]["name"] == "Alpha"
        assert len(answer["project"]["documents"]) == 1

    def test_route_rejects_a_blank_project_name(self, fake_hera, warmed):
        with pytest.raises(HTTPException) as raised:
            routes.create(_create(projectName="  "))

        assert raised.value.status_code == 400

    def test_route_rejects_an_unknown_collection(self, fake_hera, warmed):
        with pytest.raises(HTTPException) as raised:
            routes.create(_create(collection="Nowhere"))

        assert raised.value.status_code == 400


class TestCreateNotebookDocument:
    def test_route_writes_the_notebook_file(self, fake_hera, warmed, tmp_path):
        path = tmp_path / "notebooks" / "N.ipynb"

        routes.create(_create(
            collection=Collection.Cache, type="notebook", resource=str(path),
            dataFormat="JSON_dict", desc={"datasourceName": "N"},
            notebook=types.SimpleNamespace(toolkitNames=["LSM"])))

        written = json.loads(path.read_text())
        source = "".join(written["cells"][0]["source"])
        assert 'PROJECT_NAME = "Alpha"' in source
        assert 'getToolkit("LSM"' in source

    def test_route_leaves_an_existing_notebook_alone(self, fake_hera, warmed, tmp_path):
        path = tmp_path / "N.ipynb"
        path.write_text('{"mine": true}')

        routes.create(_create(
            collection=Collection.Cache, type="notebook", resource=str(path),
            notebook=types.SimpleNamespace(toolkitNames=[])))

        assert json.loads(path.read_text()) == {"mine": True}

    def test_route_rejects_a_notebook_without_a_path(self, fake_hera, warmed):
        with pytest.raises(HTTPException) as raised:
            routes.create(_create(
                resource="", notebook=types.SimpleNamespace(toolkitNames=[])))

        assert raised.value.status_code == 400

    def test_route_touches_no_file_without_a_notebook_block(self, fake_hera, warmed, tmp_path):
        path = tmp_path / "N.ipynb"

        routes.create(_create(resource=str(path)))

        assert not path.exists()


class TestNotebookFile:
    def test_writes_an_empty_notebook_with_the_project_in_it(self, tmp_path):
        path = tmp_path / "deep" / "N.ipynb"

        assert NotebookFile(str(path)).write_if_missing("Alpha", [])

        written = json.loads(path.read_text())
        assert written["nbformat"] == 4
        assert "".join(written["cells"][0]["source"]).count('PROJECT_NAME = "Alpha"') == 1

    def test_one_line_per_toolkit(self, tmp_path):
        path = tmp_path / "N.ipynb"

        NotebookFile(str(path)).write_if_missing("Alpha", ["LSM", "GIS Raster"])

        source = "".join(json.loads(path.read_text())["cells"][0]["source"])
        assert 'LSM_tk = hera.toolkit.ToolkitHome().getToolkit("LSM"' in source
        # The toolkit name is not always a legal python name.
        assert 'GIS_Raster_tk = hera.toolkit.ToolkitHome().getToolkit("GIS Raster"' in source

    def test_says_no_when_a_file_is_already_there(self, tmp_path):
        path = tmp_path / "N.ipynb"
        path.write_text("mine")

        assert not NotebookFile(str(path)).write_if_missing("Alpha", [])
        assert path.read_text() == "mine"


class TestUpdateDocument:
    def test_writes_the_fields_and_saves(self, fake_hera):
        stored = fake_hera.add("doc1", "Alpha", resource="old")

        saved = DocumentApi("doc1").update({"resource": "new", "type": "T"})

        assert stored.resource == "new"
        assert stored.saves == 1
        assert saved["resource"] == "new"

    def test_refuses_a_forbidden_field(self, fake_hera):
        fake_hera.add("doc1", "Alpha")

        with pytest.raises(Refused):
            DocumentApi("doc1").update({"resource": "new", "projectName": "Beta"})

    def test_a_refused_update_writes_nothing(self, fake_hera):
        stored = fake_hera.add("doc1", "Alpha", resource="old")

        with pytest.raises(Refused):
            DocumentApi("doc1").update({"resource": "new", "_cls": "Other"})

        assert stored.resource == "old"
        assert stored.saves == 0

    def test_an_unknown_id_is_not_found(self, fake_hera):
        with pytest.raises(NoSuchDocument):
            DocumentApi("nope").update({"resource": "new"})

    def test_route_answers_the_document_and_the_project(self, fake_hera, warmed):
        fake_hera.add("doc1", "Alpha", resource="old")
        fake_hera.add("doc2", "Alpha")

        answer = routes.update(types.SimpleNamespace(docId="doc1", fields={"resource": "new"}))

        assert answer["document"]["resource"] == "new"
        assert answer["project"]["name"] == "Alpha"
        assert len(answer["project"]["documents"]) == 2

    def test_route_rejects_a_bad_doc_id(self, fake_hera, warmed):
        with pytest.raises(HTTPException) as raised:
            routes.update(types.SimpleNamespace(docId="nope", fields={"resource": "new"}))

        assert raised.value.status_code == 404

    def test_route_rejects_an_empty_doc_id(self, fake_hera, warmed):
        with pytest.raises(HTTPException) as raised:
            routes.update(types.SimpleNamespace(docId="", fields={}))

        assert raised.value.status_code == 400

    def test_route_rejects_a_forbidden_field(self, fake_hera, warmed):
        fake_hera.add("doc1", "Alpha")

        with pytest.raises(HTTPException) as raised:
            routes.update(types.SimpleNamespace(
                docId="doc1", fields={"_id": {"$oid": "other"}}))

        assert raised.value.status_code == 400
        assert "_id" in raised.value.detail


class TestDeleteDocuments:
    def test_removes_the_document(self, fake_hera):
        fake_hera.add("doc1", "Alpha")

        DocumentApi("doc1").delete()

        assert "doc1" not in fake_hera.documents

    def test_an_unknown_id_is_not_found(self, fake_hera):
        with pytest.raises(NoSuchDocument):
            DocumentApi("nope").delete()

    def test_one_project_for_documents_of_the_same_project(self, fake_hera):
        fake_hera.add("doc1", "Alpha")
        fake_hera.add("doc2", "Alpha")

        assert DocumentApi.one_project(["doc1", "doc2"]) == "Alpha"

    def test_one_project_refuses_ids_from_two_projects(self, fake_hera):
        fake_hera.add("doc1", "Alpha")
        fake_hera.add("doc2", "Beta")

        with pytest.raises(Refused):
            DocumentApi.one_project(["doc1", "doc2"])

    def test_route_answers_the_project_left_behind(self, fake_hera, warmed):
        fake_hera.add("doc1", "Alpha")
        fake_hera.add("doc2", "Alpha")
        fake_hera.add("doc3", "Alpha")

        answer = routes.delete(types.SimpleNamespace(docIds=["doc1", "doc2"]))

        assert answer["project"]["name"] == "Alpha"
        assert len(answer["project"]["documents"]) == 1

    def test_route_rejects_an_empty_id_list(self, fake_hera, warmed):
        with pytest.raises(HTTPException) as raised:
            routes.delete(types.SimpleNamespace(docIds=[]))

        assert raised.value.status_code == 400

    def test_route_rejects_ids_from_two_projects(self, fake_hera, warmed):
        fake_hera.add("doc1", "Alpha")
        fake_hera.add("doc2", "Beta")

        with pytest.raises(HTTPException) as raised:
            routes.delete(types.SimpleNamespace(docIds=["doc1", "doc2"]))

        assert raised.value.status_code == 400
        assert "doc1" in fake_hera.documents

    def test_route_rejects_a_bad_doc_id(self, fake_hera, warmed):
        with pytest.raises(HTTPException) as raised:
            routes.delete(types.SimpleNamespace(docIds=["nope"]))

        assert raised.value.status_code == 404


# Every route reaches hera, so each one refuses while hera is still warming.
WARMING_CALLS = {
    "create": lambda: routes.create(_create()),
    "update": lambda: routes.update(types.SimpleNamespace(docId="doc1", fields={})),
    "delete": lambda: routes.delete(types.SimpleNamespace(docIds=["doc1"])),
}


@pytest.mark.parametrize("route", sorted(WARMING_CALLS))
def test_every_route_says_warming_up(route, fake_hera, monkeypatch):
    monkeypatch.setattr(dependencies.warmup, "_ready", False)

    with pytest.raises(HTTPException) as raised:
        WARMING_CALLS[route]()

    assert raised.value.status_code == 503
    assert raised.value.detail == "WARMING_UP"
