"""An in-memory stand-in for the bit of hera the project entrypoints touch."""
import json
import types
import zipfile


class FakeDatabase:
    """The little bit of hera the project entrypoints touch, in memory."""

    def __init__(self, projects, documents=None, filesDirectories=None):
        self.projects = list(projects)
        self.documents = {name: list(docs) for name, docs in (documents or {}).items()}
        self.filesDirectories = dict(filesDirectories or {})
        self.created = []
        self.deletedIds = []
        self.exports = []
        self.loads = []

    # --- hera.datalayer.All ---

    def getDocumentsAsDict(self, projectName, with_id=False):
        return {"documents": self.documents.get(projectName, [])}

    def getDocuments(self, projectName):
        docs = self.documents.get(projectName, [])
        return types.SimpleNamespace(count=lambda: len(docs))

    def deleteDocumentByID(self, docId):
        self.deletedIds.append(docId)
        for name, docs in self.documents.items():
            gone = [d for d in docs if self._idOf(d) == docId]
            self.documents[name] = [d for d in docs if self._idOf(d) != docId]
            # Losing the config document is what makes a project stop existing.
            if any(d.get("type") == f"{name}__config__" for d in gone) and name in self.projects:
                self.projects.remove(name)

    @staticmethod
    def _idOf(doc):
        return doc.get("_id", {}).get("$oid")

    # --- hera.datalayer.project ---

    def getProjectList(self):
        return list(self.projects)

    def Project(self, projectName, filesDirectory=None):
        if filesDirectory:
            self.filesDirectories[projectName] = filesDirectory
        return types.SimpleNamespace(
            filesDirectory=self.filesDirectories.get(projectName, ""),
            export=lambda path, show_progressbar=True: self._export(projectName, path, show_progressbar),
            setConfig=lambda **kwargs: self._set_config(projectName, kwargs))

    def _set_config(self, project_name, kwargs):
        if "filesDirectory" in kwargs:
            self.filesDirectories[project_name] = kwargs["filesDirectory"]

    def _export(self, projectName, path, show_progressbar):
        # Like hera: one zip entry, written to whatever path or file object it is given.
        self.exports.append({"projectName": projectName, "show_progressbar": show_progressbar})
        with zipfile.ZipFile(path, mode="w") as zf:
            zf.writestr("chunk_0", json.dumps(self.documents.get(projectName, [])))

    @staticmethod
    def _iter_pickled_docs(zf, return_batched):
        # Reads the zips _export and _zip_of write: each entry a JSON list of records.
        for name in zf.namelist():
            batch = json.loads(zf.read(name))
            if return_batched:
                yield batch
            else:
                yield from batch

    def load(self, proj, exported_project_path, is_hard_import, show_progressbar=True):
        self.loads.append({"projectName": proj, "is_hard_import": is_hard_import,
                           "show_progressbar": show_progressbar})
        with zipfile.ZipFile(exported_project_path) as zf:
            records = list(self._iter_pickled_docs(zf, return_batched=False))
        self.projects.append(proj)
        self.documents[proj] = [dict(r, projectName=proj) for r in records]

    # --- hera.datalayer.Cache_Collection ---

    def Cache_Collection(self):
        return types.SimpleNamespace(
            getDocuments=self._notebookDocuments,
            addDocument=self._addDocument)

    def _notebookDocuments(self, projectName, type):
        docs = [d for d in self.documents.get(projectName, []) if d.get("type") == type]
        return [types.SimpleNamespace(resource=d.get("resource"), desc=d.get("desc")) for d in docs]

    def _addDocument(self, projectName, resource, dataFormat, type, desc):
        doc = {
            "_id": {"$oid": f"id{len(self.documents.get(projectName, []))}"},
            "resource": resource,
            "dataFormat": dataFormat,
            "type": type,
            "desc": desc,
        }
        self.documents.setdefault(projectName, []).append(doc)

    # --- hera.utils.data.CLI ---

    def project_create(self, args):
        self.created.append(args)
        self.projects.append(args.projectName)
        self.filesDirectories[args.projectName] = args.directory
