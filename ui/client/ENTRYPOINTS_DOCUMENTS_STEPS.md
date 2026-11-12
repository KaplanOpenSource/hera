# Steps: the document entrypoints

The **Documents** part of `ENTRYPOINTS_PLAN.md` (issue #1129). Three named
entrypoints replace six `/exec` call sites. Status: **built**.

Built the same way as the project entrypoints, so read
`ENTRYPOINTS_PROJECTS_STEPS.md` first: same router shape, same warmup guard,
same `getJson` / `postJson` on the client.

## The routes

| Entrypoint | Route | Body | Answer |
|---|---|---|---|
| `document.create` | `POST /document/create` | `projectName`, `collection`, `type`, `resource`, `dataFormat`, `desc`, `notebook` | `{document, project, addedDocId}` |
| `document.update` | `POST /document/update` | `docId`, `fields` | `{document, project}` |
| `document.delete` | `POST /document/delete` | `docIds` | `{project}` |

`addedDocId` is the created document's own `_id`, given at the top level so the
caller does not have to dig for it. Delete has no document to answer with.

`project` is the same `{name, documents}` shape `/project/documents` answers, so
the client puts it straight in the store. That is what turns the
write-then-reload pairs into single calls: four of the five call sites today do
a second round trip, and only `AddDocumentButton` already gets the project back
in one.

Every route reaches hera, so each answers 503 `WARMING_UP` while it warms, like
`/project/*`. A missing or unknown `docId` is a 404. An empty `docIds` list, an
unknown `collection` or a `desc` that is not a dict is a 400 with a message.

## Not CRUD's R

There is no per-document read. `/project/documents` already brings every
document of the project in one call and the store keeps them, so a second way to
fetch the same thing would only go stale next to it.

So `document.get` from `ENTRYPOINTS_PLAN.md` is dropped. Its only caller is
`PreviewPanel`, which wants `doc.resource`. The sibling `hasPreview(doc)` is
already handed the whole document from the store, and the docid comes from an
open tab, so the document is always there. `PreviewPanel` reads the store.

## No duplicate route

`ProjectDocument` in the store already holds `_cls`, `type`, `resource`,
`dataFormat` and `desc`, for every document of the project. `_cls` is
`Metadata.Measurements` / `Metadata.Simulations` / `Metadata.Cache`, the same
three names `METADATA_CLASSES` already uses.

So duplicate is create with the fields taken from the store: the client maps
`_cls` to its collection, picks a free `datasourceName` from the store's own
list, and calls `document.create`. No read, and no copy logic on the server.

`duplicateDocuments` keeps its signature, so `TreeContextMenu` does not change.
It becomes N `document.create` calls, one per selected document. The last
answer's project goes in the store.

## Create and update stay apart

Their payloads look alike but are not. Create needs `collection` and
`projectName` and uses them. Update must read them off the document and refuse
to change them. One `save` route would mean fields that count on one path and
are ignored on the other, and only create has an id to answer with.

## No `kind`

`buildAddDocumentCode` branches on four kinds (Document, Agent, Notebook,
Workflow), but three of the four differ only in the fields they write:

- Agent: `resource` is `{"effects": {}}`, type `ToolkitDataSource`.
- Workflow: the `desc.workflow` block and an empty `parameters` index, type
  `hermesWorkflow`, Simulations collection.
- Document: the collection the user chose, no type.

Those are fields, so the client fills them and `create` stays one generic route.
`kind` never crosses the wire and the server learns nothing about the Add
dialog. `DocKind` stays in `AddDocumentButton.tsx`, client-only as it is today.

The notebook is the one real exception: it writes an `.ipynb` file on disk. So
`create` takes an optional `notebook` object, `{toolkitNames}`. When it is
there, the server writes an empty notebook at `resource` if no file is there
yet, with the toolkit lines in its first cell. When it is absent nothing touches
the disk. An existing file is always left alone, as today.

## Update sends the whole document

Today `updateDocument` compares the edited document against the loaded one and
emits one python assignment per changed top-level field. The diff goes: it was a
side effect of building python text, and it never protected much anyway, since
it is per top-level field and every `desc` edit sends the whole `desc`.

`fields` is therefore the document as the editor holds it, minus `_id`, `_cls`
and `projectName`. The server refuses those three with a 400. The client keeps
its `FORBIDDEN_FIELDS` list for the editor, but it is no longer what protects
the document.

`stringifyToPython` goes with it. None/True/False were only a problem because
the values were pasted into python text. JSON carries them.

## Server files

```
ui/server/document_api/__init__.py
ui/server/document_api/document_api.py     # class DocumentApi: one document
ui/server/document_api/notebook_file.py    # class NotebookFile: write an empty notebook
ui/server/routers/documents.py             # the three routes
```

`DocumentApi` is built from a document id and dropped when the request ends, so
it holds no state, exactly like `ProjectApi`. Methods: `update`, `delete`,
`project_name`. Creating needs no document id, so it is a static `create` on the
same class, taking the project name.

`projectName` is not sent on update or delete. The server reads it off the first
document, **before** deleting anything. A list of ids spanning two projects is a
400.

`NotebookFile` is built from a path. It writes the empty notebook with its first
cell when the path does not exist, and leaves an existing file alone.

The project in every answer comes from `ProjectApi(projectName).documents()`, so
`/document/*` and `/project/documents` cannot drift apart.

New in `api_models.py`: `CreateDocumentPayload`, `UpdateDocumentPayload`,
`DocumentIdsPayload`, and a `Collection` enum with the three values. The client
gets the same enum in shared types, so both sides name them alike.

Register the router in `server.py` next to the project one.

## Client files

- `src/io/createDocument.ts` - new: one `document.create` call, shared by the
  Add dialog and duplicate.
- `src/io/FetchDocument.tsx` - `updateDocument` calls `document.update` with the
  document as JSON, and stores the project it answers.
  `stringifyToPython` and `fetchDocument` are deleted.
- `src/io/deleteDocuments.ts` - one `postJson`; stores the returned project
  instead of calling `fetchProjectDetails`. Signature unchanged, so
  `DeleteSelectedButton` and `TreeContextMenu` stay as they are.
- `src/io/duplicateDocuments.ts` - builds each copy from the store and calls
  `createDocument`. The free-name rule moves here from the python helper.
- `src/components/details/DeleteDocumentButton.tsx` - `document.delete` with one
  id, and no reload after it.
- `src/components/project/AddDocumentButton.tsx` - builds the fields per kind,
  calls `createDocument`, takes `addedDocId` from the answer instead of diffing
  the id sets before and after.
- `src/components/details/PreviewPanel.tsx` - reads the document from the store
  by docid; no fetch, so the `useEffect` and the `resource` state go.
- `src/components/details/useShownDoc.ts` - stores the project that
  `document.update` answered, so saving is one call, not two.
- Gone: `src/components/project/buildAddDocumentCode.ts`,
  `src/components/project/buildNotebookCommand.ts`.

## Tests

- `ui/server/tests/test_documents.py` - `DocumentApi`, `NotebookFile` and the
  three routes, on the fake hera in `conftest.py`. Cover: a bad docId, a
  forbidden field, an empty id list, ids from two projects, an unknown
  collection, and that an existing notebook file is not overwritten.
- Client tests assert the entrypoint name and the parameters sent, not python
  text. `addDocument.test.tsx`, `deleteDocumentButton.test.tsx` and
  `documentOperations.test.ts` already exist and need rewriting that way.
- Add a client test for the duplicate name rule, now that it is TypeScript.
- Add a client test that saving a document makes one call, not two.
- Add a client test that the preview makes no call at all.
- Update `ui/client/tests/execPython-coverage.md`: the six sites go.

## Order of work

1. `api_models.py` payloads and the collection enum, `document_api/` classes,
   `routers/documents.py`.
2. `ui/server/tests/test_documents.py`, green before any client change.
3. `PreviewPanel` off `fetchDocument`, then delete `fetchDocument`.
4. `createDocument.ts`, `AddDocumentButton`, and delete the two builder files.
5. `duplicateDocuments.ts`, `deleteDocuments.ts`, `DeleteDocumentButton`.
6. `updateDocument` and `useShownDoc`.
7. Client tests, then the full `TEST_UI.md` checklist.

## Follow-up, not part of this

Duplicating a notebook document copies `resource` as it is, so the copy points
at the same `.ipynb` file and both documents edit one notebook. That is today's
behavior and this change keeps it.

The fix is a second shape for the `notebook` parameter, `{copyFrom}`, which
copies that file to `resource` instead of writing an empty one. Duplicate would
then set `resource` from the free name it picked. It needs its own issue, so
#1129 stays a pure refactor.

## Left over

`/exec` call sites after this: 8, all of them repositories.
