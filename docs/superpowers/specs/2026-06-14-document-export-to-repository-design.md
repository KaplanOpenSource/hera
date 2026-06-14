# Design Spec — Document Export to Repository (Issue #932)

- **Branch:** `issue932`
- **Date:** 2026-06-14
- **GitHub issue:** #932 — *Implement document export to repository*

## 1. Goal & Scope

Implement the **reverse** of the existing repository loader: take Metadata
documents that live in a Hera project and write them into a **repository JSON
file** that can later be re-loaded with the existing
`loadAllDatasourcesInRepositoryJSONToProject`.

Operations required by the issue:

1. Export a single selected document to a repository.
2. Export multiple selected documents to a repository.
3. Export **all** documents of a project to a repository.
4. Add documents to an **existing** repository, checking whether they already exist.
5. Detect duplicates using a **content hash** (primary) or the document's **ObjectId** (fallback).
6. **Override** a repository to remove duplicates.

### Decisions locked during brainstorming

| Decision | Choice |
|----------|--------|
| Resource files | **Reference only** — record the existing path, `isRelativePath:"False"`. No file copying in the MVP. |
| Duplicate identity | **Content hash primary, ObjectId as selectable fallback** (`idStrategy`). |
| Top-level toolkit key | **User-supplied** `toolkitName` (required argument). |
| Deliverables | Python API on `dataToolkit` + CLI subcommand + pytest tests. |

### Out of scope (MVP)

- Copying / bundling resource data files next to the repository JSON (deferred; a
  future `copyResources=True` flag is the natural extension point).
- Inferring the toolkit name from `desc.toolkit` / `type` prefix.
- Cross-toolkit auto-grouping (everything exported in one call goes under one
  user-supplied `toolkitName`).

## 2. Architecture & Components

```
hera/utils/data/repositoryExport.py        # NEW — pure functions over dicts (no DB / no IO)
    documentContentHash(docDict, idStrategy="contentHash") -> str
    documentToRepositoryItem(docDict, idStrategy="contentHash") -> (section, itemName, entry)
    mergeDocumentsIntoRepository(repoJSON, docDicts, toolkitName,
                                 idStrategy="contentHash", overwrite=False) -> (newRepoJSON, report)
    deduplicateRepository(repoJSON) -> (newRepoJSON, report)

hera/utils/data/toolkit.py  (class dataToolkit)   # thin facade
    exportDocumentsToRepository(...)   # query DB -> pure funcs -> write file -> (optionally) addRepository

hera/bin/hera-project  +  hera/utils/data/CLI.py  # CLI: `hera-project repository export`

hera/tests/test_repository_export.py              # tests (mostly DB-free)
```

**Principle of isolation.** All decision logic (hashing, JSON-structure merge,
duplicate detection, dedup) lives in pure functions that take and return plain
dicts. They require neither MongoDB nor a toolkit instance, so they are
unit-tested in isolation and fast. `dataToolkit` only orchestrates:
`getDocuments` → pure functions → file write.

## 3. Data Model — Repository JSON Entry

A repository JSON is `{ <toolkitName>: { <Section>: { <itemName>: <entry> } } }`.
`Section` is derived from the document collection (`_cls`): `Metadata.Measurements`
→ `Measurements`, `Metadata.Simulations` → `Simulations`, `Metadata.Cache` → `Cache`.

Each exported document becomes one entry:

```json
"<itemName>": {
  "isRelativePath": "False",
  "contentHash": "9f86d081884c7d659a2feaa0c55ad015a3bf4f1b2b0b822cd15d6c15b0f00a08",
  "sourceId": "682f1b9e4d3a2c0011aa1c3",
  "item": {
    "type": "highFreqMeteorology_HighFreqData",
    "resource": "/abs/path/slicedYamim_sonic.parquet",
    "dataFormat": "parquet",
    "desc": { "deviceType": "Sonic" }
  }
}
```

Round-trip guarantees:

- `item` contains exactly `type` / `resource` / `dataFormat` / `desc` — the fields
  the existing `_DocumentHandler` (`hera/utils/data/toolkit.py:365`) consumes when
  loading documents back into a project.
- `contentHash` and `sourceId` are stored as **siblings** of `item` (not inside
  `desc`). The existing Measurements/Simulations/Cache load path reads only
  `itemDesc["item"]`, so unknown siblings are ignored on load and never pollute
  the DB.
- `itemName` = the source ObjectId string when available, otherwise the first 16
  hex chars of `contentHash`. Guarantees a unique, stable key per document.

### Content hash

`documentContentHash(docDict, idStrategy)`:

- `idStrategy="contentHash"` (default): `sha256` over `json.dumps({type, dataFormat,
  resource, desc}, sort_keys=True, default=str)`. Deterministic across projects and
  machines; identical content → identical hash. `resource` is included, so two docs
  pointing at different files are *not* collapsed.
- `idStrategy="objectId"`: identity is the document `_id`; the hash field stores the
  ObjectId string. Two documents are duplicates only when they share an `_id`.

The "fallback" relationship: identity comparison prefers `contentHash`; when an
existing repo entry lacks a stored `contentHash` (e.g. hand-authored repos), the
comparison falls back to `sourceId`/ObjectId before treating the doc as new.

## 4. Public API

### 4.1 Pure functions (`repositoryExport.py`)

```python
def documentContentHash(docDict: dict, idStrategy: str = "contentHash") -> str:
    """Stable identity string for a document dict."""

def documentToRepositoryItem(docDict: dict, idStrategy: str = "contentHash") -> tuple[str, str, dict]:
    """Return (sectionName, itemName, entry) for one document `asDict()` result.
    Raises ValueError if `_cls` is missing/unrecognised."""

def mergeDocumentsIntoRepository(repoJSON: dict, docDicts: list[dict], toolkitName: str,
                                 idStrategy: str = "contentHash", overwrite: bool = False):
    """Merge documents under repoJSON[toolkitName]. Returns (newRepoJSON, report).
    report = {added: [...], skipped_existing: [...], overwritten: [...]}.
    A document already present (same identity, any section under toolkitName) is
    skipped unless overwrite=True."""

def deduplicateRepository(repoJSON: dict):
    """Collapse entries that share the same identity to a single entry, across all
    toolkits/sections. Returns (newRepoJSON, report) with report.removed listing
    dropped (toolkit, section, itemName) keys."""
```

All four are side-effect-free and operate on deep copies; the inputs are never
mutated.

### 4.2 Facade method (`dataToolkit`)

```python
def exportDocumentsToRepository(self, *, toolkitName, repositoryName,
                                projectName=None, documents=None,
                                idStrategy="contentHash", mode="add",
                                register=True, overwrite=False):
    """
    documents:
        - None  -> export ALL documents of `projectName`
        - a single document or id  -> export that one
        - a list of documents/ids  -> export those
    repositoryName:
        name of a registered repo OR a path to a .json file. If it resolves to an
        existing file, that file is the merge base; otherwise a new repo is created.
    mode:
        "add"      -> merge (check-exists, skip dups unless overwrite)
        "override" -> after merge, run deduplicateRepository over the whole file
    register:
        if True, (re)register the resulting file via addRepository so it shows up in
        `hera-project repository list`.
    Returns the report dict from the merge/dedup step.
    """
```

Responsibilities of the facade: resolve `documents` to a list of `asDict(with_id=True)`
dicts (querying the project / fetching by id), resolve the repository file path (load
existing JSON or start `{}`), call the pure functions, write the JSON file (pretty,
`indent=2`), and optionally `addRepository`.

### 4.3 CLI

New subcommand under the existing `repository` group in `hera/bin/hera-project`,
handler `repository_export` in `hera/utils/data/CLI.py`:

```
hera-project repository export <repositoryName|path>
    --toolkitName NAME            # required: top-level toolkit key
    --projectName NAME            # source project (default project if omitted)
    --documentId ID               # repeatable; omit to export ALL
    --idStrategy contentHash|objectId   # default contentHash
    --mode add|override           # default add
    --no-register                 # skip addRepository
    --overwrite                   # overwrite existing entries on identity match
```

Mirrors the argument style of the existing `repository add` / `repository load`
subcommands (`hera/bin/hera-project:137-155`).

## 5. Data Flow

```
exportDocumentsToRepository
  ├─ resolve source documents
  │     None        -> proj.getAllDocuments()            (Measurements+Simulations+Cache)
  │     id / [ids]  -> proj.getDocumentByID(id) ...
  │     doc / [docs]-> use as-is
  │     -> [ doc.asDict(with_id=True) for doc in docs ]
  ├─ resolve repository file
  │     registered name -> getDataSourceDocument(name).resource (path)
  │     path            -> use path
  │     existing file?  -> json.load ; else {}
  ├─ mergeDocumentsIntoRepository(repoJSON, docDicts, toolkitName, idStrategy, overwrite)
  ├─ if mode == "override": deduplicateRepository(repoJSON)
  ├─ write repoJSON -> file (indent=2)
  └─ if register: addRepository(name, path, overwrite=True)
  -> report
```

## 6. Error Handling & Edge Cases

- **Unknown `_cls`** in a document dict → `documentToRepositoryItem` raises
  `ValueError`; the facade logs and skips that document (consistent with the
  existing loader's per-item skip-and-log behaviour at `toolkit.py:270-274`).
- **Embedded resource (JSON_dict/string)** → `resource` is an object/string, not a
  path. It is hashed and stored verbatim; `isRelativePath:"False"` is correct since
  there is no path to resolve.
- **Empty document set** → write an empty/unchanged repo and return a report with
  empty lists (no error).
- **`documentId` not found** → raise `ValueError` naming the id (CLI prints a clear
  error).
- **Existing repo file is invalid JSON** → raise `ValueError` with the path; do not
  silently overwrite a corrupt-but-real file.
- **Default project guard** → exporting *from* the default project is allowed
  (read-only query); only `addRepository` writes, and it already toggles
  `_allowWritingToDefaultProject` (`toolkit.py:58-63`).
- **Idempotency** → re-running the same export in `mode="add"` without `overwrite`
  is a no-op (all docs already present → all `skipped_existing`).

## 7. Testing Plan (`hera/tests/test_repository_export.py`)

DB-free unit tests (the bulk):

- `documentContentHash`: identical content → identical hash; differing
  `resource`/`desc`/`type` → different hash; `objectId` strategy uses `_id`.
- `documentToRepositoryItem`: section mapping for each `_cls`; itemName uses
  ObjectId when present, hash prefix otherwise; entry shape (siblings + `item`);
  `ValueError` on bad `_cls`.
- `mergeDocumentsIntoRepository`: new docs added; duplicate (same identity) skipped;
  `overwrite=True` replaces; report contents correct; input not mutated.
- `deduplicateRepository`: duplicates across sections collapsed to one; report lists
  removed keys; unique entries untouched.
- **Round-trip**: build a repo via the pure functions, then load it with
  `loadRepositoryFromPath` / `loadAllDatasourcesInRepositoryJSONToProject` and assert
  the `item` fields survive (uses existing test toolkits where a DB is available;
  otherwise asserts on the resolved dict).

Integration test (DB-backed, follows `test_repository.py` fixtures): create a temp
project, add a couple of documents, `exportDocumentsToRepository`, assert the file
exists, is valid JSON, contains the expected entries, and re-loads into a fresh
project.

Follows repo conventions: `Test<Thing>` classes, session/function fixtures,
`PYTEST_PROJECT_NAME`, cleanup on teardown.

## 8. Files Touched

| File | Change |
|------|--------|
| `hera/utils/data/repositoryExport.py` | **new** — pure functions |
| `hera/utils/data/toolkit.py` | add `exportDocumentsToRepository` facade method |
| `hera/utils/data/CLI.py` | add `repository_export` handler |
| `hera/bin/hera-project` | register `repository export` subcommand |
| `hera/tests/test_repository_export.py` | **new** — tests |
| `docs/...` (optional) | mention export in repository docs if present |
