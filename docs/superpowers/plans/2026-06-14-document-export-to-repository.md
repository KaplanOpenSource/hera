# Document Export to Repository — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Export Hera project Metadata documents into a repository JSON file (reference-only), with content-hash/ObjectId duplicate detection and a dedup override mode.

**Architecture:** Approach C — a pure-function logic module (`repositoryExport.py`, no DB/IO) does all the hashing, JSON-structure merging and dedup; a thin `dataToolkit.exportDocumentsToRepository` facade queries the DB and writes the file; a CLI subcommand wraps the facade. The output round-trips through the existing `loadAllDatasourcesInRepositoryJSONToProject`.

**Tech Stack:** Python 3, mongoengine (existing data layer), pytest, argparse (existing `hera-project` CLI).

**Spec:** `docs/superpowers/specs/2026-06-14-document-export-to-repository-design.md`

---

## File Structure

| File | Responsibility |
|------|----------------|
| `hera/utils/data/repositoryExport.py` | **NEW.** Pure functions over dicts: `documentContentHash`, `documentToRepositoryItem`, `mergeDocumentsIntoRepository`, `deduplicateRepository`. No DB, no file IO. |
| `hera/utils/data/toolkit.py` | **MODIFY.** Add `exportDocumentsToRepository` method to `dataToolkit` (the facade). |
| `hera/utils/data/CLI.py` | **MODIFY.** Add `repository_export` handler. |
| `hera/bin/hera-project` | **MODIFY.** Register the `repository export` subcommand. |
| `hera/tests/test_repository_export.py` | **NEW.** Unit tests (DB-free) for the pure module + DB-backed integration test for the facade. |

### Document dict shape (input to the pure functions)

`MetadataFrame.asDict(with_id=True)` returns a dict like:

```python
{
    "_cls": "Metadata.Measurements",        # or .Simulations / .Cache
    "_id": {"$oid": "682f1b9e4d3a2c0011aa1c3"},  # may also be a plain str
    "projectName": "SOME_PROJECT",
    "type": "highFreqMeteorology_HighFreqData",
    "resource": "/abs/path/file.parquet",   # str path OR embedded object (JSON_dict)
    "dataFormat": "parquet",
    "desc": {"deviceType": "Sonic"},
}
```

`section = _cls.split(".")[1]` → `"Measurements"` / `"Simulations"` / `"Cache"`.

---

## Task 1: Pure module — `documentContentHash`

**Files:**
- Create: `hera/utils/data/repositoryExport.py`
- Test: `hera/tests/test_repository_export.py`

- [ ] **Step 1: Write the failing tests**

Create `hera/tests/test_repository_export.py`:

```python
"""Unit tests for hera.utils.data.repositoryExport (pure, DB-free)."""

import copy

import pytest

from hera.utils.data import repositoryExport as rx


# ---------------------------------------------------------------------------
# Sample document dicts (shape of MetadataFrame.asDict(with_id=True))
# ---------------------------------------------------------------------------

def _doc(_cls="Metadata.Measurements", oid="682f1b9e4d3a2c0011aa1c3",
         type_="highFreqMeteorology_HighFreqData",
         resource="/abs/path/file.parquet", dataFormat="parquet", desc=None):
    return {
        "_cls": _cls,
        "_id": {"$oid": oid},
        "projectName": "SOME_PROJECT",
        "type": type_,
        "resource": resource,
        "dataFormat": dataFormat,
        "desc": {"deviceType": "Sonic"} if desc is None else desc,
    }


class TestDocumentContentHash:
    def test_identical_content_same_hash(self):
        h1 = rx.documentContentHash(_doc())
        h2 = rx.documentContentHash(_doc(oid="DIFFERENT_ID_SAME_CONTENT"))
        assert h1 == h2  # _id is NOT part of the content hash

    def test_different_resource_differs(self):
        h1 = rx.documentContentHash(_doc(resource="/a.parquet"))
        h2 = rx.documentContentHash(_doc(resource="/b.parquet"))
        assert h1 != h2

    def test_different_desc_differs(self):
        h1 = rx.documentContentHash(_doc(desc={"deviceType": "Sonic"}))
        h2 = rx.documentContentHash(_doc(desc={"deviceType": "TRH"}))
        assert h1 != h2

    def test_desc_key_order_irrelevant(self):
        h1 = rx.documentContentHash(_doc(desc={"a": 1, "b": 2}))
        h2 = rx.documentContentHash(_doc(desc={"b": 2, "a": 1}))
        assert h1 == h2

    def test_is_hex_sha256(self):
        h = rx.documentContentHash(_doc())
        assert len(h) == 64
        int(h, 16)  # raises if not hex

    def test_objectid_strategy_uses_id(self):
        h1 = rx.documentContentHash(_doc(oid="AAA"), idStrategy="objectId")
        h2 = rx.documentContentHash(_doc(oid="BBB"), idStrategy="objectId")
        assert h1 == "AAA"
        assert h2 == "BBB"

    def test_objectid_strategy_plain_string_id(self):
        d = _doc()
        d["_id"] = "PLAIN_STRING_ID"
        assert rx.documentContentHash(d, idStrategy="objectId") == "PLAIN_STRING_ID"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest hera/tests/test_repository_export.py::TestDocumentContentHash -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'hera.utils.data.repositoryExport'`

- [ ] **Step 3: Write minimal implementation**

Create `hera/utils/data/repositoryExport.py`:

```python
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
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest hera/tests/test_repository_export.py::TestDocumentContentHash -v`
Expected: PASS (7 passed)

- [ ] **Step 5: Commit**

```bash
git add hera/utils/data/repositoryExport.py hera/tests/test_repository_export.py
git commit -m "feat: add documentContentHash for repository export (#932)"
```

---

## Task 2: Pure module — `documentToRepositoryItem`

**Files:**
- Modify: `hera/utils/data/repositoryExport.py`
- Test: `hera/tests/test_repository_export.py`

- [ ] **Step 1: Write the failing tests**

Append to `hera/tests/test_repository_export.py`:

```python
class TestDocumentToRepositoryItem:
    def test_section_mapping_measurements(self):
        section, name, entry = rx.documentToRepositoryItem(_doc("Metadata.Measurements"))
        assert section == "Measurements"

    def test_section_mapping_simulations(self):
        section, _, _ = rx.documentToRepositoryItem(_doc("Metadata.Simulations"))
        assert section == "Simulations"

    def test_section_mapping_cache(self):
        section, _, _ = rx.documentToRepositoryItem(_doc("Metadata.Cache"))
        assert section == "Cache"

    def test_itemname_is_objectid_when_present(self):
        _, name, _ = rx.documentToRepositoryItem(_doc(oid="682f1b9e4d3a2c0011aa1c3"))
        assert name == "682f1b9e4d3a2c0011aa1c3"

    def test_itemname_is_hash_prefix_without_id(self):
        d = _doc()
        del d["_id"]
        _, name, entry = rx.documentToRepositoryItem(d)
        assert name == entry["contentHash"][:16]

    def test_entry_shape(self):
        _, _, entry = rx.documentToRepositoryItem(_doc())
        assert entry["isRelativePath"] == "False"
        assert set(entry["item"].keys()) == {"type", "resource", "dataFormat", "desc"}
        assert entry["item"]["dataFormat"] == "parquet"
        assert "contentHash" in entry
        assert entry["sourceId"] == "682f1b9e4d3a2c0011aa1c3"

    def test_unknown_cls_raises(self):
        with pytest.raises(ValueError):
            rx.documentToRepositoryItem(_doc("Metadata.Bogus"))

    def test_missing_cls_raises(self):
        d = _doc()
        del d["_cls"]
        with pytest.raises(ValueError):
            rx.documentToRepositoryItem(d)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest hera/tests/test_repository_export.py::TestDocumentToRepositoryItem -v`
Expected: FAIL — `AttributeError: module ... has no attribute 'documentToRepositoryItem'`

- [ ] **Step 3: Write minimal implementation**

Append to `hera/utils/data/repositoryExport.py`:

```python
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
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest hera/tests/test_repository_export.py::TestDocumentToRepositoryItem -v`
Expected: PASS (8 passed)

- [ ] **Step 5: Commit**

```bash
git add hera/utils/data/repositoryExport.py hera/tests/test_repository_export.py
git commit -m "feat: add documentToRepositoryItem for repository export (#932)"
```

---

## Task 3: Pure module — `mergeDocumentsIntoRepository`

**Files:**
- Modify: `hera/utils/data/repositoryExport.py`
- Test: `hera/tests/test_repository_export.py`

- [ ] **Step 1: Write the failing tests**

Append to `hera/tests/test_repository_export.py`:

```python
class TestMergeDocumentsIntoRepository:
    def test_add_to_empty_repo(self):
        repo, report = rx.mergeDocumentsIntoRepository({}, [_doc()], "MeteoHighFreq")
        assert "MeteoHighFreq" in repo
        assert "Measurements" in repo["MeteoHighFreq"]
        assert len(repo["MeteoHighFreq"]["Measurements"]) == 1
        assert len(report["added"]) == 1
        assert report["skipped_existing"] == []

    def test_duplicate_skipped(self):
        repo, _ = rx.mergeDocumentsIntoRepository({}, [_doc()], "TK")
        repo2, report = rx.mergeDocumentsIntoRepository(repo, [_doc()], "TK")
        # same content hash -> not added again
        assert len(repo2["TK"]["Measurements"]) == 1
        assert len(report["skipped_existing"]) == 1
        assert report["added"] == []

    def test_overwrite_replaces(self):
        repo, _ = rx.mergeDocumentsIntoRepository({}, [_doc()], "TK")
        repo2, report = rx.mergeDocumentsIntoRepository(
            repo, [_doc()], "TK", overwrite=True
        )
        assert len(repo2["TK"]["Measurements"]) == 1
        assert len(report["overwritten"]) == 1

    def test_distinct_docs_both_added(self):
        docs = [_doc(resource="/a.parquet"), _doc(resource="/b.parquet")]
        repo, report = rx.mergeDocumentsIntoRepository({}, docs, "TK")
        assert len(repo["TK"]["Measurements"]) == 2
        assert len(report["added"]) == 2

    def test_dup_detected_across_sections(self):
        # Same identity but a Simulations doc already present under the toolkit.
        sim = _doc("Metadata.Simulations")
        repo, _ = rx.mergeDocumentsIntoRepository({}, [sim], "TK")
        # A Measurements doc whose content hash matches must still be detected.
        same = _doc("Metadata.Measurements")
        # Force identical identity by matching all hashed fields:
        same.update({k: sim[k] for k in ("type", "resource", "dataFormat", "desc")})
        repo2, report = rx.mergeDocumentsIntoRepository(repo, [same], "TK")
        assert len(report["skipped_existing"]) == 1

    def test_input_not_mutated(self):
        original = {}
        rx.mergeDocumentsIntoRepository(original, [_doc()], "TK")
        assert original == {}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest hera/tests/test_repository_export.py::TestMergeDocumentsIntoRepository -v`
Expected: FAIL — `AttributeError: ... 'mergeDocumentsIntoRepository'`

- [ ] **Step 3: Write minimal implementation**

Append to `hera/utils/data/repositoryExport.py`:

```python
def _existingIdentities(toolkitDict):
    """Return {identity: (section, itemName)} for all entries under a toolkit.

    Identity prefers the stored ``contentHash``; falls back to ``sourceId`` for
    entries (e.g. hand-authored) that lack a hash.
    """
    index = {}
    for section, entries in toolkitDict.items():
        if not isinstance(entries, dict):
            continue
        for item_name, entry in entries.items():
            if not isinstance(entry, dict):
                continue
            identity = entry.get("contentHash") or entry.get("sourceId")
            if identity is not None:
                index[identity] = (section, item_name)
    return index


def mergeDocumentsIntoRepository(repoJSON, docDicts, toolkitName,
                                 idStrategy="contentHash", overwrite=False):
    """Merge document dicts under repoJSON[toolkitName].

    Returns (newRepoJSON, report). ``report`` has keys ``added``,
    ``skipped_existing`` and ``overwritten``, each a list of itemNames. A
    document already present (matching identity anywhere under ``toolkitName``)
    is skipped unless ``overwrite=True``. The input ``repoJSON`` is never
    mutated.
    """
    repo = copy.deepcopy(repoJSON)
    toolkitDict = repo.setdefault(toolkitName, {})
    index = _existingIdentities(toolkitDict)
    report = {"added": [], "skipped_existing": [], "overwritten": []}

    for docDict in docDicts:
        section, item_name, entry = documentToRepositoryItem(docDict, idStrategy=idStrategy)
        identity = entry.get("contentHash") or entry.get("sourceId")

        if identity in index:
            if not overwrite:
                report["skipped_existing"].append(item_name)
                continue
            # Remove the previously stored entry (possibly in another section).
            old_section, old_name = index[identity]
            toolkitDict.get(old_section, {}).pop(old_name, None)
            report["overwritten"].append(item_name)
        else:
            report["added"].append(item_name)

        toolkitDict.setdefault(section, {})[item_name] = entry
        index[identity] = (section, item_name)

    return repo, report
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest hera/tests/test_repository_export.py::TestMergeDocumentsIntoRepository -v`
Expected: PASS (6 passed)

- [ ] **Step 5: Commit**

```bash
git add hera/utils/data/repositoryExport.py hera/tests/test_repository_export.py
git commit -m "feat: add mergeDocumentsIntoRepository for repository export (#932)"
```

---

## Task 4: Pure module — `deduplicateRepository`

**Files:**
- Modify: `hera/utils/data/repositoryExport.py`
- Test: `hera/tests/test_repository_export.py`

- [ ] **Step 1: Write the failing tests**

Append to `hera/tests/test_repository_export.py`:

```python
class TestDeduplicateRepository:
    def _repo_with_dup(self):
        # Two entries, same contentHash, different itemNames, same section.
        entry_a = {"isRelativePath": "False", "contentHash": "HHH", "sourceId": "A",
                   "item": {"type": "t", "resource": "/x", "dataFormat": "parquet", "desc": {}}}
        entry_b = {"isRelativePath": "False", "contentHash": "HHH", "sourceId": "B",
                   "item": {"type": "t", "resource": "/x", "dataFormat": "parquet", "desc": {}}}
        return {"TK": {"Measurements": {"A": entry_a, "B": entry_b}}}

    def test_collapses_duplicates(self):
        repo, report = rx.deduplicateRepository(self._repo_with_dup())
        assert len(repo["TK"]["Measurements"]) == 1
        assert len(report["removed"]) == 1

    def test_unique_entries_untouched(self):
        entry_a = {"isRelativePath": "False", "contentHash": "H1", "sourceId": "A",
                   "item": {"type": "t", "resource": "/x", "dataFormat": "parquet", "desc": {}}}
        entry_b = {"isRelativePath": "False", "contentHash": "H2", "sourceId": "B",
                   "item": {"type": "t", "resource": "/y", "dataFormat": "parquet", "desc": {}}}
        repo, report = rx.deduplicateRepository({"TK": {"Measurements": {"A": entry_a, "B": entry_b}}})
        assert len(repo["TK"]["Measurements"]) == 2
        assert report["removed"] == []

    def test_dedup_across_sections(self):
        entry_a = {"isRelativePath": "False", "contentHash": "HHH", "sourceId": "A",
                   "item": {"type": "t", "resource": "/x", "dataFormat": "parquet", "desc": {}}}
        entry_b = {"isRelativePath": "False", "contentHash": "HHH", "sourceId": "B",
                   "item": {"type": "t", "resource": "/x", "dataFormat": "parquet", "desc": {}}}
        repo, report = rx.deduplicateRepository(
            {"TK": {"Measurements": {"A": entry_a}, "Simulations": {"B": entry_b}}}
        )
        total = sum(len(s) for s in repo["TK"].values())
        assert total == 1
        assert len(report["removed"]) == 1

    def test_input_not_mutated(self):
        original = self._repo_with_dup()
        snapshot = copy.deepcopy(original)
        rx.deduplicateRepository(original)
        assert original == snapshot
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest hera/tests/test_repository_export.py::TestDeduplicateRepository -v`
Expected: FAIL — `AttributeError: ... 'deduplicateRepository'`

- [ ] **Step 3: Write minimal implementation**

Append to `hera/utils/data/repositoryExport.py`:

```python
def deduplicateRepository(repoJSON):
    """Collapse entries sharing the same identity to a single entry.

    Identity is ``contentHash`` (or ``sourceId`` fallback). The first occurrence
    (iterating toolkits → sections → items) is kept; later duplicates are
    removed. Returns (newRepoJSON, report) where ``report['removed']`` lists
    (toolkit, section, itemName) tuples. Input is not mutated.
    """
    repo = copy.deepcopy(repoJSON)
    report = {"removed": []}

    for toolkitName, toolkitDict in repo.items():
        if not isinstance(toolkitDict, dict):
            continue
        seen = set()
        for section, entries in toolkitDict.items():
            if not isinstance(entries, dict):
                continue
            for item_name in list(entries.keys()):
                entry = entries[item_name]
                if not isinstance(entry, dict):
                    continue
                identity = entry.get("contentHash") or entry.get("sourceId")
                if identity is None:
                    continue
                if identity in seen:
                    del entries[item_name]
                    report["removed"].append((toolkitName, section, item_name))
                else:
                    seen.add(identity)

    return repo, report
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest hera/tests/test_repository_export.py::TestDeduplicateRepository -v`
Expected: PASS (4 passed)

- [ ] **Step 5: Run the whole pure-module suite + commit**

Run: `python -m pytest hera/tests/test_repository_export.py -v`
Expected: PASS (all unit tests so far)

```bash
git add hera/utils/data/repositoryExport.py hera/tests/test_repository_export.py
git commit -m "feat: add deduplicateRepository for repository export (#932)"
```

---

## Task 5: Facade — `dataToolkit.exportDocumentsToRepository`

**Files:**
- Modify: `hera/utils/data/toolkit.py` (add method to `dataToolkit`; add `import json` is already present, add `from hera.utils.data import repositoryExport` lazily inside the method to avoid import cycles)
- Test: `hera/tests/test_repository_export.py`

- [ ] **Step 1: Write the failing integration test**

Append to `hera/tests/test_repository_export.py`:

```python
import json
import os
import tempfile

from hera.datalayer.project import Project
from hera.utils.data.toolkit import dataToolkit

EXPORT_TEST_PROJECT = "PYTEST_EXPORT_PROJECT"


@pytest.fixture(scope="module")
def export_project():
    """A temp project holding two measurements documents."""
    files_tmp = tempfile.mkdtemp(prefix="hera_export_test_")
    proj = Project(projectName=EXPORT_TEST_PROJECT)
    proj._FilesDirectory = files_tmp
    proj.setConfig(filesDirectory=files_tmp)
    proj.addMeasurementsDocument(
        resource="/abs/path/a.parquet", dataFormat="parquet",
        type="exportTest_Data", desc={"deviceType": "Sonic"},
    )
    proj.addMeasurementsDocument(
        resource="/abs/path/b.parquet", dataFormat="parquet",
        type="exportTest_Data", desc={"deviceType": "TRH"},
    )
    yield proj
    for doc in proj.getMeasurementsDocuments():
        doc.delete()


class TestExportFacade:
    def test_export_all_writes_valid_repo(self, export_project):
        tk = dataToolkit()
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "exported_repo.json")
            report = tk.exportDocumentsToRepository(
                toolkitName="ExportTK",
                repositoryName=path,
                projectName=EXPORT_TEST_PROJECT,
                register=False,
            )
            assert os.path.isfile(path)
            with open(path) as fh:
                repo = json.load(fh)
            assert "ExportTK" in repo
            assert len(repo["ExportTK"]["Measurements"]) == 2
            assert len(report["added"]) == 2

    def test_export_is_idempotent(self, export_project):
        tk = dataToolkit()
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "repo.json")
            tk.exportDocumentsToRepository(
                toolkitName="ExportTK", repositoryName=path,
                projectName=EXPORT_TEST_PROJECT, register=False,
            )
            report = tk.exportDocumentsToRepository(
                toolkitName="ExportTK", repositoryName=path,
                projectName=EXPORT_TEST_PROJECT, register=False,
            )
            assert report["added"] == []
            assert len(report["skipped_existing"]) == 2
```

This test requires a MongoDB connection (same as `test_repository.py`); it
skips automatically in environments where the DB is unavailable because the
`Project` constructor / queries will raise and pytest collection in this repo
already tolerates that pattern. Run it where Mongo is available.

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest hera/tests/test_repository_export.py::TestExportFacade -v`
Expected: FAIL — `AttributeError: 'dataToolkit' object has no attribute 'exportDocumentsToRepository'`

- [ ] **Step 3: Write the implementation**

Add this method inside the `dataToolkit` class in `hera/utils/data/toolkit.py` (place it just after `loadRepositoryFromPath`, near the end of the class):

```python
    def exportDocumentsToRepository(self, *, toolkitName, repositoryName,
                                    projectName=None, documents=None,
                                    idStrategy="contentHash", mode="add",
                                    register=True, overwrite=False):
        """
        Export project documents into a repository JSON file.

        Parameters
        ----------
        toolkitName : str
            Top-level toolkit key under which the documents are written.
        repositoryName : str
            A registered repository name OR a path to a ``.json`` file. If it
            resolves to an existing file, that file is the merge base.
        projectName : str, optional
            Source project. Defaults to the toolkit's own project.
        documents : None | doc | id | list
            None  -> export ALL documents of the project.
            A single document/id, or a list of documents/ids -> export those.
        idStrategy : {"contentHash", "objectId"}
            Duplicate-identity strategy.
        mode : {"add", "override"}
            "add" merges (skipping duplicates); "override" additionally runs a
            full deduplication pass over the resulting file.
        register : bool
            If True, register the resulting file via ``addRepository``.
        overwrite : bool
            On identity match in "add" mode, replace the existing entry.

        Returns
        -------
        dict
            The merge/dedup report.
        """
        from hera.utils.data import repositoryExport
        logger = get_classMethod_logger(self, "exportDocumentsToRepository")

        srcProjectName = projectName or self.projectName
        proj = Project(projectName=srcProjectName)

        # 1) Resolve documents -> list of asDict(with_id=True)
        docObjs = self._resolveDocumentsForExport(proj, documents)
        docDicts = [d.asDict(with_id=True) for d in docObjs]
        logger.info(f"Exporting {len(docDicts)} documents from project {srcProjectName}")

        # 2) Resolve repository file path + load existing JSON (merge base)
        repoPath = self._resolveRepositoryPath(repositoryName)
        if os.path.isfile(repoPath):
            try:
                with open(repoPath, encoding="utf-8") as fh:
                    repoJSON = json.load(fh)
            except json.JSONDecodeError as exc:
                raise ValueError(f"Existing repository file is not valid JSON: {repoPath} ({exc})")
        else:
            repoJSON = {}

        # 3) Merge (+ optional override dedup)
        repoJSON, report = repositoryExport.mergeDocumentsIntoRepository(
            repoJSON, docDicts, toolkitName, idStrategy=idStrategy, overwrite=overwrite
        )
        if mode == "override":
            repoJSON, dedupReport = repositoryExport.deduplicateRepository(repoJSON)
            report["deduplicated"] = dedupReport["removed"]

        # 4) Write the file
        os.makedirs(os.path.dirname(os.path.abspath(repoPath)), exist_ok=True)
        with open(repoPath, "w", encoding="utf-8") as fh:
            json.dump(repoJSON, fh, indent=2)
        logger.info(f"Wrote repository to {repoPath}")

        # 5) Optionally register
        if register:
            repoRegName = os.path.basename(repoPath).split(".")[0]
            self.addRepository(repositoryName=repoRegName, repositoryPath=repoPath, overwrite=True)

        return report

    def _resolveDocumentsForExport(self, proj, documents):
        """Normalise the ``documents`` argument to a list of document objects."""
        if documents is None:
            return proj.getAllDocuments()
        if not isinstance(documents, (list, tuple)):
            documents = [documents]
        resolved = []
        for d in documents:
            if isinstance(d, str):
                doc = proj.getDocumentByID(d)
                if doc is None:
                    raise ValueError(f"Document id not found in project: {d}")
                resolved.append(doc)
            else:
                resolved.append(d)
        return resolved

    def _resolveRepositoryPath(self, repositoryName):
        """Return a filesystem path for a registered repo name or a path string."""
        if repositoryName.endswith(".json") or os.path.sep in repositoryName:
            return os.path.abspath(repositoryName)
        try:
            doc = self.getDataSourceDocument(repositoryName)
        except Exception:
            doc = None
        if doc is not None and getattr(doc, "resource", None):
            return os.path.abspath(doc.resource)
        # Unknown registered name and not a path: create alongside cwd.
        return os.path.abspath(f"{repositoryName}.json")
```

Confirm the imports `json`, `os`, and `Project` are available at the top of
`hera/utils/data/toolkit.py`. `json` and `os` are already imported; add
`from hera.datalayer import Project` near the other imports if it is not present
(the module already imports from `hera.toolkit`; verify before adding).

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest hera/tests/test_repository_export.py::TestExportFacade -v`
Expected: PASS (2 passed) where Mongo is available.

- [ ] **Step 5: Commit**

```bash
git add hera/utils/data/toolkit.py hera/tests/test_repository_export.py
git commit -m "feat: add exportDocumentsToRepository facade (#932)"
```

---

## Task 6: CLI — `repository export`

**Files:**
- Modify: `hera/utils/data/CLI.py` (add `repository_export` handler)
- Modify: `hera/bin/hera-project` (register subcommand)

- [ ] **Step 1: Add the CLI handler**

Add to `hera/utils/data/CLI.py`, immediately after `repository_load` (around line 276):

```python
def repository_export(arguments):
    """
    Export project documents into a repository JSON file.
    """
    logger = logging.getLogger("hera.bin.repository_export")
    dtk = dataToolkit()

    documents = arguments.documentId if getattr(arguments, "documentId", None) else None
    projectName = getattr(arguments, "projectName", None)

    logger.info(
        f"Exporting documents from project {projectName} to repository "
        f"{arguments.repositoryName} under toolkit {arguments.toolkitName}"
    )
    report = dtk.exportDocumentsToRepository(
        toolkitName=arguments.toolkitName,
        repositoryName=arguments.repositoryName,
        projectName=projectName,
        documents=documents,
        idStrategy=arguments.idStrategy,
        mode=arguments.mode,
        register=not arguments.no_register,
        overwrite=arguments.overwrite,
    )
    print(
        f"Export complete: {len(report['added'])} added, "
        f"{len(report['skipped_existing'])} skipped, "
        f"{len(report['overwritten'])} overwritten."
    )
    if "deduplicated" in report:
        print(f"Deduplicated: {len(report['deduplicated'])} duplicate entries removed.")
```

- [ ] **Step 2: Register the subcommand**

In `hera/bin/hera-project`, immediately after the `repository load` block
(which ends near line 155), add:

```python
    # repository export
    repository_export = repository_subparsers.add_parser('export', help='Export project documents into a repository JSON file')
    repository_export.add_argument('repositoryName', type=str, help='Registered repository name or path to the output .json file')
    repository_export.add_argument('--toolkitName', type=str, required=True, help='Top-level toolkit key for the exported documents')
    repository_export.add_argument('--projectName', type=str, default=None, help='Source project name (default project if omitted)')
    repository_export.add_argument('--documentId', action='append', default=None, help='Document id to export (repeatable). Omit to export ALL documents.')
    repository_export.add_argument('--idStrategy', choices=['contentHash', 'objectId'], default='contentHash', help='Duplicate-identity strategy')
    repository_export.add_argument('--mode', choices=['add', 'override'], default='add', help='add = merge; override = merge then deduplicate the whole file')
    repository_export.add_argument('--no-register', dest='no_register', action='store_true', default=False, help='Do not register the resulting repository file')
    repository_export.add_argument('--overwrite', dest='overwrite', action='store_true', default=False, help='Overwrite existing entries on identity match')
    repository_export.set_defaults(func=CLI.repository_export)
```

- [ ] **Step 3: Smoke-test the CLI wiring**

Run: `python -c "import ast; ast.parse(open('hera/bin/hera-project').read()); print('hera-project parses')"`
Expected: `hera-project parses`

Run: `python -c "from hera.utils.data import CLI; assert hasattr(CLI, 'repository_export'); print('handler present')"`
Expected: `handler present`

- [ ] **Step 4: Commit**

```bash
git add hera/utils/data/CLI.py hera/bin/hera-project
git commit -m "feat: add 'hera-project repository export' CLI subcommand (#932)"
```

---

## Task 7: Round-trip test + full suite

**Files:**
- Test: `hera/tests/test_repository_export.py`

- [ ] **Step 1: Write the round-trip test**

Append to `hera/tests/test_repository_export.py`:

```python
class TestRoundTrip:
    def test_exported_items_match_loadRepositoryFromPath(self, export_project):
        """An exported repo file loads back with item fields intact."""
        tk = dataToolkit()
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "rt_repo.json")
            tk.exportDocumentsToRepository(
                toolkitName="ExportTK", repositoryName=path,
                projectName=EXPORT_TEST_PROJECT, register=False,
            )
            resolved = dataToolkit.loadRepositoryFromPath(path)
            assert "ExportTK" in resolved
            items = resolved["ExportTK"]["Measurements"]
            assert len(items) == 2
            for entry in items.values():
                assert set(entry["item"].keys()) == {"type", "resource", "dataFormat", "desc"}
```

- [ ] **Step 2: Run the full new suite**

Run: `python -m pytest hera/tests/test_repository_export.py -v`
Expected: PASS (all pure-module tests; DB-backed tests pass where Mongo is available)

- [ ] **Step 3: Run the existing repository suite for regressions**

Run: `python -m pytest hera/tests/test_repository.py -v`
Expected: PASS / unchanged from baseline (skips allowed where Mongo absent)

- [ ] **Step 4: Commit**

```bash
git add hera/tests/test_repository_export.py
git commit -m "test: round-trip test for repository export (#932)"
```

---

## Self-Review

**Spec coverage:**
- Export single / multiple / all → Task 5 `_resolveDocumentsForExport` (id, list, None). ✅
- Add to existing repo + check-exists → Task 3 `mergeDocumentsIntoRepository` (load base + skip dups). ✅
- Duplicate detection content-hash/ObjectId → Task 1 `documentContentHash` (`idStrategy`). ✅
- Override to remove duplicates → Task 4 `deduplicateRepository` + Task 5 `mode="override"`. ✅
- Reference-only resource handling → Task 2 entry (`isRelativePath:"False"`, raw resource). ✅
- CLI + tests → Tasks 6, 1–4, 5, 7. ✅

**Placeholder scan:** No TBD/TODO; every code step shows full code. ✅

**Type/name consistency:** `documentContentHash`, `documentToRepositoryItem`, `mergeDocumentsIntoRepository`, `deduplicateRepository`, `exportDocumentsToRepository`, `_resolveDocumentsForExport`, `_resolveRepositoryPath` — names identical across tasks. Report keys `added`/`skipped_existing`/`overwritten`/`deduplicated` consistent between Task 3, 4, 5, 6. ✅
