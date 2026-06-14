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
