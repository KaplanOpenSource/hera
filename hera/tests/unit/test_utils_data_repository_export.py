"""repositoryExport.py: pure dict/JSON helpers for exporting Metadata
documents into a repository-JSON structure. No DB, no file IO.

B91: mergeDocumentsIntoRepository appends the pre-disambiguation item_name
to report['added']/['overwritten'] before _uniqueItemName renames the key
actually written to the dict -- on a name collision between distinct
documents, the report names an item that isn't the one just added. See
TestMergeReportNameIsBrokenOnCollision.
"""
import pytest

from hera.utils.data import repositoryExport as re


def _doc(cls="Metadata.Measurements", oid="abc123", type_="myType", resource="/r", desc=None):
    return {
        "_id": {"$oid": oid} if oid else None,
        "_cls": cls,
        "type": type_,
        "resource": resource,
        "dataFormat": "parquet",
        "desc": desc or {},
    }


@pytest.mark.unit
class TestDocumentContentHash:
    def test_the_same_content_hashes_the_same(self):
        a = re.documentContentHash(_doc(oid="id1"))
        b = re.documentContentHash(_doc(oid="id2"))
        assert a == b  # _id is excluded from the hash

    def test_different_content_hashes_differently(self):
        a = re.documentContentHash(_doc(type_="A"))
        b = re.documentContentHash(_doc(type_="B"))
        assert a != b

    def test_objectid_strategy_returns_the_raw_id_string(self):
        assert re.documentContentHash(_doc(oid="myid"), idStrategy="objectId") == "myid"

    def test_objectid_strategy_raises_without_an_id(self):
        with pytest.raises(ValueError, match="no _id"):
            re.documentContentHash(_doc(oid=None), idStrategy="objectId")

    def test_an_unknown_strategy_raises(self):
        with pytest.raises(ValueError, match="Unknown idStrategy"):
            re.documentContentHash(_doc(), idStrategy="bogus")


@pytest.mark.unit
class TestDocumentToRepositoryItem:
    def test_it_maps_measurements_cls_to_the_measurements_section(self):
        section, name, entry = re.documentToRepositoryItem(_doc(cls="Metadata.Measurements"))
        assert section == "Measurements"
        assert entry["item"]["type"] == "myType"

    def test_it_maps_simulations_and_cache_sections_too(self):
        assert re.documentToRepositoryItem(_doc(cls="Metadata.Simulations"))[0] == "Simulations"
        assert re.documentToRepositoryItem(_doc(cls="Metadata.Cache"))[0] == "Cache"

    def test_a_missing_cls_raises(self):
        with pytest.raises(ValueError, match="missing or malformed"):
            re.documentToRepositoryItem(_doc(cls=None))

    def test_an_unsupported_cls_raises(self):
        with pytest.raises(ValueError, match="Unsupported"):
            re.documentToRepositoryItem(_doc(cls="Metadata.SomethingElse"))

    def test_the_item_name_prefers_the_source_object_id(self):
        _, name, _ = re.documentToRepositoryItem(_doc(oid="myid"))
        assert name == "myid"

    def test_the_item_name_falls_back_to_a_hash_prefix_without_an_id(self):
        _, name, entry = re.documentToRepositoryItem(_doc(oid=None))
        assert name == entry["contentHash"][:16]


@pytest.mark.unit
class TestMergeDocumentsIntoRepository:
    def test_a_new_document_is_added(self):
        repo, report = re.mergeDocumentsIntoRepository({}, [_doc(oid="id1")], "MyToolkit")
        assert report["added"] == ["id1"]
        assert "id1" in repo["MyToolkit"]["Measurements"]

    def test_the_input_repo_is_never_mutated(self):
        original = {}
        re.mergeDocumentsIntoRepository(original, [_doc()], "MyToolkit")
        assert original == {}

    def test_a_document_with_the_same_content_is_skipped_by_default(self):
        repo, _ = re.mergeDocumentsIntoRepository({}, [_doc(oid="id1")], "MyToolkit")
        repo2, report2 = re.mergeDocumentsIntoRepository(repo, [_doc(oid="id2")], "MyToolkit")
        assert report2["skipped_existing"] == ["id2"]
        assert "id2" not in repo2["MyToolkit"]["Measurements"]

    def test_overwrite_true_replaces_the_previous_entry(self):
        repo, _ = re.mergeDocumentsIntoRepository({}, [_doc(oid="id1", type_="old")], "MyToolkit")
        repo2, report2 = re.mergeDocumentsIntoRepository(
            repo, [_doc(oid="id2", type_="old")], "MyToolkit", overwrite=True
        )
        assert report2["overwritten"] == ["id2"]
        assert "id1" not in repo2["MyToolkit"]["Measurements"]
        assert "id2" in repo2["MyToolkit"]["Measurements"]

    def test_a_name_collision_with_different_content_is_disambiguated_in_storage(self):
        repo, _ = re.mergeDocumentsIntoRepository({}, [_doc(oid="id1", type_="A")], "MyToolkit")
        second_doc = _doc(oid="id1", type_="B")
        repo2, report2 = re.mergeDocumentsIntoRepository(repo, [second_doc], "MyToolkit")
        names = list(repo2["MyToolkit"]["Measurements"].keys())
        assert len(names) == 2
        assert names[1] != "id1"  # disambiguated in the stored dict


@pytest.mark.unit
class TestMergeReportNameIsBrokenOnCollision:
    """B91: mergeDocumentsIntoRepository appends the pre-disambiguation
    item_name to report['added']/['overwritten'] before _uniqueItemName
    renames it -- so on a name collision with distinct content, the report
    names an item_name that is not actually a key in the returned repo."""

    @pytest.mark.xfail(
        strict=True,
        reason="B91: report['added'] is recorded before _uniqueItemName "
               "disambiguates the key, so it can name an item that was "
               "never actually stored under that name. "
               "See the consolidated findings issue.",
    )
    def test_the_reported_name_should_point_to_the_just_added_document(self):
        repo, _ = re.mergeDocumentsIntoRepository({}, [_doc(oid="id1", type_="A")], "MyToolkit")
        repo2, report2 = re.mergeDocumentsIntoRepository(
            repo, [_doc(oid="id1", type_="B")], "MyToolkit"
        )
        reported_entry = repo2["MyToolkit"]["Measurements"][report2["added"][0]]
        assert reported_entry["item"]["type"] == "B"

    def test_currently_the_reported_name_belongs_to_a_different_document(self):
        """Characterisation of B91: report2['added'] == ['id1'], but the
        actual key 'id1' in repo2 is the *first* document (added by the
        earlier call); the second document actually landed at a
        disambiguated key the report never mentions."""
        repo, _ = re.mergeDocumentsIntoRepository({}, [_doc(oid="id1", type_="A")], "MyToolkit")
        repo2, report2 = re.mergeDocumentsIntoRepository(
            repo, [_doc(oid="id1", type_="B")], "MyToolkit"
        )
        assert report2["added"] == ["id1"]
        stored_under_id1 = repo2["MyToolkit"]["Measurements"]["id1"]
        assert stored_under_id1["item"]["type"] == "A"  # not the just-added "B" document


@pytest.mark.unit
class TestDeduplicateRepository:
    def test_duplicate_entries_by_content_hash_are_collapsed(self):
        repo, _ = re.mergeDocumentsIntoRepository({}, [_doc(oid="id1", type_="A")], "MyToolkit")
        repo["MyToolkit"]["Measurements"]["id2_dup"] = dict(
            repo["MyToolkit"]["Measurements"]["id1"]
        )
        deduped, report = re.deduplicateRepository(repo)
        assert len(deduped["MyToolkit"]["Measurements"]) == 1
        assert len(report["removed"]) == 1

    def test_the_input_repo_is_never_mutated(self):
        repo, _ = re.mergeDocumentsIntoRepository({}, [_doc(oid="id1")], "MyToolkit")
        repo["MyToolkit"]["Measurements"]["dup"] = dict(
            repo["MyToolkit"]["Measurements"]["id1"]
        )
        import copy

        before = copy.deepcopy(repo)
        re.deduplicateRepository(repo)
        assert repo == before

    def test_entries_with_no_identity_are_kept_as_is(self):
        repo = {"MyToolkit": {"Measurements": {"a": {"item": {}}}}}
        deduped, report = re.deduplicateRepository(repo)
        assert "a" in deduped["MyToolkit"]["Measurements"]
        assert report["removed"] == []
