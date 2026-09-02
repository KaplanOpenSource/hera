"""datalayer/document/metadataDocument.py: nonDBMetadataFrame and
parseConnectionString (pure helpers, no DB/filesystem access), plus
MetadataFrame.__str__ against a real mongomock-backed document."""
import json

import pytest

from hera.datalayer.document import parseConnectionString
from hera.datalayer.document.metadataDocument import nonDBMetadataFrame


@pytest.mark.unit
class TestNonDBMetadataFrame:
    def test_get_data_returns_the_stored_data(self):
        frame = nonDBMetadataFrame(data={"a": 1}, type="t", resource="r")
        assert frame.getData() == {"a": 1}

    def test_getitem_reads_from_the_instance_dict(self):
        frame = nonDBMetadataFrame(data={"a": 1}, type="t", resource="r")
        assert frame["type"] == "t"
        assert frame["resource"] == "r"


@pytest.mark.unit
class TestParseConnectionString:
    def test_it_splits_credentials_host_and_dbname(self):
        result = parseConnectionString("user:pass@host.com/mydb")
        assert result == {
            "username": "user",
            "password": "pass",
            "dbName": "mydb",
            "dbIP": "host.com",
        }


@pytest.mark.unit
class TestMetadataFrameStr:
    def test_it_pretty_prints_asdict_without_the_id(self, unit_project):
        doc = unit_project.addMeasurementsDocument(
            resource="somepath", dataFormat="string", type="t", desc={"a": 1},
        )
        assert str(doc) == json.dumps(doc.asDict(with_id=False), indent=4)
        assert "_id" not in json.loads(str(doc))
