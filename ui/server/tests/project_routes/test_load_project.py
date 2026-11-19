"""Loading a project from a zip."""
import io

import pytest
from fastapi import HTTPException

from load_helpers import _load, _zip_of
from project_api.project_api import ProjectApi


class TestLoadProject:
    def test_reads_the_name_from_the_first_record(self, fake_hera):
        fake_hera([])
        data = _zip_of([], [{"projectName": "Alpha", "type": "Alpha__config__"}])

        assert ProjectApi.name_in_zip(io.BytesIO(data)) == "Alpha"

    def test_no_name_in_what_is_not_a_zip(self, fake_hera):
        fake_hera([])

        assert ProjectApi.name_in_zip(io.BytesIO(b"not a zip")) is None

    def test_route_loads_under_the_exported_name(self, fake_hera, warmed):
        db = fake_hera(["Beta"], {"Beta": []})
        data = _zip_of([{"projectName": "Alpha", "type": "Alpha__config__"},
                        {"projectName": "Alpha", "type": "notebook"}])

        answer = _load(data)

        assert db.loads == [{"projectName": "Alpha", "is_hard_import": False, "show_progressbar": False}]
        assert [p["name"] for p in answer["projects"]] == ["Alpha", "Beta"]
        assert answer["project"]["name"] == "Alpha"
        assert [d["type"] for d in answer["project"]["documents"]] == ["Alpha__config__", "notebook"]

    def test_route_rejects_what_is_not_a_zip(self, fake_hera, warmed):
        db = fake_hera([])

        with pytest.raises(HTTPException) as raised:
            _load(b"not a zip")

        assert raised.value.status_code == 400
        assert db.loads == []

    def test_route_rejects_a_zip_with_no_records(self, fake_hera, warmed):
        db = fake_hera([])

        with pytest.raises(HTTPException) as raised:
            _load(_zip_of())

        assert raised.value.status_code == 400
        assert db.loads == []

    def test_route_rejects_a_name_that_is_taken(self, fake_hera, warmed):
        db = fake_hera(["Alpha"], {"Alpha": []})

        with pytest.raises(HTTPException) as raised:
            _load(_zip_of([{"projectName": "Alpha", "type": "Alpha__config__"}]))

        assert raised.value.status_code == 409
        assert db.loads == []
