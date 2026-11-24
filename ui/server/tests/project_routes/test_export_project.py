"""Exporting a project as a zip."""
import io
import json
import types
import zipfile

import pytest
from fastapi import HTTPException

from project_api.project_api import DEFAULT_PROJECT, ProjectApi
from routers import projects as routes


class TestExportProject:
    def test_answers_the_zip_hera_writes(self, fake_hera):
        db = fake_hera(["Alpha"], {"Alpha": [{"x": 1}]})

        data = ProjectApi("Alpha").export()

        with zipfile.ZipFile(io.BytesIO(data)) as zf:
            assert json.loads(zf.read("chunk_0")) == [{"x": 1}]
        assert db.exports == [{"projectName": "Alpha", "show_progressbar": False}]

    def test_route_answers_the_zip_as_a_download(self, fake_hera, warmed):
        fake_hera(["Alpha"], {"Alpha": [{"x": 1}]})

        answer = routes.export(types.SimpleNamespace(projectName="Alpha"))

        assert answer.media_type == "application/octet-stream"
        assert answer.headers["content-disposition"] == 'attachment; filename="Alpha.zip"'
        with zipfile.ZipFile(io.BytesIO(answer.body)) as zf:
            assert json.loads(zf.read("chunk_0")) == [{"x": 1}]

    def test_route_rejects_the_default_project(self, fake_hera, warmed):
        db = fake_hera([DEFAULT_PROJECT], {DEFAULT_PROJECT: []})

        with pytest.raises(HTTPException) as raised:
            routes.export(types.SimpleNamespace(projectName=DEFAULT_PROJECT))

        assert raised.value.status_code == 400
        assert db.exports == []

    def test_route_rejects_an_unknown_project(self, fake_hera, warmed):
        fake_hera(["Alpha"], {"Alpha": []})

        with pytest.raises(HTTPException) as raised:
            routes.export(types.SimpleNamespace(projectName="Nope"))

        assert raised.value.status_code == 404
