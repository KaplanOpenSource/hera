"""One project's documents."""
import pytest
from fastapi import HTTPException

import dependencies
from project_api.project_api import ProjectApi
from routers import projects as routes


class TestProjectDocuments:
    def test_answers_the_name_and_the_documents(self, fake_hera):
        fake_hera(["Alpha"], {"Alpha": [{"x": 1}]})

        assert ProjectApi("Alpha").documents() == {"name": "Alpha", "documents": [{"x": 1}]}

    def test_route_rejects_an_unknown_project(self, fake_hera, warmed):
        fake_hera(["Alpha"], {"Alpha": []})

        with pytest.raises(HTTPException) as raised:
            routes.documents("Nope")

        assert raised.value.status_code == 404

    def test_route_rejects_an_empty_name(self, fake_hera, warmed):
        fake_hera(["Alpha"], {"Alpha": []})

        with pytest.raises(HTTPException) as raised:
            routes.documents("")

        assert raised.value.status_code == 400

    def test_route_answers_while_warm(self, fake_hera, warmed):
        fake_hera(["Alpha"], {"Alpha": [{"x": 1}]})

        assert routes.documents("Alpha") == {"project": {"name": "Alpha", "documents": [{"x": 1}]}}

    def test_route_says_warming_up(self, fake_hera, monkeypatch):
        fake_hera(["Alpha"], {"Alpha": []})
        monkeypatch.setattr(dependencies.warmup, "_ready", False)

        with pytest.raises(HTTPException) as raised:
            routes.documents("Alpha")

        assert raised.value.status_code == 503
        assert raised.value.detail == "WARMING_UP"
