"""Every project route refuses while hera is warming up."""
import types

import pytest
from fastapi import HTTPException

import dependencies
from load_helpers import _load, _zip_of
from routers import projects as routes


# Every route reaches hera, so each one refuses while hera is still warming.
WARMING_CALLS = {
    "names": lambda: routes.names(),
    "folders": lambda: routes.folders(),
    "documents": lambda: routes.documents("Alpha"),
    "files": lambda: routes.files("Alpha", ""),
    "toolkits": lambda: routes.toolkits(),
    "create": lambda: routes.create(types.SimpleNamespace(
        projectName="Alpha", filesDirectory="", loadRepositories=True)),
    "delete": lambda: routes.delete(types.SimpleNamespace(
        projectName="Alpha", deleteFiles=False)),
    "load": lambda: _load(_zip_of([{"projectName": "Beta"}])),
    "detectNotebooks": lambda: routes.detect(types.SimpleNamespace(
        projectName="Alpha", filesDirectory="/tmp")),
}


@pytest.mark.parametrize("route", sorted(WARMING_CALLS))
def test_every_route_says_warming_up(route, fake_hera, monkeypatch):
    fake_hera(["Alpha"], {"Alpha": []})
    monkeypatch.setattr(dependencies.warmup, "_ready", False)

    with pytest.raises(HTTPException) as raised:
        WARMING_CALLS[route]()

    assert raised.value.status_code == 503
    assert raised.value.detail == "WARMING_UP"
