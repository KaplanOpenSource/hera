import sys
import types

import pytest

import dependencies
import init_data
from routers import startup


@pytest.fixture
def fake_hera(monkeypatch):
    """Install a fake hera so init_data can be exercised without a real DB.

    Returns a factory: ``fake_hera(projects, documents, broken=())``. Names in
    ``broken`` raise instead of answering, to test partial failure.
    """
    def _install(projects, documents=None, broken=()):
        def get_project_list():
            if "projects" in broken:
                raise RuntimeError("no mongo")
            return projects

        def get_documents_as_dict(projectName, with_id=False):
            if "project" in broken:
                raise RuntimeError("bad document")
            return {"documents": (documents or {}).get(projectName, [])}

        def get_toolkit_documents():
            if "toolkitDocs" in broken:
                raise KeyError("classpath")
            return [{"toolkit": "LSM", "desc": {"classpath": "lsm.cls"}}]

        datatypes = types.SimpleNamespace()
        datatypes.STRING = "string"
        datatypes.JSON_DICT = "JSON_dict"
        datatypes.lowercase = "skipped"
        datatypes.aFunction = lambda: "skipped"

        fake = types.ModuleType("hera")
        fake.__version__ = "9.9.9"
        fake.toolkitHome = types.SimpleNamespace(getToolkitDocuments=get_toolkit_documents)
        fake.datalayer = types.SimpleNamespace(
            datatypes=datatypes,
            All=types.SimpleNamespace(getDocumentsAsDict=get_documents_as_dict),
        )

        project_module = types.ModuleType("hera.datalayer.project")
        project_module.getProjectList = get_project_list
        datalayer_module = types.ModuleType("hera.datalayer")
        datalayer_module.All = fake.datalayer.All
        datalayer_module.datatypes = datatypes

        monkeypatch.setitem(sys.modules, "hera", fake)
        monkeypatch.setitem(sys.modules, "hera.datalayer", datalayer_module)
        monkeypatch.setitem(sys.modules, "hera.datalayer.project", project_module)

    return _install


def test_chooses_the_wanted_project(fake_hera):
    fake_hera(["Alpha", "Beta"], {"Beta": [{"x": 1}]})

    data = init_data.collect_init("Beta")

    assert data["ready"] is True
    assert data["project"] == {"name": "Beta", "documents": [{"x": 1}]}
    assert "problems" not in data


def test_falls_back_to_the_first_project(fake_hera):
    fake_hera(["Alpha", "Beta"], {"Alpha": []})

    data = init_data.collect_init("NoSuchProject")

    assert data["project"]["name"] == "Alpha"


def test_default_project_comes_first(fake_hera):
    fake_hera(["zeta", init_data.DEFAULT_PROJECT, "Alpha"], {})

    data = init_data.collect_init("")

    assert [p["name"] for p in data["projects"]] == [init_data.DEFAULT_PROJECT, "Alpha", "zeta"]
    assert data["project"]["name"] == init_data.DEFAULT_PROJECT


def test_no_projects_means_no_project(fake_hera):
    fake_hera([], {})

    data = init_data.collect_init("")

    assert data["projects"] == []
    assert "project" not in data


def test_carries_datatypes_and_server_info(fake_hera):
    fake_hera(["Alpha"], {})

    data = init_data.collect_init("Alpha")

    assert data["datatypes"] == {"STRING": "string", "JSON_DICT": "JSON_dict"}
    assert data["server"]["heraVersion"] == "9.9.9"
    assert data["server"]["username"]
    assert isinstance(data["server"]["inDocker"], bool)


def test_one_broken_piece_leaves_the_rest(fake_hera):
    fake_hera(["Alpha"], {"Alpha": []}, broken=["toolkitDocs"])

    data = init_data.collect_init("Alpha")

    assert "toolkitDocs" not in data
    assert "KeyError" in data["problems"]["toolkitDocs"]
    assert data["project"]["name"] == "Alpha"
    assert data["datatypes"]


def test_broken_project_list_still_answers(fake_hera):
    fake_hera(["Alpha"], {}, broken=["projects"])

    data = init_data.collect_init("Alpha")

    assert data["ready"] is True
    assert "projects" not in data
    assert "project" not in data
    assert "RuntimeError" in data["problems"]["projects"]


def test_route_says_not_ready_while_warming(monkeypatch):
    monkeypatch.setattr(dependencies.warmup, "_ready", False)

    assert startup.init() == {"ready": False}


def test_route_returns_the_payload_when_ready(monkeypatch, fake_hera):
    fake_hera(["Alpha"], {"Alpha": []})
    monkeypatch.setattr(dependencies.warmup, "_ready", True)

    data = startup.init("Alpha")

    assert data["ready"] is True
    assert data["project"]["name"] == "Alpha"
