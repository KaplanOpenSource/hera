"""Fixtures for the project entrypoint tests: a fake hera and a warmed-up server."""
import sys
import types

import pytest

import dependencies
from fake_database import FakeDatabase


@pytest.fixture
def fake_hera(monkeypatch):
    """Install the fake hera and return the database behind it."""
    def _install(projects, documents=None, filesDirectories=None):
        db = FakeDatabase(projects, documents, filesDirectories)

        datalayer = types.ModuleType("hera.datalayer")
        datalayer.All = types.SimpleNamespace(
            getDocumentsAsDict=db.getDocumentsAsDict,
            getDocuments=db.getDocuments,
            deleteDocumentByID=db.deleteDocumentByID)
        datalayer.Cache_Collection = db.Cache_Collection

        project_module = types.ModuleType("hera.datalayer.project")
        project_module.getProjectList = db.getProjectList
        # A function, not the bound method, so it can carry hera's static methods.
        def project_class(*args, **kwargs):
            return db.Project(*args, **kwargs)
        setattr(project_class, "_iter_pickled_docs", db._iter_pickled_docs)
        setattr(project_class, "load", db.load)
        project_module.Project = project_class

        cli_module = types.ModuleType("hera.utils.data.CLI")
        cli_module.project_create = db.project_create

        hera = types.ModuleType("hera")
        hera.__version__ = "9.9.9"
        hera.toolkitHome = types.SimpleNamespace(
            getToolkitDocuments=lambda: [{"toolkit": "LSM", "desc": {"classpath": "lsm.cls"}}])
        hera.datalayer = datalayer

        monkeypatch.setitem(sys.modules, "hera", hera)
        monkeypatch.setitem(sys.modules, "hera.datalayer", datalayer)
        monkeypatch.setitem(sys.modules, "hera.datalayer.project", project_module)
        monkeypatch.setitem(sys.modules, "hera.utils", types.ModuleType("hera.utils"))
        monkeypatch.setitem(sys.modules, "hera.utils.data", types.ModuleType("hera.utils.data"))
        monkeypatch.setitem(sys.modules, "hera.utils.data.CLI", cli_module)
        return db

    return _install


@pytest.fixture
def warmed(monkeypatch):
    monkeypatch.setattr(dependencies.warmup, "_ready", True)
