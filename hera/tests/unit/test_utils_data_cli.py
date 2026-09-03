"""utils/data/CLI.py: the lazy-import bootstrap plus the argparse command
handlers behind ``hera-project``.

Covered here: ``_setup``, ``_lazy_setup`` (and its ``wrapper``),
``project_list``, ``project_create``, ``project_dump``, ``project_load``,
``repository_list``, ``repository_add``, ``repository_remove``,
``repository_show``, ``repository_load``, ``repository_export``,
``add_toolkit``, ``display_datasource_versions``,
``update_datasource_default_version``, ``update``, ``populate``,
``db_list``, ``db_create``, ``db_remove``, ``toolkit_list``,
``toolkit_register``, ``toolkit_load``, ``toolkit_default_repo_show``,
``toolkit_default_repo_set``, ``toolkit_import_json`` and
``project_measurements_list``.

The two pure helpers of the same module -- ``_parse_query_value`` and
``load_project_name`` -- already have their own file
(``test_utils_data_cli_pure.py``) and are not re-tested here.

Every handler is exercised against fakes injected over the module globals
that ``_setup`` would otherwise fill in (with ``_setup`` itself stubbed
out), so the tests assert on the CLI plumbing -- which collaborator was
built, which arguments were forwarded, what was printed, which branch was
taken, what was raised -- and never touch MongoDB or the data layer.
``TestTheFakesDoNotDrift`` checks each faked method really exists on the
class it stands in for, so a rename in production has to break a test.

Deliberately not tested: the argparse wiring itself (``hera/bin/hera-project``
is an extensionless script whose parser is built inside ``main()`` and
cannot be imported), and the behaviour of the collaborators the handlers
delegate to (``dataToolkit.addRepository``, ``ToolkitHome.registerToolkit``,
``Project.addDocumentFromDict``, ...), which belong to their own units.
The exact pandas/tabulate table rendering is asserted only through the
values that must appear in it, never character by character.

Three bugs are pinned here:

* B179 -- ``toolkit_register`` calls ``registerToolkit`` with the
  undefined names ``toolkit_name``/``toolkit_path``, so the command can
  never register anything.
* B180 -- ``update_datasource_default_version`` calls
  ``setDataSourceDefaultVersion`` on a ``Project``, but that method only
  exists on ``abstractToolkit``.
* B181 -- ``add_toolkit`` reports a non-object ``--params`` JSON as
  invalid JSON, hiding the real reason.
"""
import json
import logging
import os
from argparse import Namespace

import pandas
import pytest
from tabulate import tabulate

from hera.utils.data import CLI


# ---------------------------------------------------------------------------
# fakes
# ---------------------------------------------------------------------------

class _Recorder:
    """Common call-recording machinery for every fake below."""

    def __init__(self):
        self.calls = []

    def _record(self, callName, /, *args, **kwargs):
        # positional-only, so a recorded call may itself take a `callName`
        # or `name` keyword -- toolkit_list passes name=None.
        self.calls.append((callName, args, kwargs))

    def names(self):
        return [entry[0] for entry in self.calls]

    def call(self, name):
        """The single recorded call named ``name``."""
        matches = [entry for entry in self.calls if entry[0] == name]
        assert len(matches) == 1, f"expected exactly one {name} call, got {self.calls}"
        return matches[0]

    def kwargs(self, name):
        return self.call(name)[2]


class _FakeDocument:
    """Stands in for a datalayer document: ``.desc`` plus loose attributes."""

    def __init__(self, desc=None, **attributes):
        self.desc = dict(desc or {})
        self.asDictCalls = []
        for key, value in attributes.items():
            setattr(self, key, value)

    def asDict(self, **kwargs):
        self.asDictCalls.append(kwargs)
        return dict(self.desc, _id="the-id") if kwargs.get("with_id") else dict(self.desc)


class _FakeDataToolkit(_Recorder):
    """Stands in for hera.utils.data.toolkit.dataToolkit."""

    def __init__(self, harness, connectionName=None):
        super().__init__()
        self.harness = harness
        self.connectionName = connectionName

    def getRepositoryTable(self):
        self._record("getRepositoryTable")
        return self.harness.repositoryTable

    def addRepository(self, **kwargs):
        self._record("addRepository", **kwargs)

    def deleteDataSource(self, **kwargs):
        self._record("deleteDataSource", **kwargs)

    def getDataSourceData(self, **kwargs):
        self._record("getDataSourceData", **kwargs)
        return self.harness.repositoryData

    def getDataSourceList(self):
        self._record("getDataSourceList")
        return self.harness.datasourceList

    def loadAllDatasourcesInRepositoryJSONToProject(self, **kwargs):
        self._record("loadAllDatasourcesInRepositoryJSONToProject", **kwargs)

    def loadAllDatasourcesInAllRepositoriesToProject(self, **kwargs):
        self._record("loadAllDatasourcesInAllRepositoriesToProject", **kwargs)
        if kwargs.get("projectName") in self.harness.unpopulatableProjects:
            raise RuntimeError(f"cannot populate {kwargs.get('projectName')}")

    def exportDocumentsToRepository(self, **kwargs):
        self._record("exportDocumentsToRepository", **kwargs)
        return self.harness.exportReport


def _makeToolkitHomeClass(harness):
    """A ToolkitHome stand-in whose instances register with ``harness``."""

    class _FakeToolkitHome(_Recorder):
        def __init__(self, projectName=None, **kwargs):
            super().__init__()
            self.projectName = projectName
            self.buildArguments = kwargs
            harness.homes.append(self)

        def getToolkitDocuments(self, **kwargs):
            self._record("getToolkitDocuments", **kwargs)
            return harness.toolkitDocuments

        def registerToolkit(self, **kwargs):
            self._record("registerToolkit", **kwargs)
            if harness.registerToolkitError is not None:
                raise harness.registerToolkitError

        def getToolkit(self, **kwargs):
            self._record("getToolkit", **kwargs)
            if harness.getToolkitError is not None:
                raise harness.getToolkitError
            return harness.loadedToolkit

        def getDefaultRepository(self, **kwargs):
            self._record("getDefaultRepository", **kwargs)
            return harness.defaultRepository

        def setDefaultRepository(self, **kwargs):
            self._record("setDefaultRepository", **kwargs)

        def import_toolkits_from_json(self, **kwargs):
            self._record("import_toolkits_from_json", **kwargs)
            return harness.importedToolkits

        def import_experiments_from_json(self, **kwargs):
            self._record("import_experiments_from_json", **kwargs)
            return harness.importedExperiments

    return _FakeToolkitHome


def _makeProjectClass(harness):
    """A Project stand-in exposing exactly the methods the handlers call.

    It deliberately does NOT define ``setDataSourceDefaultVersion``,
    because the real ``hera.datalayer.Project`` does not either -- see
    B180 and ``test_the_real_project_class_has_no_...`` below.
    """

    class _FakeProject(_Recorder):
        DEFAULTPROJECT = "defaultProject"

        def __init__(self, projectName=None, **kwargs):
            super().__init__()
            self.projectName = projectName
            self.buildArguments = kwargs
            harness.projects.append(self)

        def getCacheDocuments(self, **kwargs):
            self._record("getCacheDocuments", **kwargs)
            return harness.cacheDocuments

        def getMeasurementsDocuments(self, **kwargs):
            self._record("getMeasurementsDocuments", **kwargs)
            return harness.measurementsDocuments

        def getSimulationsDocuments(self, **kwargs):
            self._record("getSimulationsDocuments", **kwargs)
            return harness.simulationsDocuments

        def getMeasurementsDocumentsAsDict(self, **kwargs):
            self._record("getMeasurementsDocumentsAsDict", **kwargs)
            return dict(documents=harness.documentsAsDict)

        def getConfig(self, **kwargs):
            self._record("getConfig", **kwargs)
            return harness.projectConfig

        def addDocumentFromDict(self, doc):
            self._record("addDocumentFromDict", doc)

    return _FakeProject


class _FakeAll(_Recorder):
    """Stands in for hera.datalayer.All."""

    def __init__(self, harness):
        super().__init__()
        self.harness = harness

    def getDocuments(self, **kwargs):
        self._record("getDocuments", **kwargs)
        return self.harness.allDocuments


class _Harness(_Recorder):
    """Everything ``_setup`` would have bound, plus the canned answers."""

    def __init__(self):
        super().__init__()
        # collaborators built by the handlers
        self.toolkits = []
        self.homes = []
        self.projects = []
        self.All = _FakeAll(self)
        self.Project = _makeProjectClass(self)
        self.ToolkitHome = _makeToolkitHomeClass(self)

        # canned answers
        self.projectList = []
        self.mongoJSON = {}
        self.jsonPayload = None
        self.allDocuments = []
        self.repositoryTable = pandas.DataFrame()
        self.repositoryData = {}
        self.datasourceList = []
        self.exportReport = dict(added=[], skipped_existing=[], overwritten=[])
        self.toolkitDocuments = []
        self.loadedToolkit = None
        self.defaultRepository = None
        self.importedToolkits = None
        self.importedExperiments = None
        self.cacheDocuments = []
        self.measurementsDocuments = []
        self.simulationsDocuments = []
        self.documentsAsDict = []
        self.projectConfig = {}

        # failure injection
        self.unpopulatableProjects = set()
        self.registerToolkitError = None
        self.getToolkitError = None

    # -- the module-level callables --------------------------------------
    def dataToolkit(self, **kwargs):
        self._record("dataToolkit", **kwargs)
        toolkit = _FakeDataToolkit(self, **kwargs)
        self.toolkits.append(toolkit)
        return toolkit

    def getProjectList(self, **kwargs):
        self._record("getProjectList", **kwargs)
        return self.projectList

    def createProjectDirectory(self, **kwargs):
        self._record("createProjectDirectory", **kwargs)

    def removeConnection(self, *args, **kwargs):
        self._record("removeConnection", *args, **kwargs)

    def addOrUpdateDatabase(self, **kwargs):
        self._record("addOrUpdateDatabase", **kwargs)

    def getMongoJSON(self):
        self._record("getMongoJSON")
        return self.mongoJSON

    def loadJSON(self, *args, **kwargs):
        self._record("loadJSON", *args, **kwargs)
        return self.jsonPayload

    # -- accessors -------------------------------------------------------
    @property
    def toolkit(self):
        assert len(self.toolkits) == 1, f"expected one dataToolkit, got {len(self.toolkits)}"
        return self.toolkits[0]

    @property
    def home(self):
        assert len(self.homes) == 1, f"expected one ToolkitHome, got {len(self.homes)}"
        return self.homes[0]

    @property
    def project(self):
        assert len(self.projects) == 1, f"expected one Project, got {len(self.projects)}"
        return self.projects[0]


_GLOBALS = (
    "getProjectList", "createProjectDirectory", "Project", "removeConnection",
    "addOrUpdateDatabase", "getMongoJSON", "datalayer_All", "loadJSON",
    "dataToolkit", "ToolkitHome", "pandas", "tabulate", "toolkitHome",
)


@pytest.fixture()
def wired(monkeypatch):
    """Neutralise the deferred imports and inject the fakes.

    ``_lazy_setup``'s wrapper looks ``_setup`` up in the module globals at
    call time, so replacing it here keeps the real hera imports from
    overwriting the fakes.
    """
    harness = _Harness()
    monkeypatch.setattr(CLI, "_setup", lambda: None)
    monkeypatch.setattr(CLI, "dataToolkit", harness.dataToolkit)
    monkeypatch.setattr(CLI, "ToolkitHome", harness.ToolkitHome)
    monkeypatch.setattr(CLI, "Project", harness.Project)
    monkeypatch.setattr(CLI, "getProjectList", harness.getProjectList)
    monkeypatch.setattr(CLI, "createProjectDirectory", harness.createProjectDirectory)
    monkeypatch.setattr(CLI, "removeConnection", harness.removeConnection)
    monkeypatch.setattr(CLI, "addOrUpdateDatabase", harness.addOrUpdateDatabase)
    monkeypatch.setattr(CLI, "getMongoJSON", harness.getMongoJSON)
    monkeypatch.setattr(CLI, "datalayer_All", harness.All)
    monkeypatch.setattr(CLI, "loadJSON", harness.loadJSON)
    monkeypatch.setattr(CLI, "pandas", pandas)
    monkeypatch.setattr(CLI, "tabulate", tabulate)
    return harness


class _ListHandler(logging.Handler):
    def __init__(self, records):
        super().__init__()
        self.records = records

    def emit(self, record):
        self.records.append(record)


@pytest.fixture()
def captureLog():
    """Capture one named logger's messages.

    ``caplog`` is useless here: hera configures the ``hera`` logger with
    ``propagate`` off, so records never reach the root handler that
    ``caplog`` installs (they surface on stdout instead, which would make
    the assertions depend on hera's logging configuration).  Attaching a
    handler to the named logger directly is both hermetic and precise.
    """
    attached = []

    def _capture(name, level=logging.DEBUG):
        logger = logging.getLogger(name)
        messages = []
        handler = _ListHandler(messages)
        attached.append((logger, handler, logger.level))
        logger.addHandler(handler)
        logger.setLevel(level)
        return messages

    yield _capture

    for logger, handler, previousLevel in attached:
        logger.removeHandler(handler)
        logger.setLevel(previousLevel)


def _messages(records):
    return "\n".join(record.getMessage() for record in records)


@pytest.fixture()
def uninitialized(monkeypatch):
    """Rewind the module to its pre-``_setup`` state, and restore it after."""
    monkeypatch.setattr(CLI, "_initialized", False)
    for name in _GLOBALS:
        monkeypatch.setattr(CLI, name, None)


# ---------------------------------------------------------------------------
# the bootstrap
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestSetup:
    def test_it_binds_every_deferred_name(self, uninitialized):
        CLI._setup()

        assert CLI.pandas is pandas
        assert CLI.tabulate is tabulate
        for name in _GLOBALS:
            assert getattr(CLI, name) is not None, f"{name} was left unbound"
        assert CLI._initialized is True

    def test_it_binds_the_names_to_the_real_hera_objects(self, uninitialized):
        from hera import toolkitHome as realToolkitHome
        from hera.datalayer import All, Project, getProjectList
        from hera.toolkit import ToolkitHome
        from hera.utils.data.toolkit import dataToolkit

        CLI._setup()

        assert CLI.toolkitHome is realToolkitHome
        assert CLI.Project is Project
        assert CLI.getProjectList is getProjectList
        assert CLI.datalayer_All is All
        assert CLI.ToolkitHome is ToolkitHome
        assert CLI.dataToolkit is dataToolkit

    def test_a_second_call_does_not_rebind_anything(self, monkeypatch):
        monkeypatch.setattr(CLI, "_initialized", True)
        monkeypatch.setattr(CLI, "pandas", "sentinel")
        CLI._setup()
        assert CLI.pandas == "sentinel"

    def test_a_fresh_import_of_the_module_binds_nothing(self):
        """The whole point of the bootstrap: importing must stay instant.

        Executed into a throwaway module object rather than through
        ``importlib.reload``, so the module every other test in this file
        shares is left exactly as it was.
        """
        import importlib.util

        spec = importlib.util.spec_from_file_location("_probeDataCLI", CLI.__file__)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        assert module._initialized is False
        for name in _GLOBALS:
            assert getattr(module, name) is None, f"{name} was bound at import time"


@pytest.mark.unit
class TestLazySetup:
    def test_the_wrapper_runs_setup_before_the_handler(self, monkeypatch):
        events = []
        monkeypatch.setattr(CLI, "_setup", lambda: events.append("setup"))

        @CLI._lazy_setup
        def probe(arguments):
            events.append(("probe", arguments))
            return "returned"

        assert probe("ARGS") == "returned"
        assert events == ["setup", ("probe", "ARGS")]

    def test_it_forwards_positional_and_keyword_arguments(self, monkeypatch):
        monkeypatch.setattr(CLI, "_setup", lambda: None)

        @CLI._lazy_setup
        def probe(*args, **kwargs):
            return args, kwargs

        assert probe(1, 2, key="value") == ((1, 2), {"key": "value"})

    def test_the_wrapper_does_not_swallow_handler_exceptions(self, monkeypatch):
        monkeypatch.setattr(CLI, "_setup", lambda: None)

        @CLI._lazy_setup
        def probe(arguments):
            raise KeyError("boom")

        with pytest.raises(KeyError, match="boom"):
            probe(None)

    def test_the_decorated_handlers_keep_their_own_identity(self):
        assert CLI.project_list.__name__ == "project_list"
        assert "List all the projects" in CLI.project_list.__doc__
        assert CLI.project_list.__wrapped__.__name__ == "project_list"

    @pytest.mark.parametrize("name", [
        "project_list", "project_create", "project_dump", "project_load",
        "repository_list", "repository_add", "repository_remove",
        "repository_show", "repository_load", "repository_export",
        "add_toolkit", "display_datasource_versions",
        "update_datasource_default_version", "update", "populate",
        "db_list", "db_create", "db_remove", "toolkit_list",
        "toolkit_register", "toolkit_load", "toolkit_default_repo_show",
        "toolkit_default_repo_set", "toolkit_import_json",
        "project_measurements_list",
    ])
    def test_every_command_handler_is_wrapped_in_the_bootstrap(self, name):
        assert hasattr(getattr(CLI, name), "__wrapped__"), \
            f"{name} is not decorated with @_lazy_setup"


# ---------------------------------------------------------------------------
# project commands
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestProjectList:
    def test_it_prints_one_row_per_project(self, wired, capsys):
        wired.projectList = ["BETA", "ALPHA"]
        CLI.project_list(Namespace(connectionName="conn", fulldetails=False))

        out = capsys.readouterr().out
        assert "Projects in the connection conn" in out
        assert "ALPHA" in out and "BETA" in out

    def test_the_projects_are_sorted_by_name(self, wired, capsys):
        wired.projectList = ["ZULU", "ALPHA", "MIKE"]
        CLI.project_list(Namespace(connectionName="conn", fulldetails=False))

        out = capsys.readouterr().out
        assert out.index("ALPHA") < out.index("MIKE") < out.index("ZULU")

    def test_the_connection_name_is_forwarded_to_getprojectlist(self, wired):
        CLI.project_list(Namespace(connectionName="myConnection", fulldetails=False))
        assert wired.kwargs("getProjectList") == dict(connectionName="myConnection")

    def test_a_missing_connection_name_falls_back_to_the_os_user(self, wired, monkeypatch):
        import getpass

        monkeypatch.setattr(getpass, "getuser", lambda: "theUser")
        CLI.project_list(Namespace(connectionName=None, fulldetails=False))
        assert wired.kwargs("getProjectList") == dict(connectionName="theUser")

    def test_without_full_details_no_project_is_opened(self, wired):
        wired.projectList = ["ALPHA"]
        CLI.project_list(Namespace(connectionName="conn", fulldetails=False))
        assert wired.projects == []

    def test_full_details_adds_the_three_document_counts(self, wired, capsys):
        wired.projectList = ["ALPHA"]
        wired.cacheDocuments = [1, 2]
        wired.measurementsDocuments = [1, 2, 3]
        wired.simulationsDocuments = [1]
        CLI.project_list(Namespace(connectionName="conn", fulldetails=True))

        out = capsys.readouterr().out
        assert "Cache documents" in out
        assert "Measurements documents" in out
        assert "Simulation documents" in out
        assert wired.project.names() == ["getCacheDocuments",
                                         "getMeasurementsDocuments",
                                         "getSimulationsDocuments"]

    def test_full_details_opens_each_project_on_the_same_connection(self, wired):
        wired.projectList = ["ALPHA", "BETA"]
        CLI.project_list(Namespace(connectionName="conn", fulldetails=True))

        assert [project.projectName for project in wired.projects] == ["ALPHA", "BETA"]
        assert all(project.buildArguments == dict(connectionName="conn")
                   for project in wired.projects)

    def test_an_empty_project_list_prints_an_empty_frame_without_sorting(self, wired, capsys):
        CLI.project_list(Namespace(connectionName="conn", fulldetails=True))
        assert "Empty DataFrame" in capsys.readouterr().out


@pytest.mark.unit
class TestProjectCreate:
    def _arguments(self, **overrides):
        arguments = Namespace(projectName="NEW_PROJECT", directory=None,
                              loadRepositories=False, overwrite=False)
        for key, value in overrides.items():
            setattr(arguments, key, value)
        return arguments

    def test_an_explicit_directory_is_used_verbatim(self, wired, capsys, tmp_path):
        CLI.project_create(self._arguments(directory=str(tmp_path)))

        assert wired.kwargs("createProjectDirectory") == dict(
            outputPath=str(tmp_path), projectName="NEW_PROJECT")
        assert f"Created project NEW_PROJECT in directory {tmp_path}" in capsys.readouterr().out

    def test_a_missing_directory_becomes_the_project_name_under_the_cwd(self, wired, monkeypatch, tmp_path):
        monkeypatch.chdir(tmp_path)
        CLI.project_create(self._arguments())

        assert wired.kwargs("createProjectDirectory")["outputPath"] == \
            os.path.join(str(tmp_path), "NEW_PROJECT")

    def test_the_directory_is_not_created_by_the_handler_itself(self, wired, monkeypatch, tmp_path):
        """Creating it is ``createProjectDirectory``'s job, not the CLI's."""
        monkeypatch.chdir(tmp_path)
        CLI.project_create(self._arguments())
        assert not (tmp_path / "NEW_PROJECT").exists()

    def test_no_repositories_are_loaded_when_the_flag_is_off(self, wired, tmp_path):
        CLI.project_create(self._arguments(directory=str(tmp_path), loadRepositories=False))
        assert wired.toolkits == []

    def test_the_repositories_are_loaded_into_the_new_project(self, wired, tmp_path):
        CLI.project_create(self._arguments(directory=str(tmp_path),
                                           loadRepositories=True, overwrite=True))

        assert wired.toolkit.kwargs("loadAllDatasourcesInAllRepositoriesToProject") == \
            dict(projectName="NEW_PROJECT", overwrite=True)

    def test_the_datatoolkit_is_built_without_a_connection_name(self, wired, tmp_path):
        CLI.project_create(self._arguments(directory=str(tmp_path), loadRepositories=True))
        assert wired.kwargs("dataToolkit") == {}


@pytest.mark.unit
class TestProjectDump:
    def _arguments(self, **overrides):
        arguments = Namespace(projectName="DUMP_PROJECT", query=[],
                              fileName=None, outputFormat="json")
        for key, value in overrides.items():
            setattr(arguments, key, value)
        return arguments

    def test_the_project_name_is_always_part_of_the_query(self, wired):
        CLI.project_dump(self._arguments())
        assert wired.All.kwargs("getDocuments") == dict(projectName="DUMP_PROJECT")

    def test_each_query_element_is_split_on_the_first_equals_sign(self, wired):
        CLI.project_dump(self._arguments(query=["type=Cache", "desc=a=b"]))
        assert wired.All.kwargs("getDocuments") == dict(
            projectName="DUMP_PROJECT", type="Cache", desc="a=b")

    def test_query_keys_are_stripped_and_values_are_typed(self, wired):
        CLI.project_dump(self._arguments(query=["  height = 10 ", "flag=true",
                                                "ratio=0.5", "missing=none"]))
        forwarded = wired.All.kwargs("getDocuments")
        assert forwarded["height"] == 10
        assert forwarded["flag"] is True
        assert forwarded["ratio"] == pytest.approx(0.5)
        assert forwarded["missing"] is None

    def test_a_query_element_without_a_value_becomes_an_empty_string(self, wired):
        """``partition`` yields '' for a bare key, which is forwarded as-is."""
        CLI.project_dump(self._arguments(query=["type"]))
        assert wired.All.kwargs("getDocuments")["type"] == ""

    def test_json_output_is_the_indented_document_list(self, wired, capsys):
        wired.allDocuments = [_FakeDocument(dict(type="Cache"))]
        CLI.project_dump(self._arguments(outputFormat="json"))

        printed = json.loads(capsys.readouterr().out)
        assert printed == [dict(type="Cache", _id="the-id")]

    def test_every_document_is_dumped_with_its_id(self, wired):
        documents = [_FakeDocument(dict(type="Cache")), _FakeDocument(dict(type="Sim"))]
        wired.allDocuments = documents
        CLI.project_dump(self._arguments())

        assert [document.asDictCalls for document in documents] == \
            [[dict(with_id=True)], [dict(with_id=True)]]

    def test_table_output_prints_a_dataframe_instead_of_json(self, wired, capsys):
        wired.allDocuments = [_FakeDocument(dict(type="Cache", resource="/a"))]
        CLI.project_dump(self._arguments(outputFormat="table"))

        out = capsys.readouterr().out
        assert "type" in out and "Cache" in out
        assert not out.lstrip().startswith("[")

    def test_an_unknown_output_format_prints_a_hint_and_nothing_else(self, wired, capsys):
        wired.allDocuments = [_FakeDocument(dict(type="Cache"))]
        CLI.project_dump(self._arguments(outputFormat="yaml"))

        out = capsys.readouterr().out
        assert out.strip() == "The outputFormat must be 'json' or 'table'"

    def test_the_documents_are_also_written_to_the_requested_file(self, wired, tmp_path):
        wired.allDocuments = [_FakeDocument(dict(type="Cache"))]
        target = tmp_path / "dump.json"
        CLI.project_dump(self._arguments(fileName=str(target)))

        assert json.loads(target.read_text()) == [dict(type="Cache", _id="the-id")]

    def test_an_empty_project_writes_an_empty_json_list(self, wired, tmp_path):
        target = tmp_path / "dump.json"
        CLI.project_dump(self._arguments(fileName=str(target)))
        assert json.loads(target.read_text()) == []

    def test_a_bare_filename_is_written_relative_to_the_cwd(self, wired, monkeypatch, tmp_path):
        """There is no --directory for this command, so the CWD is the only
        sensible base; asserted so a later change cannot pass unnoticed."""
        monkeypatch.chdir(tmp_path)
        CLI.project_dump(self._arguments(fileName="dump.json"))
        assert (tmp_path / "dump.json").exists()

    def test_writing_to_a_file_still_prints_the_json(self, wired, capsys, tmp_path):
        CLI.project_dump(self._arguments(fileName=str(tmp_path / "dump.json")))
        assert capsys.readouterr().out.strip() == "[]"


@pytest.mark.unit
class TestProjectLoad:
    def _arguments(self, path="dump.json", **overrides):
        arguments = Namespace(projectName="LOAD_PROJECT", file=path)
        for key, value in overrides.items():
            setattr(arguments, key, value)
        return arguments

    def test_the_file_name_is_handed_to_loadjson(self, wired):
        wired.jsonPayload = []
        CLI.project_load(self._arguments("/some/dump.json"))
        assert wired.call("loadJSON")[1] == ("/some/dump.json",)

    def test_every_document_is_added_to_the_project(self, wired):
        wired.jsonPayload = [dict(a=1), dict(b=2)]
        CLI.project_load(self._arguments())

        added = [entry[1][0] for entry in wired.project.calls
                 if entry[0] == "addDocumentFromDict"]
        assert added == [dict(a=1), dict(b=2)]
        assert wired.project.projectName == "LOAD_PROJECT"

    def test_progress_is_printed_one_line_per_document(self, wired, capsys):
        wired.jsonPayload = [dict(a=1), dict(b=2)]
        CLI.project_load(self._arguments())

        assert capsys.readouterr().out.splitlines() == [
            "Loading document 1/2", "Loading document 2/2"]

    def test_a_json_object_instead_of_a_list_is_refused(self, wired, capsys):
        wired.jsonPayload = dict(a=1)
        CLI.project_load(self._arguments())

        assert "ERROR: Expected a JSON list of documents, got dict" in capsys.readouterr().out
        assert wired.project.names() == []

    def test_a_json_scalar_instead_of_a_list_is_refused(self, wired, capsys):
        wired.jsonPayload = "not a list"
        CLI.project_load(self._arguments())
        assert "got str" in capsys.readouterr().out

    def test_an_empty_list_loads_nothing_and_prints_nothing(self, wired, capsys):
        wired.jsonPayload = []
        CLI.project_load(self._arguments())

        assert capsys.readouterr().out == ""
        assert wired.project.names() == []

    def test_the_project_is_opened_before_the_payload_is_validated(self, wired):
        """Characterisation: the refusal path has already built a Project."""
        wired.jsonPayload = dict(a=1)
        CLI.project_load(self._arguments())
        assert len(wired.projects) == 1


# ---------------------------------------------------------------------------
# repository commands
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestRepositoryList:
    def test_an_empty_table_prints_a_friendly_message(self, wired, capsys):
        wired.repositoryTable = pandas.DataFrame()
        CLI.repository_list(Namespace())

        assert capsys.readouterr().out.strip() == "The user does not have repositories."

    def test_listing_prints_one_row_per_repository(self, wired, capsys):
        wired.repositoryTable = pandas.DataFrame([
            dict(repository="repoA", path="/a/repoA.json"),
            dict(repository="repoB", path="/b/repoB.json"),
        ])
        CLI.repository_list(Namespace())

        out = capsys.readouterr().out
        assert "repoA" in out and "repoB" in out
        assert "/a/repoA.json" in out

    def test_it_asks_the_data_toolkit_for_the_table(self, wired):
        CLI.repository_list(Namespace())
        assert wired.toolkit.names() == ["getRepositoryTable"]

    def test_the_command_ignores_its_arguments(self, wired, capsys):
        """The handler's parameter is spelled ``argumets`` and never read."""
        CLI.repository_list(Namespace(anything="at all"))
        assert capsys.readouterr().out.strip() == "The user does not have repositories."


@pytest.mark.unit
class TestRepositoryAdd:
    def test_the_repository_name_is_derived_from_the_file_name(self, wired):
        CLI.repository_add(Namespace(repositoryName="/data/repos/myRepo.json",
                                     overwrite=False))

        assert wired.toolkit.kwargs("addRepository") == dict(
            repositoryName="myRepo",
            repositoryPath="/data/repos/myRepo.json",
            overwrite=False)

    def test_the_overwrite_flag_is_forwarded(self, wired):
        CLI.repository_add(Namespace(repositoryName="myRepo.json", overwrite=True))
        assert wired.toolkit.kwargs("addRepository")["overwrite"] is True

    def test_a_name_without_a_suffix_is_used_as_is(self, wired):
        CLI.repository_add(Namespace(repositoryName="/data/plainRepo", overwrite=False))
        assert wired.toolkit.kwargs("addRepository")["repositoryName"] == "plainRepo"

    def test_a_dotted_file_name_is_truncated_at_the_first_dot(self, wired):
        """``basename(...).split('.')[0]`` is not ``splitext``: a repository
        file called ``meteo.v2.json`` registers under the name ``meteo``,
        so ``meteo.v3.json`` would collide with it."""
        CLI.repository_add(Namespace(repositoryName="/data/meteo.v2.json", overwrite=False))
        assert wired.toolkit.kwargs("addRepository")["repositoryName"] == "meteo"

    def test_the_full_path_is_kept_as_the_repository_path(self, wired):
        CLI.repository_add(Namespace(repositoryName="../rel/repo.json", overwrite=False))
        assert wired.toolkit.kwargs("addRepository")["repositoryPath"] == "../rel/repo.json"

    def test_it_logs_the_registration(self, wired, captureLog):
        records = captureLog("hera.bin.repository_add")
        CLI.repository_add(Namespace(repositoryName="/d/myRepo.json", overwrite=False))

        assert "Adding the repository /d/myRepo.json as name myRepo" in _messages(records)


@pytest.mark.unit
class TestRepositoryRemove:
    def test_the_repository_name_is_removed_as_a_datasource(self, wired):
        CLI.repository_remove(Namespace(repositoryName="myRepo"))
        assert wired.toolkit.kwargs("deleteDataSource") == dict(datasourceName="myRepo")

    def test_it_does_not_ask_for_confirmation(self, wired, capsys):
        CLI.repository_remove(Namespace(repositoryName="myRepo"))
        assert capsys.readouterr().out == ""


@pytest.mark.unit
class TestRepositoryShow:
    def _repository(self, **overrides):
        repository = {
            "MeteoLowFreq": {
                "DataSource": {
                    "YAVNEEL": dict(isRelativePath="True",
                                    item=dict(resource="a/b.parquet",
                                              dataFormat="parquet")),
                },
            },
        }
        repository.update(overrides)
        return repository

    def test_it_asks_for_the_repository_by_name(self, wired):
        wired.repositoryData = self._repository()
        CLI.repository_show(Namespace(repositoryName="meteo_data_v1"))
        assert wired.toolkit.kwargs("getDataSourceData") == dict(datasourceName="meteo_data_v1")

    def test_the_toolkit_name_heads_the_output(self, wired, capsys):
        wired.repositoryData = self._repository()
        CLI.repository_show(Namespace(repositoryName="meteo_data_v1"))

        out = capsys.readouterr().out
        assert "Toolkit:" in out and "MeteoLowFreq" in out

    def test_all_four_datatype_sections_are_always_printed(self, wired, capsys):
        wired.repositoryData = self._repository()
        CLI.repository_show(Namespace(repositoryName="meteo_data_v1"))

        out = capsys.readouterr().out
        for datatype in ("DataSource", "Measurements", "Cache", "Simulations"):
            assert datatype in out

    def test_each_datasource_shows_its_name_path_flag_and_items(self, wired, capsys):
        wired.repositoryData = self._repository()
        CLI.repository_show(Namespace(repositoryName="meteo_data_v1"))

        out = capsys.readouterr().out
        assert "YAVNEEL" in out
        assert "Is it relative path? True" in out
        assert "a/b.parquet" in out and "parquet" in out

    def test_a_datatype_outside_the_known_four_is_silently_ignored(self, wired, capsys):
        wired.repositoryData = {"MeteoLowFreq": {"Config": dict(stationType="IMS")}}
        CLI.repository_show(Namespace(repositoryName="meteo_data_v1"))
        assert "IMS" not in capsys.readouterr().out

    def test_an_empty_repository_prints_nothing_at_all(self, wired, capsys):
        wired.repositoryData = {}
        CLI.repository_show(Namespace(repositoryName="empty"))
        assert capsys.readouterr().out == ""

    def test_a_datasource_without_isrelativepath_raises_keyerror(self, wired):
        wired.repositoryData = {"T": {"DataSource": {"D": dict(item={})}}}
        with pytest.raises(KeyError, match="isRelativePath"):
            CLI.repository_show(Namespace(repositoryName="broken"))


@pytest.mark.unit
class TestRepositoryLoad:
    def test_the_repository_json_is_read_and_forwarded(self, wired, tmp_path):
        wired.jsonPayload = {"MeteoLowFreq": {}}
        path = str(tmp_path / "repo.json")
        CLI.repository_load(Namespace(repositoryName=path, projectName="P", overwrite=False))

        assert wired.call("loadJSON")[1] == (path,)
        forwarded = wired.toolkit.kwargs("loadAllDatasourcesInRepositoryJSONToProject")
        assert forwarded["repositoryJSON"] == {"MeteoLowFreq": {}}
        assert forwarded["projectName"] == "P"
        assert forwarded["overwrite"] is False

    def test_the_basedir_is_the_repository_files_own_directory(self, wired, tmp_path):
        path = str(tmp_path / "sub" / "repo.json")
        CLI.repository_load(Namespace(repositoryName=path, projectName="P", overwrite=False))

        forwarded = wired.toolkit.kwargs("loadAllDatasourcesInRepositoryJSONToProject")
        assert forwarded["basedir"] == str(tmp_path / "sub")

    def test_a_relative_repository_path_gets_an_absolute_basedir(self, wired, monkeypatch, tmp_path):
        monkeypatch.chdir(tmp_path)
        CLI.repository_load(Namespace(repositoryName="repo.json", projectName="P", overwrite=False))

        forwarded = wired.toolkit.kwargs("loadAllDatasourcesInRepositoryJSONToProject")
        assert forwarded["basedir"] == str(tmp_path)

    def test_without_a_projectname_attribute_the_default_project_is_used(self, wired, tmp_path):
        CLI.repository_load(Namespace(repositoryName=str(tmp_path / "repo.json"),
                                      overwrite=False))
        forwarded = wired.toolkit.kwargs("loadAllDatasourcesInRepositoryJSONToProject")
        assert forwarded["projectName"] is None

    def test_the_overwrite_flag_is_forwarded(self, wired, tmp_path):
        CLI.repository_load(Namespace(repositoryName=str(tmp_path / "repo.json"),
                                      projectName="P", overwrite=True))
        forwarded = wired.toolkit.kwargs("loadAllDatasourcesInRepositoryJSONToProject")
        assert forwarded["overwrite"] is True


@pytest.mark.unit
class TestRepositoryExport:
    def _arguments(self, **overrides):
        arguments = Namespace(repositoryName="out.json", toolkitName="MeteoLowFreq",
                              projectName="SRC", documentId=None,
                              idStrategy="contentHash", mode="add",
                              no_register=False, overwrite=False)
        for key, value in overrides.items():
            setattr(arguments, key, value)
        return arguments

    def test_every_argument_reaches_the_export_call(self, wired):
        CLI.repository_export(self._arguments(documentId=["id1", "id2"],
                                              idStrategy="objectId", mode="override",
                                              overwrite=True))

        assert wired.toolkit.kwargs("exportDocumentsToRepository") == dict(
            toolkitName="MeteoLowFreq", repositoryName="out.json",
            projectName="SRC", documents=["id1", "id2"],
            idStrategy="objectId", mode="override",
            register=True, overwrite=True)

    def test_no_register_is_forwarded_inverted_as_register(self, wired):
        CLI.repository_export(self._arguments(no_register=True))
        assert wired.toolkit.kwargs("exportDocumentsToRepository")["register"] is False

    def test_an_empty_document_id_list_exports_everything(self, wired):
        """An empty --documentId list is falsy, so it means 'all documents'."""
        CLI.repository_export(self._arguments(documentId=[]))
        assert wired.toolkit.kwargs("exportDocumentsToRepository")["documents"] is None

    def test_a_missing_documentid_attribute_exports_everything(self, wired):
        arguments = self._arguments()
        del arguments.documentId
        CLI.repository_export(arguments)
        assert wired.toolkit.kwargs("exportDocumentsToRepository")["documents"] is None

    def test_a_missing_projectname_attribute_becomes_none(self, wired):
        arguments = self._arguments()
        del arguments.projectName
        CLI.repository_export(arguments)
        assert wired.toolkit.kwargs("exportDocumentsToRepository")["projectName"] is None

    def test_the_three_counters_are_summarised(self, wired, capsys):
        wired.exportReport = dict(added=["a", "b"], skipped_existing=["c"],
                                  overwritten=["d", "e", "f"])
        CLI.repository_export(self._arguments())

        assert "Export complete: 2 added, 1 skipped, 3 overwritten." in capsys.readouterr().out

    def test_the_deduplication_line_appears_only_in_override_mode_reports(self, wired, capsys):
        wired.exportReport = dict(added=[], skipped_existing=[], overwritten=[],
                                  deduplicated=["x", "y"])
        CLI.repository_export(self._arguments(mode="override"))

        assert "Deduplicated: 2 duplicate entries removed." in capsys.readouterr().out

    def test_a_report_without_deduplication_prints_only_the_summary(self, wired, capsys):
        CLI.repository_export(self._arguments())
        assert "Deduplicated" not in capsys.readouterr().out


# ---------------------------------------------------------------------------
# addToolkit
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestAddToolkit:
    def _arguments(self, **overrides):
        arguments = Namespace(toolkit_name="MyToolkit",
                              toolkit_path="mypkg.mymod.MyToolkit",
                              params=None, version="0,0,1", overwrite=False)
        for key, value in overrides.items():
            setattr(arguments, key, value)
        return arguments

    def test_every_argument_reaches_registertoolkit(self, wired):
        CLI.add_toolkit(self._arguments(overwrite=True))

        assert wired.home.kwargs("registerToolkit") == dict(
            toolkit_name="MyToolkit", toolkit_path="mypkg.mymod.MyToolkit",
            params=None, version=(0, 0, 1), overwrite=True)

    def test_the_toolkithome_is_built_without_a_project(self, wired):
        CLI.add_toolkit(self._arguments())
        assert wired.home.projectName is None

    def test_the_version_string_becomes_a_three_int_tuple(self, wired):
        CLI.add_toolkit(self._arguments(version="1, 2 ,3"))
        assert wired.home.kwargs("registerToolkit")["version"] == (1, 2, 3)

    @pytest.mark.parametrize("version", ["0,0", "0,0,1,2", "a,b,c", "", "0.0.1"])
    def test_a_malformed_version_is_rejected(self, wired, version):
        with pytest.raises(ValueError, match="Invalid version format"):
            CLI.add_toolkit(self._arguments(version=version))
        assert wired.homes == []

    def test_a_json_object_params_string_is_parsed(self, wired):
        CLI.add_toolkit(self._arguments(params='{"a": 1, "b": "two"}'))
        assert wired.home.kwargs("registerToolkit")["params"] == dict(a=1, b="two")

    def test_an_empty_params_string_is_treated_as_absent(self, wired):
        CLI.add_toolkit(self._arguments(params=""))
        assert wired.home.kwargs("registerToolkit")["params"] is None

    def test_unparseable_params_are_rejected(self, wired):
        with pytest.raises(ValueError, match="Invalid JSON passed to --params"):
            CLI.add_toolkit(self._arguments(params="{not json"))
        assert wired.homes == []

    def test_a_registration_failure_is_logged_and_re_raised(self, wired, captureLog):
        records = captureLog("hera.bin.hera_projects")
        wired.registerToolkitError = RuntimeError("already registered")

        with pytest.raises(RuntimeError, match="already registered"):
            CLI.add_toolkit(self._arguments())

        assert "Failed to register toolkit 'MyToolkit': already registered" in _messages(records)

    def test_a_successful_registration_is_logged(self, wired, captureLog):
        records = captureLog("hera.bin.hera_projects")
        CLI.add_toolkit(self._arguments())

        assert "registered successfully" in _messages(records)

    @pytest.mark.xfail(
        strict=True,
        reason="B181: add_toolkit validates --params with `raise "
               "ValueError('Params must be a JSON object.')` inside the same "
               "try block that catches JSON decoding errors, and the bare "
               "`except Exception` re-raises everything as 'Invalid JSON "
               "passed to --params'. A perfectly valid JSON array or scalar "
               "is therefore reported as malformed JSON, so the user is told "
               "to fix their quoting instead of to pass an object. See the "
               "consolidated findings issue.",
    )
    def test_valid_json_that_is_not_an_object_should_say_so(self, wired):
        with pytest.raises(ValueError, match="must be a JSON object"):
            CLI.add_toolkit(self._arguments(params="[1, 2]"))

    def test_valid_json_that_is_not_an_object_is_reported_as_invalid_json(self, wired):
        """Characterisation of B181."""
        with pytest.raises(ValueError) as caught:
            CLI.add_toolkit(self._arguments(params="[1, 2]"))

        assert str(caught.value) == "Invalid JSON passed to --params: [1, 2]"
        assert isinstance(caught.value.__cause__, ValueError)
        assert str(caught.value.__cause__) == "Params must be a JSON object."


# ---------------------------------------------------------------------------
# datasource versions
# ---------------------------------------------------------------------------

def _versionDocument(toolkit="MeteoLowFreq", datasourceName="YAVNEEL", version=(0, 0, 1)):
    return dict(desc=dict(toolkit=toolkit, datasourceName=datasourceName, version=version))


@pytest.mark.unit
class TestDisplayDatasourceVersions:
    def _arguments(self, **overrides):
        arguments = Namespace(projectName="VERSIONS", datasource=None, default=False)
        for key, value in overrides.items():
            setattr(arguments, key, value)
        return arguments

    def test_every_datasource_version_is_tabulated(self, wired, capsys):
        wired.documentsAsDict = [_versionDocument(),
                                 _versionDocument(datasourceName="TEL_AVIV",
                                                  version=(0, 1, 0))]
        CLI.display_datasource_versions(self._arguments())

        out = capsys.readouterr().out
        assert "toolkit" in out and "datasourceName" in out and "version" in out
        assert "YAVNEEL" in out and "TEL_AVIV" in out

    def test_the_project_is_opened_by_name(self, wired):
        CLI.display_datasource_versions(self._arguments())
        assert wired.project.projectName == "VERSIONS"

    def test_a_datasource_filter_keeps_only_the_matching_rows(self, wired, capsys):
        wired.documentsAsDict = [_versionDocument(),
                                 _versionDocument(datasourceName="TEL_AVIV")]
        CLI.display_datasource_versions(self._arguments(datasource="YAVNEEL"))

        out = capsys.readouterr().out
        assert "YAVNEEL" in out and "TEL_AVIV" not in out

    def test_documents_without_the_expected_desc_keys_are_skipped(self, wired, capsys):
        wired.documentsAsDict = [dict(desc=dict(other="thing")), _versionDocument()]
        CLI.display_datasource_versions(self._arguments())

        out = capsys.readouterr().out
        assert "YAVNEEL" in out
        assert "thing" not in out

    def test_an_empty_project_says_so(self, wired, capsys):
        CLI.display_datasource_versions(self._arguments())
        assert "Are you sure project VERSIONS exists?" in capsys.readouterr().out

    def test_an_unmatched_datasource_filter_names_it_in_the_message(self, wired, capsys):
        wired.documentsAsDict = [_versionDocument()]
        CLI.display_datasource_versions(self._arguments(datasource="NOPE"))

        out = capsys.readouterr().out
        assert "Are you sure datasource NOPE and project VERSIONS exists?" in out

    def test_the_default_flag_reads_the_versions_from_the_project_config(self, wired, capsys):
        wired.documentsAsDict = [_versionDocument()]
        wired.projectConfig = {"YAVNEEL_defaultVersion": (0, 0, 2)}
        CLI.display_datasource_versions(self._arguments(default=True))

        out = capsys.readouterr().out
        assert "DEFAULT_VERSION" in out
        assert "(0, 0, 2)" in out
        assert wired.project.names().count("getConfig") == 1

    def test_the_default_flag_drops_datasources_without_a_default(self, wired, capsys):
        wired.documentsAsDict = [_versionDocument()]
        wired.projectConfig = {}
        CLI.display_datasource_versions(self._arguments(default=True))

        assert "No data to display" in capsys.readouterr().out

    def test_the_default_listing_shows_no_document_version_column(self, wired, capsys):
        wired.documentsAsDict = [_versionDocument()]
        wired.projectConfig = {"YAVNEEL_defaultVersion": (0, 0, 2)}
        CLI.display_datasource_versions(self._arguments(default=True))

        columns = [cell.strip() for cell in
                   capsys.readouterr().out.splitlines()[1].strip("|").split("|")]
        assert columns == ["toolkit", "datasourceName", "DEFAULT_VERSION"]

    def test_a_datasource_filter_combines_with_the_default_flag(self, wired, capsys):
        wired.documentsAsDict = [_versionDocument(),
                                 _versionDocument(datasourceName="TEL_AVIV")]
        wired.projectConfig = {"YAVNEEL_defaultVersion": (0, 0, 2),
                               "TEL_AVIV_defaultVersion": (1, 0, 0)}
        CLI.display_datasource_versions(self._arguments(default=True, datasource="YAVNEEL"))

        out = capsys.readouterr().out
        assert "YAVNEEL" in out and "TEL_AVIV" not in out


@pytest.mark.unit
class TestUpdateDatasourceDefaultVersion:
    def _arguments(self, **overrides):
        arguments = Namespace(projectName="VERSIONS", datasource="YAVNEEL",
                              version="0,0,2")
        for key, value in overrides.items():
            setattr(arguments, key, value)
        return arguments

    def test_the_version_string_is_parsed_into_a_tuple_of_ints(self, wired):
        arguments = self._arguments(version=" 1 , 2 , 3 ")
        with pytest.raises(AttributeError):  # B180, pinned below
            CLI.update_datasource_default_version(arguments)
        assert arguments.version == (1, 2, 3)

    def test_a_malformed_version_is_rejected_before_the_project_is_opened(self, wired):
        with pytest.raises(ValueError):
            CLI.update_datasource_default_version(self._arguments(version="zero,one"))
        assert wired.projects == []

    def test_the_project_is_opened_by_name(self, wired):
        with pytest.raises(AttributeError):  # B180
            CLI.update_datasource_default_version(self._arguments())
        assert wired.project.projectName == "VERSIONS"

    @pytest.mark.xfail(
        strict=True,
        reason="B180: update_datasource_default_version builds a "
               "hera.datalayer.Project and calls "
               "proj.setDataSourceDefaultVersion(...), but that method is "
               "defined on hera.toolkit.abstractToolkit, not on Project -- "
               "Project has no __getattr__ fallback either. So "
               "`hera-project project version update <proj> <ds> <d,d,d>` "
               "always dies with AttributeError and no default version can "
               "ever be set from the CLI. See the consolidated findings issue.",
    )
    def test_the_default_version_should_reach_the_project(self, wired):
        CLI.update_datasource_default_version(self._arguments())
        assert wired.project.kwargs("setDataSourceDefaultVersion") == dict(
            datasourceName="YAVNEEL", version=(0, 0, 2))

    def test_it_currently_raises_attributeerror(self, wired):
        """Characterisation of B180."""
        with pytest.raises(AttributeError, match="setDataSourceDefaultVersion"):
            CLI.update_datasource_default_version(self._arguments())

    def test_the_real_project_class_has_no_setdatasourcedefaultversion(self):
        """Characterisation of B180: where the method actually lives."""
        from hera.datalayer import Project
        from hera.toolkit import abstractToolkit

        assert not hasattr(Project, "setDataSourceDefaultVersion")
        assert hasattr(abstractToolkit, "setDataSourceDefaultVersion")


# ---------------------------------------------------------------------------
# updateRepositories / populate
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestUpdate:
    def test_an_explicit_project_name_is_loaded_directly(self, wired):
        CLI.update(Namespace(projectName="MY_PROJECT", overwrite=False))

        assert wired.toolkit.kwargs("loadAllDatasourcesInAllRepositoriesToProject") == \
            dict(projectName="MY_PROJECT", overwrite=False)
        assert wired.names() == ["dataToolkit"]

    def test_the_overwrite_flag_is_forwarded(self, wired):
        CLI.update(Namespace(projectName="MY_PROJECT", overwrite=True))
        assert wired.toolkit.kwargs("loadAllDatasourcesInAllRepositoriesToProject")["overwrite"] is True

    def test_without_a_project_name_the_case_configuration_supplies_it(self, wired, monkeypatch, tmp_path):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "caseConfiguration.json").write_text(json.dumps(dict(projectName="FROM_FILE")))
        wired.jsonPayload = dict(projectName="FROM_FILE")

        arguments = Namespace(projectName=None, overwrite=False)
        CLI.update(arguments)

        assert arguments.projectName == "FROM_FILE"
        assert wired.toolkit.kwargs("loadAllDatasourcesInAllRepositoriesToProject")["projectName"] == "FROM_FILE"
        assert wired.call("loadJSON")[1] == (os.path.join(str(tmp_path), "caseConfiguration.json"),)

    def test_without_a_project_name_and_without_the_file_it_refuses(self, wired, monkeypatch, tmp_path):
        monkeypatch.chdir(tmp_path)
        with pytest.raises(ValueError, match="caseConfiguration json file must be in folder"):
            CLI.update(Namespace(projectName=None, overwrite=False))
        assert wired.toolkits == []

    def test_a_case_configuration_without_a_project_name_is_refused(self, wired, monkeypatch, tmp_path):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "caseConfiguration.json").write_text("{}")
        wired.jsonPayload = {}

        with pytest.raises(ValueError, match="the key 'projectName' does not exist"):
            CLI.update(Namespace(projectName=None, overwrite=False))

    def test_an_empty_project_name_is_treated_as_missing(self, wired, monkeypatch, tmp_path):
        monkeypatch.chdir(tmp_path)
        with pytest.raises(ValueError, match="caseConfiguration"):
            CLI.update(Namespace(projectName="", overwrite=False))

    def test_a_missing_projectname_attribute_is_not_tolerated(self, wired):
        """Characterisation: the later ``'projectName' in arguments`` guard is
        dead code, because line 497 already dereferences the attribute."""
        with pytest.raises(AttributeError, match="projectName"):
            CLI.update(Namespace(overwrite=False))


@pytest.mark.unit
class TestPopulate:
    def _arguments(self, **overrides):
        arguments = Namespace(projectName=None, overwrite=False)
        for key, value in overrides.items():
            setattr(arguments, key, value)
        return arguments

    def test_without_any_repository_it_tells_the_user_what_to_run(self, wired, capsys):
        wired.datasourceList = []
        CLI.populate(self._arguments())

        out = capsys.readouterr().out
        assert "No repositories registered. Use 'hera-project repository add' first." in out
        assert "loadAllDatasourcesInAllRepositoriesToProject" not in wired.toolkit.names()

    def test_a_named_project_is_the_only_one_populated(self, wired, capsys):
        wired.datasourceList = ["repoA"]
        wired.projectList = ["OTHER"]
        CLI.populate(self._arguments(projectName="ONLY_THIS"))

        loaded = [entry[2]["projectName"] for entry in wired.toolkit.calls
                  if entry[0] == "loadAllDatasourcesInAllRepositoriesToProject"]
        assert loaded == ["ONLY_THIS"]
        assert "getProjectList" not in wired.names()

    def test_without_a_named_project_every_project_is_populated(self, wired):
        wired.datasourceList = ["repoA"]
        wired.projectList = ["ALPHA", "BETA"]
        CLI.populate(self._arguments())

        loaded = [entry[2]["projectName"] for entry in wired.toolkit.calls
                  if entry[0] == "loadAllDatasourcesInAllRepositoriesToProject"]
        assert loaded == ["ALPHA", "BETA"]

    def test_the_default_project_is_never_populated(self, wired, capsys):
        wired.datasourceList = ["repoA"]
        wired.projectList = ["ALPHA", wired.Project.DEFAULTPROJECT]
        CLI.populate(self._arguments())

        loaded = [entry[2]["projectName"] for entry in wired.toolkit.calls
                  if entry[0] == "loadAllDatasourcesInAllRepositoriesToProject"]
        assert loaded == ["ALPHA"]
        assert "Success: 1, Failed: 0" in capsys.readouterr().out

    def test_the_banner_lists_the_repositories_projects_and_overwrite(self, wired, capsys):
        wired.datasourceList = ["repoA", "repoB"]
        wired.projectList = ["ALPHA", "BETA"]
        CLI.populate(self._arguments(overwrite=True))

        out = capsys.readouterr().out
        assert "Repositories: repoA, repoB" in out
        assert "Projects: 2" in out
        assert "Overwrite: True" in out

    def test_the_overwrite_flag_is_forwarded_to_every_project(self, wired):
        wired.datasourceList = ["repoA"]
        wired.projectList = ["ALPHA"]
        CLI.populate(self._arguments(overwrite=True))
        assert wired.toolkit.kwargs("loadAllDatasourcesInAllRepositoriesToProject")["overwrite"] is True

    def test_a_failing_project_is_counted_and_does_not_stop_the_run(self, wired, capsys, captureLog):
        records = captureLog("hera.bin.populate")
        wired.datasourceList = ["repoA"]
        wired.projectList = ["ALPHA", "BROKEN", "GAMMA"]
        wired.unpopulatableProjects = {"BROKEN"}

        CLI.populate(self._arguments())

        assert "Done. Success: 2, Failed: 1" in capsys.readouterr().out
        assert "Failed to populate BROKEN" in _messages(records)

    def test_each_project_is_announced_before_it_is_populated(self, wired, capsys):
        wired.datasourceList = ["repoA"]
        wired.projectList = ["ALPHA"]
        CLI.populate(self._arguments())
        assert "Populating ALPHA..." in capsys.readouterr().out

    def test_missing_overwrite_and_projectname_attributes_are_tolerated(self, wired, capsys):
        wired.datasourceList = ["repoA"]
        wired.projectList = ["ALPHA"]
        CLI.populate(Namespace())

        assert wired.toolkit.kwargs("loadAllDatasourcesInAllRepositoriesToProject") == \
            dict(projectName="ALPHA", overwrite=False)
        assert "Overwrite: False" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# db commands
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestDbList:
    def test_it_prints_one_row_per_connection(self, wired, capsys):
        wired.mongoJSON = {"alice": dict(dbIP="1.2.3.4", dbName="alice_db"),
                           "bob": dict(dbIP="5.6.7.8", dbName="bob_db")}
        CLI.db_list(Namespace(fulldetails=False))

        out = capsys.readouterr().out
        assert "alice" in out and "bob" in out

    def test_full_details_renames_dbip_and_dbname_for_display(self, wired, capsys):
        wired.mongoJSON = {"alice": dict(dbIP="1.2.3.4", dbName="alice_db",
                                         username="alice", password="secret")}
        CLI.db_list(Namespace(fulldetails=True))

        out = capsys.readouterr().out
        assert "IP" in out and "databaseName" in out
        assert "dbIP" not in out and "dbName" not in out
        assert "1.2.3.4" in out and "alice_db" in out

    def test_without_full_details_only_the_connection_name_is_shown(self, wired, capsys):
        wired.mongoJSON = {"alice": dict(dbIP="1.2.3.4", dbName="alice_db")}
        CLI.db_list(Namespace(fulldetails=False))

        out = capsys.readouterr().out
        assert "alice" in out
        assert "1.2.3.4" not in out

    def test_full_details_shows_the_stored_password(self, wired, capsys):
        """Characterisation: --onlyName is the only way to hide credentials."""
        wired.mongoJSON = {"alice": dict(password="secret")}
        CLI.db_list(Namespace(fulldetails=True))
        assert "secret" in capsys.readouterr().out

    def test_no_configured_connection_prints_an_empty_frame(self, wired, capsys):
        wired.mongoJSON = {}
        CLI.db_list(Namespace(fulldetails=True))
        assert "Empty DataFrame" in capsys.readouterr().out


@pytest.mark.unit
class TestDbCreate:
    def test_every_argument_reaches_addorupdatedatabase(self, wired):
        CLI.db_create(Namespace(connectionName="conn", username="user",
                                password="pass", IP="1.2.3.4",
                                databaseName="theDB"))

        assert wired.kwargs("addOrUpdateDatabase") == dict(
            connectionName="conn", username="user", password="pass",
            databaseIP="1.2.3.4", databaseName="theDB")

    def test_the_ip_argument_is_renamed_to_databaseip(self, wired):
        CLI.db_create(Namespace(connectionName="c", username="u", password="p",
                                IP="9.9.9.9", databaseName="d"))
        assert "IP" not in wired.kwargs("addOrUpdateDatabase")
        assert wired.kwargs("addOrUpdateDatabase")["databaseIP"] == "9.9.9.9"

    def test_it_prints_nothing_on_success(self, wired, capsys):
        CLI.db_create(Namespace(connectionName="c", username="u", password="p",
                                IP="1", databaseName="d"))
        assert capsys.readouterr().out == ""


@pytest.mark.unit
class TestDbRemove:
    def test_the_connection_name_is_passed_positionally(self, wired):
        CLI.db_remove(Namespace(connectionName="conn"))
        assert wired.call("removeConnection") == ("removeConnection", ("conn",), {})

    def test_it_prints_nothing_on_success(self, wired, capsys):
        CLI.db_remove(Namespace(connectionName="conn"))
        assert capsys.readouterr().out == ""


# ---------------------------------------------------------------------------
# toolkit commands
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestToolkitList:
    def test_a_document_dict_becomes_one_table_row(self, wired, capsys):
        wired.toolkitDocuments = [dict(toolkit="MyToolkit",
                                       desc=dict(classpath="a.b.C", source="db",
                                                 type="dynamic",
                                                 repositoryName="repoA",
                                                 version=(0, 0, 1)))]
        CLI.toolkit_list(Namespace(project="P"))

        out = capsys.readouterr().out
        assert "MyToolkit" in out and "a.b.C" in out and "repoA" in out
        assert "toolkit" in out and "cls" in out and "source" in out

    def test_the_project_is_forwarded_to_both_the_home_and_the_query(self, wired):
        CLI.toolkit_list(Namespace(project="MY_PROJECT"))

        assert wired.home.projectName == "MY_PROJECT"
        assert wired.home.kwargs("getToolkitDocuments") == dict(name=None,
                                                                projectName="MY_PROJECT")

    def test_a_missing_source_defaults_to_internal(self, wired, capsys):
        wired.toolkitDocuments = [dict(toolkit="MyToolkit", desc={})]
        CLI.toolkit_list(Namespace(project="P"))
        assert "internal" in capsys.readouterr().out

    def test_an_object_document_is_read_through_getattr(self, wired, capsys):
        wired.toolkitDocuments = [_FakeDocument(dict(classpath="x.y.Z"),
                                                toolkit="ObjToolkit")]
        CLI.toolkit_list(Namespace(project="P"))

        out = capsys.readouterr().out
        assert "ObjToolkit" in out and "x.y.Z" in out

    def test_no_documents_prints_an_explanation(self, wired, capsys):
        wired.toolkitDocuments = []
        CLI.toolkit_list(Namespace(project="P"))
        assert "No toolkits found (static + DB are empty)." in capsys.readouterr().out

    def test_a_none_answer_is_treated_as_no_documents(self, wired, capsys):
        wired.toolkitDocuments = None
        CLI.toolkit_list(Namespace(project="P"))
        assert "No toolkits found" in capsys.readouterr().out

    def test_a_failure_inside_the_home_is_reported_not_raised(self, wired, capsys, monkeypatch):
        def _explode(projectName=None, **kwargs):
            raise RuntimeError("no database")

        monkeypatch.setattr(CLI, "ToolkitHome", _explode)
        CLI.toolkit_list(Namespace(project="P"))

        assert "[ERROR] no database" in capsys.readouterr().out

    def test_a_missing_project_argument_is_not_caught(self, wired):
        """``arguments.project`` is read before the try block."""
        with pytest.raises(AttributeError, match="project"):
            CLI.toolkit_list(Namespace())


@pytest.mark.unit
class TestToolkitRegisterIsBroken:
    LOCATABLE = "hera.utils.data.CLI.load_project_name"

    def _arguments(self, **overrides):
        arguments = Namespace(project="P", cls=self.LOCATABLE, name="MyToolkit",
                              repository="repoA", version="1,2,3")
        for key, value in overrides.items():
            setattr(arguments, key, value)
        return arguments

    @pytest.mark.xfail(
        strict=True,
        reason="B179: toolkit_register calls "
               "th.registerToolkit(toolkit_name=toolkit_name, "
               "toolkit_path=toolkit_path, ...) but neither name is ever "
               "defined -- the handler's own values are called cls_path and "
               "name. The NameError is swallowed by the bare `except "
               "Exception`, which prints '[ERROR] name \"toolkit_name\" is "
               "not defined' and returns 0, so the command silently "
               "registers nothing and its --version/--repository parsing is "
               "unreachable. See the consolidated findings issue.",
    )
    def test_registering_a_toolkit_should_reach_toolkithome(self, wired):
        CLI.toolkit_register(self._arguments())
        assert wired.home.kwargs("registerToolkit")["toolkit_name"] == "MyToolkit"

    def test_it_currently_prints_a_nameerror_and_registers_nothing(self, wired, capsys):
        """Characterisation of B179."""
        CLI.toolkit_register(self._arguments())

        out = capsys.readouterr().out
        assert "[ERROR]" in out
        assert "toolkit_name" in out and "is not defined" in out
        assert wired.home.names() == []

    def test_the_module_really_has_no_toolkit_name_global(self):
        """Characterisation of B179: nothing supplies the missing names."""
        assert not hasattr(CLI, "toolkit_name")
        assert not hasattr(CLI, "toolkit_path")

    def test_an_unlocatable_classpath_is_reported_before_the_nameerror(self, wired, capsys):
        CLI.toolkit_register(self._arguments(cls="no.such.module.Class"))

        out = capsys.readouterr().out
        assert "[ERROR] Cannot locate class by classpath: no.such.module.Class" in out
        assert wired.home.names() == []

    def test_the_project_still_reaches_the_toolkithome_constructor(self, wired, capsys):
        CLI.toolkit_register(self._arguments(project="MY_PROJECT"))
        capsys.readouterr()
        assert wired.home.projectName == "MY_PROJECT"

    def test_a_malformed_version_does_not_raise(self, wired, capsys):
        """The (0,0,1) fallback keeps a bad --version from surfacing at all."""
        CLI.toolkit_register(self._arguments(version="not,a,version"))
        assert "[ERROR]" in capsys.readouterr().out


@pytest.mark.unit
class TestToolkitLoad:
    def test_the_toolkit_is_requested_by_name_for_the_project(self, wired, capsys):
        wired.loadedToolkit = _FakeDocument(name="theLoadedToolkit")
        CLI.toolkit_load(Namespace(project="MY_PROJECT", name="MyToolkit"))

        assert wired.home.projectName == "MY_PROJECT"
        assert wired.home.kwargs("getToolkit") == dict(toolkitName="MyToolkit")
        assert "Loaded toolkit: theLoadedToolkit" in capsys.readouterr().out

    def test_a_toolkit_without_a_name_attribute_falls_back_to_the_argument(self, wired, capsys):
        wired.loadedToolkit = object()
        CLI.toolkit_load(Namespace(project="P", name="MyToolkit"))
        assert "Loaded toolkit: MyToolkit" in capsys.readouterr().out

    def test_a_failure_is_reported_not_raised(self, wired, capsys):
        wired.getToolkitError = ValueError("unknown toolkit")
        CLI.toolkit_load(Namespace(project="P", name="Nope"))
        assert "[ERROR] unknown toolkit" in capsys.readouterr().out


@pytest.mark.unit
class TestToolkitDefaultRepoShow:
    def test_the_default_repository_is_printed(self, wired, capsys):
        wired.defaultRepository = "repoA"
        CLI.toolkit_default_repo_show(Namespace(project="MY_PROJECT"))

        assert wired.home.kwargs("getDefaultRepository") == dict(projectName="MY_PROJECT")
        assert capsys.readouterr().out.strip() == "repoA"

    def test_no_default_repository_prints_a_placeholder(self, wired, capsys):
        wired.defaultRepository = None
        CLI.toolkit_default_repo_show(Namespace(project="P"))
        assert capsys.readouterr().out.strip() == "<no default repository>"

    @pytest.mark.parametrize("arguments", [Namespace(), Namespace(project=None)])
    def test_a_missing_project_falls_back_to_defaultproject(self, wired, capsys, arguments):
        CLI.toolkit_default_repo_show(arguments)

        assert wired.home.projectName == "DefaultProject"
        assert wired.home.kwargs("getDefaultRepository") == dict(projectName="DefaultProject")

    def test_a_failure_is_reported_not_raised(self, wired, capsys, monkeypatch):
        monkeypatch.setattr(CLI, "ToolkitHome",
                            lambda **kwargs: (_ for _ in ()).throw(RuntimeError("boom")))
        CLI.toolkit_default_repo_show(Namespace(project="P"))
        assert "[ERROR] boom" in capsys.readouterr().out


@pytest.mark.unit
class TestToolkitDefaultRepoSet:
    def test_the_repository_is_set_for_the_project(self, wired, capsys):
        CLI.toolkit_default_repo_set(Namespace(project="MY_PROJECT", repository="repoA"))

        assert wired.home.kwargs("setDefaultRepository") == dict(
            projectName="MY_PROJECT", repositoryName="repoA")
        assert "Default repository set to 'repoA' for project 'MY_PROJECT'." in \
            capsys.readouterr().out

    @pytest.mark.parametrize("repository", [None, ""])
    def test_a_missing_repository_is_refused_before_any_work(self, wired, capsys, repository):
        CLI.toolkit_default_repo_set(Namespace(project="P", repository=repository))

        assert capsys.readouterr().out.strip() == "[ERROR] --repository is required"
        assert wired.homes == []

    def test_a_missing_repository_attribute_is_refused_too(self, wired, capsys):
        CLI.toolkit_default_repo_set(Namespace(project="P"))
        assert "[ERROR] --repository is required" in capsys.readouterr().out

    def test_a_missing_project_falls_back_to_defaultproject(self, wired, capsys):
        CLI.toolkit_default_repo_set(Namespace(repository="repoA"))
        assert wired.home.kwargs("setDefaultRepository")["projectName"] == "DefaultProject"


@pytest.mark.unit
class TestToolkitImportJson:
    def _arguments(self, **overrides):
        arguments = Namespace(project="MY_PROJECT", file="/repo/toolkits.json",
                              no_experiments=False)
        for key, value in overrides.items():
            setattr(arguments, key, value)
        return arguments

    def test_the_registered_toolkits_are_printed(self, wired, capsys):
        wired.importedToolkits = ["MyToolkit", "Other"]
        CLI.toolkit_import_json(self._arguments())

        assert wired.home.kwargs("import_toolkits_from_json") == dict(
            projectName="MY_PROJECT", json_path="/repo/toolkits.json")
        assert "Registered toolkits: ['MyToolkit', 'Other']" in capsys.readouterr().out

    def test_an_empty_result_says_there_were_no_toolkits(self, wired, capsys):
        wired.importedToolkits = []
        CLI.toolkit_import_json(self._arguments())
        assert "No toolkits in JSON." in capsys.readouterr().out

    def test_the_experiments_are_imported_by_default(self, wired, capsys):
        wired.importedExperiments = ["expA"]
        CLI.toolkit_import_json(self._arguments())

        assert wired.home.kwargs("import_experiments_from_json") == dict(
            projectName="MY_PROJECT", json_path="/repo/toolkits.json")
        assert "Loaded experiments data: ['expA']" in capsys.readouterr().out

    def test_no_experiments_skips_the_second_import(self, wired, capsys):
        CLI.toolkit_import_json(self._arguments(no_experiments=True))
        assert "import_experiments_from_json" not in wired.home.names()

    def test_an_empty_experiment_result_is_not_announced(self, wired, capsys):
        wired.importedExperiments = []
        CLI.toolkit_import_json(self._arguments())

        assert "import_experiments_from_json" in wired.home.names()
        assert "Loaded experiments data" not in capsys.readouterr().out

    def test_missing_attributes_become_none_and_false(self, wired):
        CLI.toolkit_import_json(Namespace())

        assert wired.home.projectName is None
        assert wired.home.kwargs("import_toolkits_from_json") == dict(projectName=None,
                                                                      json_path=None)
        assert "import_experiments_from_json" in wired.home.names()

    def test_a_failure_is_reported_not_raised(self, wired, capsys, monkeypatch):
        monkeypatch.setattr(CLI, "ToolkitHome",
                            lambda **kwargs: (_ for _ in ()).throw(RuntimeError("bad json")))
        CLI.toolkit_import_json(self._arguments())
        assert "[ERROR] bad json" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# project measurements / simulations / cache list
# ---------------------------------------------------------------------------

@pytest.fixture()
def measurementDocuments(wired, monkeypatch):
    """``project_measurements_list`` imports Project inside the function, so
    the patch has to land on ``hera.datalayer`` rather than on the CLI
    module's global.
    """
    import hera.datalayer

    monkeypatch.setattr(hera.datalayer, "Project", wired.Project)
    return wired


@pytest.mark.unit
class TestProjectMeasurementsList:
    def _document(self, **overrides):
        fields = dict(type="ToolkitDataSource", datasourceName="YAVNEEL",
                      resource="/data/yavneel.parquet", dataFormat="parquet")
        fields.update(overrides)
        desc = dict(version=(0, 0, 1), repository="meteo_data_v1",
                    datasourceName=fields["datasourceName"])
        return _FakeDocument(desc, **fields)

    def _arguments(self, **overrides):
        arguments = Namespace(project="MY_PROJECT", type=None, shortcut=None,
                              contains=None)
        for key, value in overrides.items():
            setattr(arguments, key, value)
        return arguments

    def test_it_prints_a_markdown_table_of_the_documents(self, measurementDocuments, capsys):
        measurementDocuments.measurementsDocuments = [self._document()]
        CLI.project_measurements_list(self._arguments())

        out = capsys.readouterr().out
        assert "|" in out
        assert "datasourceName" in out and "YAVNEEL" in out
        assert "/data/yavneel.parquet" in out and "meteo_data_v1" in out

    def test_the_project_is_opened_by_name(self, measurementDocuments):
        CLI.project_measurements_list(self._arguments())
        assert measurementDocuments.project.buildArguments == {}
        assert measurementDocuments.project.projectName == "MY_PROJECT"

    def test_without_a_project_the_default_project_is_opened(self, measurementDocuments):
        CLI.project_measurements_list(self._arguments(project=None))
        assert measurementDocuments.project.projectName is None

    def test_no_type_filter_queries_everything(self, measurementDocuments):
        CLI.project_measurements_list(self._arguments())
        assert measurementDocuments.project.kwargs("getMeasurementsDocuments") == {}

    def test_an_explicit_type_becomes_the_query_filter(self, measurementDocuments):
        CLI.project_measurements_list(self._arguments(type="Experiment_rawData"))
        assert measurementDocuments.project.kwargs("getMeasurementsDocuments") == \
            dict(type="Experiment_rawData")

    @pytest.mark.parametrize("shortcut,expected", [
        ("ds", "ToolkitDataSource"),
        ("exp", "Experiment_rawData"),
        ("sim", "Simulation_rawData"),
        ("cache", "Cache_rawData"),
    ])
    def test_each_shortcut_maps_to_its_document_type(self, measurementDocuments,
                                                     shortcut, expected):
        CLI.project_measurements_list(self._arguments(shortcut=shortcut))
        assert measurementDocuments.project.kwargs("getMeasurementsDocuments") == \
            dict(type=expected)

    def test_the_all_shortcut_removes_the_type_filter(self, measurementDocuments):
        CLI.project_measurements_list(self._arguments(shortcut="all", type="Whatever"))
        assert measurementDocuments.project.kwargs("getMeasurementsDocuments") == {}

    def test_a_shortcut_overrides_an_explicit_type(self, measurementDocuments):
        CLI.project_measurements_list(self._arguments(shortcut="ds", type="Whatever"))
        assert measurementDocuments.project.kwargs("getMeasurementsDocuments") == \
            dict(type="ToolkitDataSource")

    def test_an_unknown_shortcut_is_refused_and_lists_the_valid_ones(self, measurementDocuments, capsys):
        CLI.project_measurements_list(self._arguments(shortcut="nope"))

        out = capsys.readouterr().out
        assert "Unknown shortcut 'nope'. Valid: ds, exp, sim, cache, all" in out
        assert measurementDocuments.projects == []

    def test_the_contains_filter_matches_the_datasource_name(self, measurementDocuments, capsys):
        measurementDocuments.measurementsDocuments = [
            self._document(datasourceName="YAVNEEL"),
            self._document(datasourceName="TEL_AVIV"),
        ]
        CLI.project_measurements_list(self._arguments(contains="AVIV"))

        out = capsys.readouterr().out
        assert "TEL_AVIV" in out and "YAVNEEL" not in out

    def test_the_contains_filter_also_matches_the_resource(self, measurementDocuments, capsys):
        measurementDocuments.measurementsDocuments = [
            self._document(datasourceName="A", resource="/data/keep.parquet"),
            self._document(datasourceName="B", resource="/data/drop.parquet"),
        ]
        CLI.project_measurements_list(self._arguments(contains="keep"))

        out = capsys.readouterr().out
        assert "keep.parquet" in out and "drop.parquet" not in out

    def test_the_contains_filter_is_case_sensitive(self, measurementDocuments, capsys):
        measurementDocuments.measurementsDocuments = [
            self._document(datasourceName="YAVNEEL", resource="/data/STATION.parquet")]
        CLI.project_measurements_list(self._arguments(contains="yavneel"))

        assert "No measurements found for given filters." in capsys.readouterr().out

    def test_the_name_falls_back_to_the_description_when_the_attribute_is_empty(
            self, measurementDocuments, capsys):
        document = _FakeDocument(dict(datasourceName="FROM_DESC"), type="T",
                                 datasourceName="", resource="", dataFormat="")
        measurementDocuments.measurementsDocuments = [document]
        CLI.project_measurements_list(self._arguments(contains="FROM_DESC"))

        assert "FROM_DESC" in capsys.readouterr().out

    def test_an_empty_project_says_no_measurements_were_found(self, measurementDocuments, capsys):
        measurementDocuments.measurementsDocuments = []
        CLI.project_measurements_list(self._arguments())
        assert "No measurements found for given filters." in capsys.readouterr().out

    def test_missing_attributes_are_all_tolerated(self, measurementDocuments, capsys):
        measurementDocuments.measurementsDocuments = [self._document()]
        CLI.project_measurements_list(Namespace())

        assert measurementDocuments.project.projectName is None
        assert "YAVNEEL" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# fake-drift guards
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestTheFakesDoNotDrift:
    def test_the_data_toolkit_really_has_every_faked_method(self):
        from hera.utils.data.toolkit import dataToolkit

        for name in ("getRepositoryTable", "addRepository", "deleteDataSource",
                     "getDataSourceData", "getDataSourceList",
                     "loadAllDatasourcesInRepositoryJSONToProject",
                     "loadAllDatasourcesInAllRepositoriesToProject",
                     "exportDocumentsToRepository"):
            assert hasattr(dataToolkit, name), f"dataToolkit lost {name}"

    def test_the_toolkit_home_really_has_every_faked_method(self):
        from hera.toolkit import ToolkitHome

        for name in ("getToolkitDocuments", "registerToolkit", "getToolkit",
                     "getDefaultRepository", "setDefaultRepository",
                     "import_toolkits_from_json", "import_experiments_from_json"):
            assert hasattr(ToolkitHome, name), f"ToolkitHome lost {name}"

    def test_the_project_really_has_every_faked_method(self):
        from hera.datalayer import Project

        for name in ("getCacheDocuments", "getMeasurementsDocuments",
                     "getSimulationsDocuments", "getMeasurementsDocumentsAsDict",
                     "getConfig", "addDocumentFromDict"):
            assert hasattr(Project, name), f"Project lost {name}"

    def test_the_default_project_name_matches_the_real_one(self):
        from hera.datalayer import Project

        assert _makeProjectClass(_Harness()).DEFAULTPROJECT == Project.DEFAULTPROJECT

    def test_the_datalayer_all_collection_really_has_getdocuments(self):
        from hera.datalayer import All

        assert hasattr(All, "getDocuments")
