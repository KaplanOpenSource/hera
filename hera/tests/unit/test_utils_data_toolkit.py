"""utils/data/toolkit.py: dataToolkit -- the repository registry that backs
``hera-data``, plus the repository-JSON loader that fans a repository out into
a project (Config / Datasource / Measurements / Simulations / Cache / Function
sections) and the reverse export path.

Construction strategy
---------------------
``dataToolkit.__init__`` hardcodes ``projectName=self.DEFAULTPROJECT`` and
``filesDirectory=None`` and accepts no ``filesDirectory`` argument, so a naive
construction in a normal process would reach the developer's real MongoDB
default project.  Here it is built *directly* -- ``dataToolkit()`` -- and that
is safe for two independent reasons, both verified by TestConstruction below:

1. This directory's conftest reroutes ``hera.datalayer`` at mongomock before
   any test module is imported, so ``"defaultProject"`` is a project in an
   in-memory database that ``_reset_unit_database`` drops after every test.
   Nothing reaches a real mongod -- ``_block_network`` fails any socket call.
2. ``Project.__init__`` short-circuits every disk operation when
   ``projectName == DEFAULTPROJECT`` (hera/datalayer/project.py:262-284): it
   neither reads nor writes the project config and never runs the
   ``os.makedirs(~/.hera/<projectName>)`` branch, so ``FilesDirectory`` stays
   ``None`` and no directory is created.  HOME is a temp directory anyway, and
   the conftest's ``_no_home_writes`` guard fails the test if ``~/.hera``
   appears.

The one place that *would* write to a real home is the loader, which builds
``ToolkitHome(projectName=<target>)`` and then a target toolkit -- a non-default
project, so ``~/.hera/<projectName>`` would be created.  ``ToolkitHome`` is
therefore monkeypatched **in the module under test's own namespace**
(``hera.utils.data.toolkit.ToolkitHome``, a module attribute -- never a method
on a live class, and never a singleton instance) with a stand-in that hands the
loader a real ``abstractToolkit`` bound to ``unit_files_directory``.  One test,
``TestRealToolkitHomeWiring``, uses the *real* ``ToolkitHome`` (through a
subclass that only injects ``filesDirectory``) to show the production wiring is
genuinely reached.  Likewise the export path's ``Project`` is patched at
``hera.datalayer.project.Project`` -- the module attribute its local import
reads -- to default ``filesDirectory``.

Covered
-------
``__init__``, ``addRepository``, ``getRepositoryTable``, ``getRepository``,
``loadAllDatasourcesInAllRepositoriesToProject``,
``loadAllDatasourcesInRepositoryToProject``, ``getToolkitDocument``,
``loadAllDatasourcesInRepositoryJSONToProject``, ``_handle_Config``,
``_handle_DataSource``, ``_DocumentHandler``, ``_handle_Function``,
``_makeItemPathAbsolute``, ``resolveDataSourcePaths``,
``loadRepositoryFromPath``, ``exportDocumentsToRepository``,
``_resolveDocumentsForExport``, ``_resolveRepositoryPath``.

Bugs pinned here (xfail(strict=True) for the intended behaviour + a passing
characterisation of what actually happens):

* B151 ``addRepository`` decides the ``.json`` extension with
  ``"json" not in repositoryPath`` -- a substring test over the whole path.
* B152 ``addRepository`` re-arms the default project's read-only guard only
  on the success path; a raising ``addDataSource`` leaves it disarmed.
* B153 the ``Cache`` entry of the handler table is a lambda whose third
  parameter is named ``itemDesc`` while the dispatcher calls every handler with
  ``docTypeDict=`` -- Cache sections can never load.
* B154 ``_makeItemPathAbsolute`` coerces ``isRelativePath`` with ``bool()``,
  so the string ``"False"`` that the sibling handler *requires* is truthy.
* B159 ``_DocumentHandler`` never removes the ``isRelativePath`` control key
  from the item, so it is forwarded as a keyword to ``add<Type>Document`` and
  raises ``TypeError`` -- which is also why B154 is only reachable by calling
  ``_makeItemPathAbsolute`` directly.
* B155 ``_DocumentHandler`` guards ``getattr(toolkit, ...)`` with
  ``if retrieveFunc is None`` -- unreachable, ``getattr`` raises first.
* B156 ``_handle_Function`` resolves the method name on ``self`` (the
  dataToolkit, on the read-only default project) instead of on the toolkit the
  section belongs to.
* B157 ``_handle_DataSource`` pops ``resourceFilePath`` inside the ``open()``
  call, so a failed read leaves the item with neither key.
* B158 ``_resolveDocumentsForExport`` guards ``getDocumentByID`` with
  ``if doc is None`` -- unreachable, mongoengine raises ``DoesNotExist``.

Deliberately not covered
------------------------
* The dynamic-toolkit auto-registration path of
  ``loadAllDatasourcesInRepositoryJSONToProject`` -- the ``auto_register_missing``
  argument is documented but its implementation was deleted (the numbered
  comments jump from step 1 to step 3), so there is nothing to exercise beyond
  the "unknown toolkit is skipped" test that is present.
* ``exportDocumentsToRepository(mode="override")``'s dedup pass and the
  ``idStrategy`` variants: those live in ``repositoryExport`` and are already
  covered by test_utils_data_repository_export.py; only the orchestration is
  tested here.
* Real data payloads (parquet/netcdf resources are strings on disk that are
  never opened by this module).
"""
import json
import os

import pytest

from hera.datalayer import Project
from hera.toolkit import abstractToolkit
from hera.utils.data import toolkit as dataToolkitModule
from hera.utils.data.toolkit import dataToolkit

TARGET_TOOLKIT_NAME = "known"


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def make_target_toolkit(unit_files_directory):
    """Build a real abstractToolkit in a named project, files under tmp_path."""

    def _build(projectName, toolkitName=TARGET_TOOLKIT_NAME):
        return abstractToolkit(
            toolkitName=toolkitName,
            projectName=projectName,
            filesDirectory=unit_files_directory,
        )

    return _build


@pytest.fixture()
def dtk(monkeypatch, make_target_toolkit):
    """A dataToolkit whose loader hands out real, temp-backed toolkits.

    ToolkitHome is replaced as a *module attribute* of the module under test
    (that is where the loader looks it up), not as an attribute of any live
    object, so monkeypatch's teardown restores the real class cleanly.
    """

    class _StandInToolkitHome:
        def __init__(self, projectName=None, filesDirectory=None):
            self.projectName = projectName

        def getToolkit(self, toolkitName, filesDirectory=None, **kwargs):
            if toolkitName != TARGET_TOOLKIT_NAME:
                raise ValueError(f"toolkit {toolkitName} is not registered")
            return make_target_toolkit(self.projectName, toolkitName=toolkitName)

    monkeypatch.setattr(dataToolkitModule, "ToolkitHome", _StandInToolkitHome)
    return dataToolkit()


@pytest.fixture()
def patched_datalayer_project(monkeypatch, unit_files_directory):
    """Default filesDirectory for the Project the export path builds itself.

    ``exportDocumentsToRepository`` does ``from hera.datalayer.project import
    Project`` inside the method body, so the name is read off the module at
    call time -- patching the module attribute is enough and is undone on
    teardown.
    """
    import hera.datalayer.project as projectModule

    realProject = projectModule.Project

    def _projectWithFilesDirectory(projectName=None, **kwargs):
        kwargs.setdefault("filesDirectory", unit_files_directory)
        return realProject(projectName=projectName, **kwargs)

    monkeypatch.setattr(projectModule, "Project", _projectWithFilesDirectory)
    return _projectWithFilesDirectory


def _writeRepository(directory, name, payload):
    """Write a repository JSON file and return its path."""
    path = os.path.join(str(directory), f"{name}.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle)
    return path


def _datasourceSection(resource, isRelativePath="False", dataFormat="string", **extra):
    return {"isRelativePath": isRelativePath,
            "item": dict(resource=resource, dataFormat=dataFormat, **extra)}


# ---------------------------------------------------------------------------
# construction
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestConstruction:
    def test_it_always_binds_to_the_default_project(self):
        assert dataToolkit().projectName == Project.DEFAULTPROJECT

    def test_its_toolkit_name_is_heradata(self):
        assert dataToolkit().toolkitName == "heradata"

    def test_the_default_project_gets_no_files_directory_on_disk(self):
        """The safety property this whole module leans on: the default
        project never triggers the ~/.hera/<projectName> makedirs branch."""
        assert dataToolkit().FilesDirectory is None

    def test_the_default_project_starts_read_only(self):
        assert dataToolkit()._allowWritingToDefaultProject is False

    def test_writing_to_the_default_project_is_refused_by_default(self):
        with pytest.raises(RuntimeError, match="read-only"):
            dataToolkit().addMeasurementsDocument(
                resource="x", dataFormat="string", type="T", desc={}
            )


# ---------------------------------------------------------------------------
# addRepository / getRepositoryTable / getRepository
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestAddRepository:
    def test_it_registers_the_repository_as_a_json_dict_data_source(self, dtk, tmp_path):
        path = _writeRepository(tmp_path, "repo", {})
        dtk.addRepository("myRepo", path)
        document = dtk.getDataSourceDocument("myRepo")
        assert document.dataFormat == "JSON_dict"
        assert document.desc["datasourceName"] == "myRepo"

    def test_the_resource_is_stored_as_an_absolute_path(self, dtk, tmp_path, monkeypatch):
        _writeRepository(tmp_path, "repo", {})
        monkeypatch.chdir(tmp_path)
        dtk.addRepository("myRepo", "repo.json")
        assert dtk.getDataSourceDocument("myRepo").resource == str(tmp_path / "repo.json")

    def test_a_missing_extension_is_appended(self, dtk, tmp_path):
        """The directory must not itself contain the substring 'json' --
        see B151 below; pytest's own tmp_path is named after the test, so
        naming this test '..._json_...' would break it."""
        _writeRepository(tmp_path, "repo", {})
        dtk.addRepository("myRepo", str(tmp_path / "repo"))
        assert dtk.getDataSourceDocument("myRepo").resource.endswith("repo.json")

    def test_re_registering_without_overwrite_raises(self, dtk, tmp_path):
        path = _writeRepository(tmp_path, "repo", {})
        dtk.addRepository("myRepo", path)
        with pytest.raises(ValueError, match="already exists"):
            dtk.addRepository("myRepo", path)

    def test_overwrite_replaces_the_registration(self, dtk, tmp_path):
        first = _writeRepository(tmp_path, "one", {})
        second = _writeRepository(tmp_path, "two", {})
        dtk.addRepository("myRepo", first)
        dtk.addRepository("myRepo", second, overwrite=True)
        assert dtk.getDataSourceDocument("myRepo").resource == second
        assert dtk.getDataSourceList() == ["myRepo"]

    def test_it_re_arms_the_default_project_guard_afterwards(self, dtk, tmp_path):
        dtk.addRepository("myRepo", _writeRepository(tmp_path, "repo", {}))
        assert dtk._allowWritingToDefaultProject is False


@pytest.mark.unit
class TestAddRepositoryExtensionHeuristic:
    """B151: the extension is decided with a substring test on the path."""

    @pytest.mark.xfail(
        strict=True,
        reason="B151: addRepository appends '.json' only when the substring "
               "'json' is absent from the WHOLE path, so a repository living "
               "under a directory whose name contains 'json' is registered "
               "with no extension and its resource points at a file that does "
               "not exist. The test should be on the suffix "
               "(os.path.splitext / endswith('.json')). "
               "See the consolidated findings issue.",
    )
    def test_the_extension_is_appended_regardless_of_the_directory_name(self, dtk, tmp_path):
        directory = tmp_path / "jsonRepositories"
        directory.mkdir()
        _writeRepository(directory, "repo", {})
        dtk.addRepository("myRepo", str(directory / "repo"))
        assert dtk.getDataSourceDocument("myRepo").resource.endswith("repo.json")

    def test_a_json_in_the_directory_name_suppresses_the_extension(self, dtk, tmp_path):
        """Characterisation of B151."""
        directory = tmp_path / "jsonRepositories"
        directory.mkdir()
        _writeRepository(directory, "repo", {})
        dtk.addRepository("myRepo", str(directory / "repo"))
        resource = dtk.getDataSourceDocument("myRepo").resource
        assert resource == str(directory / "repo")
        assert not os.path.exists(resource)


@pytest.mark.unit
class TestAddRepositoryLeavesTheDefaultProjectWritable:
    """B152: the read-only guard is re-armed only on the success path."""

    @pytest.mark.xfail(
        strict=True,
        reason="B152: addRepository sets _allowWritingToDefaultProject=True, "
               "calls addDataSource, then sets it back to False. There is no "
               "try/finally, so when addDataSource raises (a duplicate "
               "repository name without overwrite=True) the flag stays True "
               "and the read-only default project accepts arbitrary writes for "
               "the rest of the instance's life. "
               "See the consolidated findings issue.",
    )
    def test_a_failed_registration_still_re_arms_the_guard(self, dtk, tmp_path):
        path = _writeRepository(tmp_path, "repo", {})
        dtk.addRepository("myRepo", path)
        with pytest.raises(ValueError):
            dtk.addRepository("myRepo", path)
        assert dtk._allowWritingToDefaultProject is False

    def test_a_failed_registration_disarms_the_guard_permanently(self, dtk, tmp_path):
        """Characterisation of B152."""
        path = _writeRepository(tmp_path, "repo", {})
        dtk.addRepository("myRepo", path)
        with pytest.raises(ValueError):
            dtk.addRepository("myRepo", path)
        assert dtk._allowWritingToDefaultProject is True
        # and the guard that raised RuntimeError in TestConstruction is gone
        document = dtk.addMeasurementsDocument(
            resource="sneaked-in", dataFormat="string", type="NOT_A_DATASOURCE", desc={}
        )
        assert document is not None


@pytest.mark.unit
class TestGetRepositoryTable:
    def test_an_empty_registry_yields_an_empty_table(self, dtk):
        assert len(dtk.getRepositoryTable()) == 0

    def test_each_repository_becomes_one_row(self, dtk, tmp_path):
        dtk.addRepository("one", _writeRepository(tmp_path, "one", {}))
        dtk.addRepository("two", _writeRepository(tmp_path, "two", {}))
        table = dtk.getRepositoryTable()
        assert sorted(table["datasourceName"]) == ["one", "two"]

    def test_the_table_reports_the_json_dict_format(self, dtk, tmp_path):
        dtk.addRepository("one", _writeRepository(tmp_path, "one", {}))
        assert list(dtk.getRepositoryTable()["dataFormat"]) == ["JSON_dict"]


@pytest.mark.unit
class TestGetRepository:
    def test_it_returns_the_parsed_json_of_the_registered_file(self, dtk, tmp_path):
        payload = {TARGET_TOOLKIT_NAME: {"Config": {"stationType": "IMS"}}}
        dtk.addRepository("myRepo", _writeRepository(tmp_path, "repo", payload))
        assert dtk.getRepository("myRepo") == payload

    def test_an_unregistered_name_raises_valueerror_from_the_json_loader(self, dtk):
        """getDataSourceData returns None, which loadJSON rejects -- the error
        names NoneType rather than the missing repository."""
        with pytest.raises(ValueError, match="NoneType"):
            dtk.getRepository("noSuchRepository")


@pytest.mark.unit
class TestGetToolkitDocument:
    def test_a_registered_repository_is_found_by_its_name(self, dtk, tmp_path):
        dtk.addRepository("myRepo", _writeRepository(tmp_path, "repo", {}))
        found = dtk.getToolkitDocument("myRepo")
        assert found is not None
        assert found.desc["datasourceName"] == "myRepo"

    def test_an_unknown_name_returns_none(self, dtk):
        assert dtk.getToolkitDocument("noSuchToolkit") is None


# ---------------------------------------------------------------------------
# loadAllDatasourcesInRepositoryJSONToProject: input handling
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestRepositoryJSONInputHandling:
    def test_a_path_like_string_is_skipped(self, dtk):
        assert dtk.loadAllDatasourcesInRepositoryJSONToProject(
            projectName="P", repositoryJSON="/some/where/repo.json"
        ) is None

    def test_a_json_string_is_parsed_and_loaded(self, dtk, make_target_toolkit):
        dtk.loadAllDatasourcesInRepositoryJSONToProject(
            projectName="JSON_STRING_PROJECT",
            repositoryJSON=json.dumps({TARGET_TOOLKIT_NAME: {"Config": {"a": 1}}}),
        )
        assert make_target_toolkit("JSON_STRING_PROJECT").getConfig()["a"] == 1

    def test_an_unparseable_string_is_reported_and_ignored(self, dtk):
        assert dtk.loadAllDatasourcesInRepositoryJSONToProject(
            projectName="P", repositoryJSON="this is not json"
        ) is None

    def test_a_non_dict_payload_is_ignored(self, dtk):
        assert dtk.loadAllDatasourcesInRepositoryJSONToProject(
            projectName="P", repositoryJSON=["a", "list"]
        ) is None

    def test_an_empty_dict_is_ignored(self, dtk):
        assert dtk.loadAllDatasourcesInRepositoryJSONToProject(
            projectName="P", repositoryJSON={}
        ) is None

    def test_a_toolkit_that_cannot_be_resolved_is_skipped_quietly(self, dtk):
        """An unresolvable toolkit key short-circuits before section
        dispatch, so even a bogus section name raises nothing."""
        assert dtk.loadAllDatasourcesInRepositoryJSONToProject(
            projectName="P", repositoryJSON={"noSuchToolkit": {"Bogus": {}}}
        ) is None

    def test_an_unknown_section_name_raises(self, dtk):
        with pytest.raises(ValueError, match="Handler"):
            dtk.loadAllDatasourcesInRepositoryJSONToProject(
                projectName="P", repositoryJSON={TARGET_TOOLKIT_NAME: {"Bogus": {}}}
            )

    def test_section_names_are_matched_case_insensitively_via_title(self, dtk, make_target_toolkit):
        """'DataSource'.title() == 'Datasource', which is the table key."""
        dtk.loadAllDatasourcesInRepositoryJSONToProject(
            projectName="TITLE_PROJECT",
            repositoryJSON={TARGET_TOOLKIT_NAME: {
                "DataSource": {"A": _datasourceSection("/absolute/resource")}
            }},
        )
        names = make_target_toolkit("TITLE_PROJECT").getDataSourceList()
        assert names == ["A"]

    def test_a_failing_section_is_swallowed_and_only_logged(self, dtk, make_target_toolkit):
        """Every handler call is wrapped in `except Exception: logger.error`,
        so a broken section leaves the loader looking successful."""
        dtk.loadAllDatasourcesInRepositoryJSONToProject(
            projectName="SWALLOW_PROJECT",
            repositoryJSON={TARGET_TOOLKIT_NAME: {
                # no isRelativePath -> the handler's assert fires
                "Datasource": {"A": {"item": {"resource": "/r", "dataFormat": "string"}}}
            }},
        )
        assert make_target_toolkit("SWALLOW_PROJECT").getDataSourceList() == []


# ---------------------------------------------------------------------------
# _handle_Config
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestHandleConfig:
    def test_a_config_section_lands_in_the_target_project_config(self, dtk, make_target_toolkit):
        dtk.loadAllDatasourcesInRepositoryJSONToProject(
            projectName="CONFIG_PROJECT",
            repositoryJSON={TARGET_TOOLKIT_NAME: {"Config": {"stationType": "IMS"}}},
        )
        assert make_target_toolkit("CONFIG_PROJECT").getConfig()["stationType"] == "IMS"

    def test_it_forwards_every_key_as_a_keyword_argument(self, dtk, make_target_toolkit):
        target = make_target_toolkit("CONFIG_PROJECT_2")
        dtk._handle_Config(target, "Config", {"a": 1, "b": "two"}, False, "")
        config = target.getConfig()
        assert (config["a"], config["b"]) == (1, "two")


# ---------------------------------------------------------------------------
# _handle_DataSource
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestHandleDataSource:
    def test_a_relative_resource_is_joined_onto_the_basedir(self, dtk, make_target_toolkit, tmp_path):
        target = make_target_toolkit("DS_PROJECT_1")
        dtk._handle_DataSource(
            target, "Datasource",
            {"A": _datasourceSection("nested/file.parquet", isRelativePath="True",
                                     dataFormat="parquet")},
            False, str(tmp_path),
        )
        assert target.getDataSourceDocument("A").resource == str(tmp_path / "nested/file.parquet")

    def test_an_absolute_resource_is_left_alone(self, dtk, make_target_toolkit, tmp_path):
        target = make_target_toolkit("DS_PROJECT_2")
        dtk._handle_DataSource(
            target, "Datasource", {"A": _datasourceSection("/already/absolute")},
            False, str(tmp_path),
        )
        assert target.getDataSourceDocument("A").resource == "/already/absolute"

    def test_a_boolean_true_flag_is_accepted_as_well_as_the_string(self, dtk, make_target_toolkit, tmp_path):
        target = make_target_toolkit("DS_PROJECT_3")
        dtk._handle_DataSource(
            target, "Datasource", {"A": _datasourceSection("rel", isRelativePath=True)},
            False, str(tmp_path),
        )
        assert target.getDataSourceDocument("A").resource == str(tmp_path / "rel")

    def test_a_missing_isrelativepath_flag_is_rejected(self, dtk, make_target_toolkit):
        target = make_target_toolkit("DS_PROJECT_4")
        with pytest.raises(AssertionError, match="isRelativePath"):
            dtk._handle_DataSource(
                target, "Datasource",
                {"A": {"item": {"resource": "/r", "dataFormat": "string"}}}, False, "",
            )

    def test_extra_item_keys_become_description_fields(self, dtk, make_target_toolkit):
        target = make_target_toolkit("DS_PROJECT_5")
        dtk._handle_DataSource(
            target, "Datasource",
            {"A": _datasourceSection("/r", desc={"station": "YAVNEEL"})}, False, "",
        )
        assert target.getDataSourceDocument("A").desc["desc"] == {"station": "YAVNEEL"}

    def test_a_duplicate_without_overwrite_keeps_the_first_registration(self, dtk, make_target_toolkit):
        target = make_target_toolkit("DS_PROJECT_6")
        dtk._handle_DataSource(target, "Datasource", {"A": _datasourceSection("/one")}, False, "")
        dtk._handle_DataSource(target, "Datasource", {"A": _datasourceSection("/two")}, False, "")
        assert [d["resource"] for d in target.getDataSourceMap()] == ["/one"]

    def test_overwrite_replaces_the_existing_registration(self, dtk, make_target_toolkit):
        target = make_target_toolkit("DS_PROJECT_7")
        dtk._handle_DataSource(target, "Datasource", {"A": _datasourceSection("/one")}, False, "")
        dtk._handle_DataSource(target, "Datasource", {"A": _datasourceSection("/two")}, True, "")
        assert [d["resource"] for d in target.getDataSourceMap()] == ["/two"]

    def test_resourcefilepath_stores_the_file_contents_as_the_resource(self, dtk, make_target_toolkit, tmp_path):
        with open(tmp_path / "payload.json", "w", encoding="utf-8") as handle:
            json.dump({"k": "v"}, handle)
        target = make_target_toolkit("DS_PROJECT_8")
        dtk._handle_DataSource(
            target, "Datasource",
            {"A": {"isRelativePath": "True",
                   "item": {"resourceFilePath": "payload.json", "dataFormat": "JSON_dict"}}},
            False, str(tmp_path),
        )
        assert target.getDataSourceDocument("A").resource == {"k": "v"}

    def test_resource_wins_when_resourcefilepath_is_also_given(self, dtk, make_target_toolkit, tmp_path):
        target = make_target_toolkit("DS_PROJECT_9")
        dtk._handle_DataSource(
            target, "Datasource",
            {"A": {"isRelativePath": "False",
                   "item": {"resource": "/the/resource",
                            "resourceFilePath": "/the/other/file.json",
                            "dataFormat": "string"}}},
            False, str(tmp_path),
        )
        assert target.getDataSourceDocument("A").resource == "/the/resource"


@pytest.mark.unit
class TestHandleDataSourcePopsResourceFilePathTooEarly:
    """B157: the pop happens inside the open() call, before the read."""

    @pytest.mark.xfail(
        strict=True,
        reason="B157: _handle_DataSource does "
               "`with open(theItem.pop('resourceFilePath')) as f:` -- the key "
               "is removed as the arguments are evaluated, so when the file "
               "cannot be read the `except` logs the real cause and execution "
               "falls through to addDataSource(**theItem) with neither "
               "'resource' nor 'resourceFilePath' present. The caller sees "
               "TypeError: missing 'resource' instead of the FileNotFoundError. "
               "See the consolidated findings issue.",
    )
    def test_an_unreadable_resource_file_reports_the_read_failure(self, dtk, make_target_toolkit, tmp_path):
        target = make_target_toolkit("DS_POP_1")
        with pytest.raises(FileNotFoundError):
            dtk._handle_DataSource(
                target, "Datasource",
                {"A": {"isRelativePath": "True",
                       "item": {"resourceFilePath": "missing.json", "dataFormat": "JSON_dict"}}},
                False, str(tmp_path),
            )

    def test_an_unreadable_resource_file_surfaces_as_a_missing_argument(self, dtk, make_target_toolkit, tmp_path):
        """Characterisation of B157."""
        target = make_target_toolkit("DS_POP_2")
        with pytest.raises(TypeError, match="resource"):
            dtk._handle_DataSource(
                target, "Datasource",
                {"A": {"isRelativePath": "True",
                       "item": {"resourceFilePath": "missing.json", "dataFormat": "JSON_dict"}}},
                False, str(tmp_path),
            )


# ---------------------------------------------------------------------------
# _DocumentHandler (Measurements / Simulations / Cache)
# ---------------------------------------------------------------------------

def _documentSection(resource, documentType="MY_TYPE", desc=None, dataFormat="parquet"):
    return {"D": {"item": {"resource": resource, "dataFormat": dataFormat,
                           "type": documentType, "desc": dict(desc or {})}}}


@pytest.mark.unit
class TestDocumentHandler:
    def test_a_measurements_section_adds_a_measurements_document(self, dtk, make_target_toolkit, tmp_path):
        dtk.loadAllDatasourcesInRepositoryJSONToProject(
            projectName="DOC_PROJECT_1",
            repositoryJSON={TARGET_TOOLKIT_NAME: {
                "Measurements": _documentSection("data/x.parquet", desc={"station": "A"})}},
            basedir=str(tmp_path),
        )
        documents = make_target_toolkit("DOC_PROJECT_1").getMeasurementsDocuments(type="MY_TYPE")
        assert len(documents) == 1
        assert documents[0].resource == str(tmp_path / "data/x.parquet")

    def test_a_simulations_section_adds_a_simulations_document(self, dtk, make_target_toolkit):
        target = make_target_toolkit("DOC_PROJECT_2")
        dtk._DocumentHandler(target, "Simulations", _documentSection("/s.parquet"),
                             False, "Simulations", "")
        assert len(target.getSimulationsDocuments(type="MY_TYPE")) == 1

    def test_the_toolkit_name_is_tagged_onto_the_description(self, dtk, make_target_toolkit):
        target = make_target_toolkit("DOC_PROJECT_3")
        dtk._DocumentHandler(target, "Measurements", _documentSection("/m.parquet"),
                             False, "Measurements", "")
        assert target.getMeasurementsDocuments(type="MY_TYPE")[0].desc["toolkit"] == TARGET_TOOLKIT_NAME

    def test_a_duplicate_without_overwrite_is_not_added_again(self, dtk, make_target_toolkit):
        target = make_target_toolkit("DOC_PROJECT_4")
        dtk._DocumentHandler(target, "Measurements", _documentSection("/one.parquet"),
                             False, "Measurements", "")
        dtk._DocumentHandler(target, "Measurements", _documentSection("/two.parquet"),
                             False, "Measurements", "")
        documents = target.getMeasurementsDocuments(type="MY_TYPE")
        assert [d.resource for d in documents] == ["/one.parquet"]

    def test_overwrite_updates_the_resource_of_the_existing_document(self, dtk, make_target_toolkit):
        target = make_target_toolkit("DOC_PROJECT_5")
        dtk._DocumentHandler(target, "Measurements", _documentSection("/one.parquet"),
                             False, "Measurements", "")
        dtk._DocumentHandler(target, "Measurements", _documentSection("/two.parquet"),
                             True, "Measurements", "")
        documents = target.getMeasurementsDocuments(type="MY_TYPE")
        assert [d.resource for d in documents] == ["/two.parquet"]

    def test_documents_differing_in_description_are_stored_side_by_side(self, dtk, make_target_toolkit):
        """The existence query is built from the item's own desc, so a
        different desc is a different document rather than a duplicate."""
        target = make_target_toolkit("DOC_PROJECT_6")
        dtk._DocumentHandler(target, "Measurements",
                             _documentSection("/a.parquet", desc={"station": "A"}),
                             False, "Measurements", "")
        dtk._DocumentHandler(target, "Measurements",
                             _documentSection("/b.parquet", desc={"station": "B"}),
                             False, "Measurements", "")
        assert len(target.getMeasurementsDocuments(type="MY_TYPE")) == 2


@pytest.mark.unit
class TestDocumentHandlerRejectsUnknownDocumentTypes:
    """B155: the `if retrieveFunc is None` guard is unreachable."""

    @pytest.mark.xfail(
        strict=True,
        reason="B155: _DocumentHandler does "
               "`retrieveFunc = getattr(toolkit, f'get{documentType}Documents')` "
               "and only then checks `if retrieveFunc is None: raise ValueError`. "
               "getattr without a default raises AttributeError first, so the "
               "guard is dead code and the helpful ValueError naming the legal "
               "document types can never be produced. "
               "See the consolidated findings issue.",
    )
    def test_an_unknown_document_type_raises_the_documented_valueerror(self, dtk, make_target_toolkit):
        target = make_target_toolkit("DOC_TYPE_1")
        with pytest.raises(ValueError, match="must be"):
            dtk._DocumentHandler(target, "Bogus", _documentSection("/x.parquet"),
                                 False, "Bogus", "")

    def test_an_unknown_document_type_raises_attributeerror_instead(self, dtk, make_target_toolkit):
        """Characterisation of B155."""
        target = make_target_toolkit("DOC_TYPE_2")
        with pytest.raises(AttributeError, match="getBogusDocuments"):
            dtk._DocumentHandler(target, "Bogus", _documentSection("/x.parquet"),
                                 False, "Bogus", "")


@pytest.mark.unit
class TestCacheSectionIsUnreachable:
    """B153: the Cache lambda's parameter name does not match the call."""

    @pytest.mark.xfail(
        strict=True,
        reason="B153: in the handler table the Measurements/Simulations "
               "lambdas take (toolkit, itemName, docTypeDict, overwrite, "
               "basedir) but the Cache lambda takes (toolkit, itemName, "
               "itemDesc, overwrite, basedir), while the dispatcher always "
               "calls handler(toolkit=..., itemName=..., docTypeDict=..., "
               "overwrite=..., basedir=...). The Cache lambda therefore raises "
               "TypeError (unexpected keyword 'docTypeDict') on every call, and "
               "the surrounding `except Exception: logger.error` hides it -- a "
               "Cache section in a repository JSON can never be loaded. "
               "See the consolidated findings issue.",
    )
    def test_a_cache_section_adds_a_cache_document(self, dtk, make_target_toolkit):
        dtk.loadAllDatasourcesInRepositoryJSONToProject(
            projectName="CACHE_PROJECT_1",
            repositoryJSON={TARGET_TOOLKIT_NAME: {"Cache": _documentSection("/c.parquet")}},
        )
        assert len(make_target_toolkit("CACHE_PROJECT_1").getCacheDocuments(type="MY_TYPE")) == 1

    def test_a_cache_section_silently_loads_nothing(self, dtk, make_target_toolkit):
        """Characterisation of B153."""
        dtk.loadAllDatasourcesInRepositoryJSONToProject(
            projectName="CACHE_PROJECT_2",
            repositoryJSON={TARGET_TOOLKIT_NAME: {"Cache": _documentSection("/c.parquet")}},
        )
        assert len(make_target_toolkit("CACHE_PROJECT_2").getCacheDocuments(type="MY_TYPE")) == 0

    def test_the_cache_handler_rejects_the_keyword_the_dispatcher_uses(self, dtk, make_target_toolkit):
        """Characterisation of B153: the same failure, one layer down,
        with the swallowing try/except taken out of the picture."""
        target = make_target_toolkit("CACHE_PROJECT_3")

        # rebuild exactly the lambda the module installs for 'Cache'
        cacheHandler = (lambda toolkit, itemName, itemDesc, overwrite, basedir:
                        dtk._DocumentHandler(toolkit, itemName, itemDesc, overwrite, "Cache", basedir))
        with pytest.raises(TypeError, match="docTypeDict"):
            cacheHandler(toolkit=target, itemName="Cache",
                         docTypeDict=_documentSection("/c.parquet"),
                         overwrite=False, basedir="")


# ---------------------------------------------------------------------------
# _handle_Function
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestHandleFunction:
    def test_a_dict_value_is_one_call_with_overwrite_appended(self, dtk, tmp_path):
        path = _writeRepository(tmp_path, "inner", {})
        dtk._handle_Function(
            None, "Function",
            {"addRepository": {"repositoryName": "fromFunction", "repositoryPath": path}},
            False, "",
        )
        assert dtk.getDataSourceList() == ["fromFunction"]

    def test_a_list_value_is_one_call_per_element(self, dtk, tmp_path):
        first = _writeRepository(tmp_path, "first", {})
        second = _writeRepository(tmp_path, "second", {})
        dtk._handle_Function(
            None, "Function",
            {"addRepository": [
                {"repositoryName": "one", "repositoryPath": first},
                {"repositoryName": "two", "repositoryPath": second},
            ]},
            False, "",
        )
        assert sorted(dtk.getDataSourceList()) == ["one", "two"]

    def test_a_non_dict_element_in_the_list_is_skipped(self, dtk, tmp_path):
        path = _writeRepository(tmp_path, "inner", {})
        dtk._handle_Function(
            None, "Function",
            {"addRepository": ["nonsense", {"repositoryName": "ok", "repositoryPath": path}]},
            False, "",
        )
        assert dtk.getDataSourceList() == ["ok"]

    def test_a_scalar_value_raises(self, dtk):
        with pytest.raises(ValueError, match="must be dict"):
            dtk._handle_Function(None, "Function", {"addRepository": "nonsense"}, False, "")

    def test_an_unknown_function_name_raises_attributeerror(self, dtk):
        with pytest.raises(AttributeError):
            dtk._handle_Function(None, "Function", {"noSuchMethod": {}}, False, "")


@pytest.mark.unit
class TestHandleFunctionTargetsTheWrongObject:
    """B156: methods are resolved on the dataToolkit, not on the toolkit."""

    @pytest.mark.xfail(
        strict=True,
        reason="B156: _handle_Function does `getattr(self, itemName)` where "
               "self is the dataToolkit (bound to the read-only default "
               "project), not `getattr(toolkit, itemName)`. A Function section "
               "sits under a toolkit key in the repository JSON, so its calls "
               "are meant for that toolkit in the target project; instead the "
               "toolkit argument is unused and only dataToolkit's own methods "
               "can ever be invoked -- and any of them that touch the config "
               "blow up because the default project forbids configuration. "
               "See the consolidated findings issue.",
    )
    def test_a_function_section_calls_the_method_on_the_target_toolkit(self, dtk, make_target_toolkit):
        dtk.loadAllDatasourcesInRepositoryJSONToProject(
            projectName="FUNCTION_PROJECT_1",
            repositoryJSON={TARGET_TOOLKIT_NAME: {"Function": {"setConfig": {"answer": 42}}}},
        )
        assert make_target_toolkit("FUNCTION_PROJECT_1").getConfig()["answer"] == 42

    def test_a_function_section_calls_the_method_on_the_datatoolkit(self, dtk, make_target_toolkit):
        """Characterisation of B156."""
        dtk.loadAllDatasourcesInRepositoryJSONToProject(
            projectName="FUNCTION_PROJECT_2",
            repositoryJSON={TARGET_TOOLKIT_NAME: {"Function": {"setConfig": {"answer": 42}}}},
        )
        # nothing reached the target project ...
        assert "answer" not in make_target_toolkit("FUNCTION_PROJECT_2").getConfig()
        # ... because setConfig ran on the default project, which refuses it
        with pytest.raises(ValueError, match="[Dd]efault project"):
            dtk.setConfig(answer=42)


# ---------------------------------------------------------------------------
# _makeItemPathAbsolute
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestMakeItemPathAbsolute:
    def test_a_relative_resource_is_joined_onto_the_basedir_by_default(self, dtk):
        assert dtk._makeItemPathAbsolute({"resource": "a/b.parquet"}, "/base") == "/base/a/b.parquet"

    def test_a_boolean_false_flag_leaves_the_resource_untouched(self, dtk):
        item = {"resource": "a/b.parquet", "isRelativePath": False}
        assert dtk._makeItemPathAbsolute(item, "/base") == "a/b.parquet"

    def test_an_absolute_resource_survives_the_join(self, dtk):
        item = {"resource": "/already/absolute", "isRelativePath": True}
        assert dtk._makeItemPathAbsolute(item, "/base") == "/already/absolute"


@pytest.mark.unit
class TestMakeItemPathAbsoluteMisreadsTheStringFlag:
    """B154: isRelativePath is coerced with bool(), so "False" is truthy."""

    @pytest.mark.xfail(
        strict=True,
        reason='B154: _makeItemPathAbsolute does '
               'bool(theItem.get("isRelativePath", True)). Repository JSON '
               'carries this flag as the STRING "True"/"False" -- the sibling '
               '_handle_DataSource asserts exactly that -- and bool("False") '
               'is True, so a resource explicitly marked as not relative is '
               'still joined onto basedir. Affects every Measurements / '
               'Simulations / Cache item, unlike _handle_DataSource which '
               'compares against the strings correctly. '
               "See the consolidated findings issue.",
    )
    def test_the_string_false_leaves_the_resource_untouched(self, dtk):
        item = {"resource": "a/b.parquet", "isRelativePath": "False"}
        assert dtk._makeItemPathAbsolute(item, "/base") == "a/b.parquet"

    def test_the_string_false_is_treated_as_relative(self, dtk):
        """Characterisation of B154."""
        item = {"resource": "a/b.parquet", "isRelativePath": "False"}
        assert dtk._makeItemPathAbsolute(item, "/base") == "/base/a/b.parquet"

    def test_only_a_direct_call_can_reach_the_misreading(self, dtk):
        """B154 is only observable through a direct call, because the one
        caller of _makeItemPathAbsolute cannot pass the flag at all -- see
        B159."""
        assert dtk._makeItemPathAbsolute({"resource": "x", "isRelativePath": "False"}, "/base") \
            != dtk._makeItemPathAbsolute({"resource": "x", "isRelativePath": False}, "/base")


@pytest.mark.unit
class TestDocumentHandlerForwardsTheControlFlagAsData:
    """B159: the isRelativePath key is never removed from the item."""

    def _section(self, isRelativePath):
        return {"D": {"item": {"resource": "a/b.parquet", "dataFormat": "parquet",
                               "type": "MY_TYPE", "isRelativePath": isRelativePath}}}

    @pytest.mark.xfail(
        strict=True,
        reason="B159: _DocumentHandler reads isRelativePath out of the item "
               "(via _makeItemPathAbsolute) but never pops it, then calls "
               "add<Type>Document(**theItem). add<Type>Document takes only "
               "resource/dataFormat/type/desc, so any Measurements, "
               "Simulations or Cache entry that declares isRelativePath dies "
               "with TypeError -- silently, since the loader logs and "
               "continues. Note the two handlers also disagree about where "
               "the flag lives: _handle_DataSource reads it from the wrapper "
               "next to 'item', _DocumentHandler from inside 'item'. "
               "See the consolidated findings issue.",
    )
    def test_an_item_may_declare_isrelativepath(self, dtk, make_target_toolkit, tmp_path):
        target = make_target_toolkit("REL_FLAG_PROJECT_1")
        dtk._DocumentHandler(target, "Measurements", self._section("True"),
                             False, "Measurements", str(tmp_path))
        assert len(target.getMeasurementsDocuments(type="MY_TYPE")) == 1

    def test_declaring_isrelativepath_raises_a_typeerror(self, dtk, make_target_toolkit, tmp_path):
        """Characterisation of B159."""
        target = make_target_toolkit("REL_FLAG_PROJECT_2")
        with pytest.raises(TypeError, match="isRelativePath"):
            dtk._DocumentHandler(target, "Measurements", self._section("False"),
                                 False, "Measurements", str(tmp_path))

    def test_through_the_loader_such_an_item_silently_vanishes(self, dtk, make_target_toolkit, tmp_path):
        """Characterisation of B159."""
        dtk.loadAllDatasourcesInRepositoryJSONToProject(
            projectName="REL_FLAG_PROJECT_3",
            repositoryJSON={TARGET_TOOLKIT_NAME: {"Measurements": self._section("True")}},
            basedir=str(tmp_path),
        )
        target = make_target_toolkit("REL_FLAG_PROJECT_3")
        assert len(target.getMeasurementsDocuments(type="MY_TYPE")) == 0


# ---------------------------------------------------------------------------
# repository -> project, via the DB
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestLoadAllDatasourcesInRepositoryToProject:
    def test_relative_resources_resolve_against_the_repository_file_directory(
        self, dtk, make_target_toolkit, tmp_path
    ):
        repositoryDirectory = tmp_path / "repositories"
        repositoryDirectory.mkdir()
        path = _writeRepository(repositoryDirectory, "repo", {TARGET_TOOLKIT_NAME: {
            "Datasource": {"A": _datasourceSection("data/x.parquet", isRelativePath="True",
                                                   dataFormat="parquet")}}})
        dtk.addRepository("myRepo", path)
        dtk.loadAllDatasourcesInRepositoryToProject(projectName="REPO_PROJECT_1",
                                                    repositoryName="myRepo")
        document = make_target_toolkit("REPO_PROJECT_1").getDataSourceDocument("A")
        assert document.resource == str(repositoryDirectory / "data/x.parquet")

    def test_an_unregistered_repository_name_raises(self, dtk):
        with pytest.raises(AttributeError):
            dtk.loadAllDatasourcesInRepositoryToProject(projectName="P",
                                                        repositoryName="noSuchRepository")


@pytest.mark.unit
class TestLoadAllDatasourcesInAllRepositoriesToProject:
    def test_every_registered_repository_is_loaded(self, dtk, make_target_toolkit, tmp_path):
        dtk.addRepository("one", _writeRepository(tmp_path, "one", {TARGET_TOOLKIT_NAME: {
            "Datasource": {"A": _datasourceSection("/a")}}}))
        dtk.addRepository("two", _writeRepository(tmp_path, "two", {TARGET_TOOLKIT_NAME: {
            "Datasource": {"B": _datasourceSection("/b")}}}))
        dtk.loadAllDatasourcesInAllRepositoriesToProject(projectName="ALL_PROJECT_1")
        names = make_target_toolkit("ALL_PROJECT_1").getDataSourceList()
        assert sorted(names) == ["A", "B"]

    def test_a_repository_that_raises_valueerror_does_not_stop_the_others(
        self, dtk, make_target_toolkit, tmp_path
    ):
        dtk.addRepository("broken", _writeRepository(tmp_path, "broken", {TARGET_TOOLKIT_NAME: {
            "Bogus": {}}}))
        dtk.addRepository("good", _writeRepository(tmp_path, "good", {TARGET_TOOLKIT_NAME: {
            "Config": {"loaded": True}}}))
        dtk.loadAllDatasourcesInAllRepositoriesToProject(projectName="ALL_PROJECT_2")
        assert make_target_toolkit("ALL_PROJECT_2").getConfig()["loaded"] is True

    def test_an_empty_registry_is_a_no_op(self, dtk):
        assert dtk.loadAllDatasourcesInAllRepositoriesToProject(projectName="P") is None


@pytest.mark.unit
class TestRealToolkitHomeWiring:
    """One end-to-end pass through the production ToolkitHome, to prove the
    stand-in used elsewhere in this module is standing in for something real.

    ToolkitHome is subclassed rather than patched: the subclass only injects
    filesDirectory so the target toolkit writes under tmp_path instead of
    ~/.hera/<projectName>. Nothing on the real class is mutated.
    """

    @pytest.fixture()
    def dtkWithRealHome(self, monkeypatch, unit_files_directory):
        from hera.toolkit import ToolkitHome

        temporaryFilesDirectory = unit_files_directory

        class _TempFilesToolkitHome(ToolkitHome):
            def __init__(self, projectName=None, filesDirectory=None):
                super().__init__(projectName=projectName,
                                 filesDirectory=filesDirectory or temporaryFilesDirectory)

            def getToolkit(self, toolkitName, filesDirectory=None, **kwargs):
                return super().getToolkit(
                    toolkitName,
                    filesDirectory=filesDirectory or temporaryFilesDirectory,
                    **kwargs,
                )

        monkeypatch.setattr(dataToolkitModule, "ToolkitHome", _TempFilesToolkitHome)
        return dataToolkit()

    def test_a_config_section_reaches_a_real_static_toolkit(self, dtkWithRealHome, unit_files_directory):
        from hera import toolkitHome

        dtkWithRealHome.loadAllDatasourcesInRepositoryJSONToProject(
            projectName="REAL_HOME_PROJECT",
            repositoryJSON={toolkitHome.METEOROLOGY_LOWFREQ: {"Config": {"stationType": "IMS"}}},
        )
        target = toolkitHome.getToolkit(
            toolkitHome.METEOROLOGY_LOWFREQ,
            projectName="REAL_HOME_PROJECT",
            filesDirectory=unit_files_directory,
        )
        assert target.getConfig()["stationType"] == "IMS"


# ---------------------------------------------------------------------------
# the DB-free helpers
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestResolveDataSourcePaths:
    def test_a_relative_resource_under_an_item_wrapper_becomes_absolute(self):
        source = {"T": {"Datasource": {"A": {"isRelativePath": "True",
                                             "item": {"resource": "x.parquet"}}}}}
        resolved = dataToolkit.resolveDataSourcePaths(source, basedir="/base")
        assert resolved["T"]["Datasource"]["A"]["item"]["resource"] == "/base/x.parquet"

    def test_the_flag_may_live_on_the_item_itself(self):
        source = {"T": {"S": {"A": {"item": {"resource": "x", "isRelativePath": True}}}}}
        resolved = dataToolkit.resolveDataSourcePaths(source, basedir="/base")
        assert resolved["T"]["S"]["A"]["item"]["resource"] == "/base/x"

    def test_an_entry_with_no_item_wrapper_is_resolved_in_place(self):
        source = {"T": {"S": {"A": {"resource": "x", "isRelativePath": "True"}}}}
        resolved = dataToolkit.resolveDataSourcePaths(source, basedir="/base")
        assert resolved["T"]["S"]["A"]["resource"] == "/base/x"

    def test_a_resource_not_marked_relative_is_left_alone(self):
        source = {"T": {"S": {"A": {"isRelativePath": "False", "item": {"resource": "x"}}}}}
        resolved = dataToolkit.resolveDataSourcePaths(source, basedir="/base")
        assert resolved["T"]["S"]["A"]["item"]["resource"] == "x"

    def test_non_dict_values_at_every_level_are_skipped(self):
        source = {"toolkitIsAString": "nope",
                  "sectionIsAString": {"S": "nope"},
                  "itemIsAString": {"S": {"A": "nope"}},
                  "itemHasNoResource": {"S": {"A": {"item": {"dataFormat": "parquet"}}}}}
        assert dataToolkit.resolveDataSourcePaths(source, basedir="/base") == source

    def test_the_input_is_deep_copied_not_mutated(self):
        source = {"T": {"S": {"A": {"isRelativePath": "True", "item": {"resource": "x"}}}}}
        dataToolkit.resolveDataSourcePaths(source, basedir="/base")
        assert source["T"]["S"]["A"]["item"]["resource"] == "x"


@pytest.mark.unit
class TestLoadRepositoryFromPath:
    def test_it_resolves_relative_resources_against_the_files_own_directory(self, tmp_path):
        path = _writeRepository(tmp_path, "repo", {"T": {"Datasource": {
            "A": {"isRelativePath": "True", "item": {"resource": "x.parquet"}}}}})
        loaded = dataToolkit.loadRepositoryFromPath(path)
        assert loaded["T"]["Datasource"]["A"]["item"]["resource"] == str(tmp_path / "x.parquet")

    def test_a_missing_file_raises_filenotfounderror(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="Repository JSON not found"):
            dataToolkit.loadRepositoryFromPath(str(tmp_path / "absent.json"))

    def test_a_directory_is_reported_as_not_a_file(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="Repository JSON not found"):
            dataToolkit.loadRepositoryFromPath(str(tmp_path))

    def test_it_is_reachable_without_constructing_the_toolkit(self, tmp_path):
        """Both direct-load helpers are staticmethods on purpose -- no DB."""
        path = _writeRepository(tmp_path, "repo", {})
        assert dataToolkit.loadRepositoryFromPath(path) == {}


# ---------------------------------------------------------------------------
# export
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestResolveRepositoryPath:
    def test_a_json_path_is_made_absolute(self, dtk, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        assert dtk._resolveRepositoryPath("sub/repo.json") == str(tmp_path / "sub/repo.json")

    def test_a_registered_name_resolves_to_its_stored_resource(self, dtk, tmp_path):
        path = _writeRepository(tmp_path, "repo", {})
        dtk.addRepository("myRepo", path)
        assert dtk._resolveRepositoryPath("myRepo") == path

    def test_an_unknown_bare_name_lands_in_the_current_directory(self, dtk, tmp_path, monkeypatch):
        """Documented ("create alongside cwd") but worth stating: an
        unregistered bare name makes export write into whatever directory the
        process happens to be in."""
        monkeypatch.chdir(tmp_path)
        assert dtk._resolveRepositoryPath("unknownName") == str(tmp_path / "unknownName.json")


@pytest.mark.unit
class TestResolveDocumentsForExport:
    def test_all_documents_of_the_project_are_returned_by_default(self, dtk, patched_datalayer_project):
        project = patched_datalayer_project(projectName="EXPORT_PROJECT_1")
        project.addMeasurementsDocument(resource="/a", dataFormat="parquet", type="T", desc={})
        project.addSimulationsDocument(resource="/b", dataFormat="parquet", type="T", desc={})
        assert len(dtk._resolveDocumentsForExport(project, None)) == 2

    def test_the_projects_own_config_document_is_excluded(self, dtk, patched_datalayer_project):
        project = patched_datalayer_project(projectName="EXPORT_PROJECT_2")
        project.setConfig(filesDirectory="/wherever")
        project.addMeasurementsDocument(resource="/a", dataFormat="parquet", type="T", desc={})
        types = [document.type for document in dtk._resolveDocumentsForExport(project, None)]
        assert types == ["T"]

    def test_a_single_document_is_wrapped_in_a_list(self, dtk, patched_datalayer_project):
        project = patched_datalayer_project(projectName="EXPORT_PROJECT_3")
        document = project.addMeasurementsDocument(resource="/a", dataFormat="parquet",
                                                   type="T", desc={})
        assert dtk._resolveDocumentsForExport(project, document) == [document]

    def test_a_document_id_string_is_looked_up(self, dtk, patched_datalayer_project):
        project = patched_datalayer_project(projectName="EXPORT_PROJECT_4")
        document = project.addMeasurementsDocument(resource="/a", dataFormat="parquet",
                                                   type="T", desc={})
        resolved = dtk._resolveDocumentsForExport(project, str(document.id))
        assert [d.resource for d in resolved] == ["/a"]


@pytest.mark.unit
class TestResolveDocumentsForExportRejectsUnknownIds:
    """B158: the `if doc is None` guard after getDocumentByID is unreachable."""

    @pytest.mark.xfail(
        strict=True,
        reason="B158: _resolveDocumentsForExport does "
               "`doc = proj.getDocumentByID(d)` then "
               "`if doc is None: raise ValueError('Document id not found ...')`. "
               "getDocumentByID goes through mongoengine's QuerySet.get(), "
               "which raises DoesNotExist rather than returning None, so the "
               "guard is dead code and callers get a mongoengine exception "
               "instead of the intended ValueError. "
               "See the consolidated findings issue.",
    )
    def test_an_unknown_document_id_raises_the_documented_valueerror(self, dtk, patched_datalayer_project):
        project = patched_datalayer_project(projectName="EXPORT_PROJECT_5")
        with pytest.raises(ValueError, match="Document id not found"):
            dtk._resolveDocumentsForExport(project, "0" * 24)

    def test_an_unknown_document_id_raises_doesnotexist(self, dtk, patched_datalayer_project):
        """Characterisation of B158."""
        from mongoengine.errors import DoesNotExist

        project = patched_datalayer_project(projectName="EXPORT_PROJECT_6")
        with pytest.raises(DoesNotExist):
            dtk._resolveDocumentsForExport(project, "0" * 24)


@pytest.mark.unit
class TestExportDocumentsToRepository:
    def test_it_writes_a_repository_json_keyed_by_the_toolkit_name(
        self, dtk, patched_datalayer_project, tmp_path
    ):
        project = patched_datalayer_project(projectName="EXPORT_RUN_1")
        project.addMeasurementsDocument(resource="/a.parquet", dataFormat="parquet",
                                        type="T1", desc={"station": "A"})
        destination = tmp_path / "exported.json"
        dtk.exportDocumentsToRepository(toolkitName="myToolkit",
                                        repositoryName=str(destination),
                                        projectName="EXPORT_RUN_1", register=False)
        with open(destination, encoding="utf-8") as handle:
            written = json.load(handle)
        assert list(written) == ["myToolkit"]
        assert list(written["myToolkit"]) == ["Measurements"]

    def test_the_report_names_the_added_documents(self, dtk, patched_datalayer_project, tmp_path):
        project = patched_datalayer_project(projectName="EXPORT_RUN_2")
        project.addMeasurementsDocument(resource="/a.parquet", dataFormat="parquet",
                                        type="T1", desc={})
        report = dtk.exportDocumentsToRepository(toolkitName="myToolkit",
                                                 repositoryName=str(tmp_path / "e.json"),
                                                 projectName="EXPORT_RUN_2", register=False)
        assert len(report["added"]) == 1
        assert report["skipped_existing"] == []

    def test_re_exporting_the_same_document_is_skipped_as_a_duplicate(
        self, dtk, patched_datalayer_project, tmp_path
    ):
        project = patched_datalayer_project(projectName="EXPORT_RUN_3")
        project.addMeasurementsDocument(resource="/a.parquet", dataFormat="parquet",
                                        type="T1", desc={})
        destination = str(tmp_path / "e.json")
        dtk.exportDocumentsToRepository(toolkitName="myToolkit", repositoryName=destination,
                                        projectName="EXPORT_RUN_3", register=False)
        report = dtk.exportDocumentsToRepository(toolkitName="myToolkit",
                                                 repositoryName=destination,
                                                 projectName="EXPORT_RUN_3", register=False)
        assert report["added"] == []
        assert len(report["skipped_existing"]) == 1

    def test_override_mode_adds_a_deduplication_count_to_the_report(
        self, dtk, patched_datalayer_project, tmp_path
    ):
        project = patched_datalayer_project(projectName="EXPORT_RUN_4")
        project.addMeasurementsDocument(resource="/a.parquet", dataFormat="parquet",
                                        type="T1", desc={})
        report = dtk.exportDocumentsToRepository(toolkitName="myToolkit",
                                                 repositoryName=str(tmp_path / "e.json"),
                                                 projectName="EXPORT_RUN_4", register=False,
                                                 mode="override")
        assert "deduplicated" in report

    def test_register_true_registers_the_written_file_as_a_repository(
        self, dtk, patched_datalayer_project, tmp_path
    ):
        project = patched_datalayer_project(projectName="EXPORT_RUN_5")
        project.addMeasurementsDocument(resource="/a.parquet", dataFormat="parquet",
                                        type="T1", desc={})
        destination = tmp_path / "myExport.json"
        dtk.exportDocumentsToRepository(toolkitName="myToolkit",
                                        repositoryName=str(destination),
                                        projectName="EXPORT_RUN_5", register=True)
        assert dtk.getDataSourceList() == ["myExport"]
        assert dtk.getDataSourceDocument("myExport").resource == str(destination)

    def test_an_existing_file_that_is_not_json_is_rejected(
        self, dtk, patched_datalayer_project, tmp_path
    ):
        patched_datalayer_project(projectName="EXPORT_RUN_6")
        destination = tmp_path / "broken.json"
        destination.write_text("this is not json", encoding="utf-8")
        with pytest.raises(ValueError, match="not valid JSON"):
            dtk.exportDocumentsToRepository(toolkitName="myToolkit",
                                            repositoryName=str(destination),
                                            projectName="EXPORT_RUN_6", register=False)

    def test_an_explicit_document_list_limits_what_is_exported(
        self, dtk, patched_datalayer_project, tmp_path
    ):
        project = patched_datalayer_project(projectName="EXPORT_RUN_7")
        first = project.addMeasurementsDocument(resource="/a.parquet", dataFormat="parquet",
                                                type="T1", desc={})
        project.addMeasurementsDocument(resource="/b.parquet", dataFormat="parquet",
                                        type="T1", desc={})
        report = dtk.exportDocumentsToRepository(toolkitName="myToolkit",
                                                 repositoryName=str(tmp_path / "e.json"),
                                                 projectName="EXPORT_RUN_7", register=False,
                                                 documents=[first])
        assert len(report["added"]) == 1
