"""The repository entrypoints: the registered list, one file, a folder, loading."""
import copy
import json
import logging
import os
import sys
import types

import pytest
from fastapi import HTTPException

import dependencies
from repository_api.errors import NoSuchRepository
from repository_api.load_log import LoadLog
from repository_api.repo_folder import RepoFolder
from repository_api.repo_json import RepositoryJson
from repository_api.repo_json_file import RepoJsonFile
from repository_api.repository_table import RepositoryTable
from repository_api.toolkit_check import ToolkitCheck
from routers import repositories as routes
from test_documents import FakeDatabase

# A repository tree with one relative datasource, as a repository file holds it.
REPO_TREE = {
    "MeteoLowFreq": {
        "Datasource": {
            "YAVNEEL": {
                "isRelativePath": "True",
                "item": {"resource": "meteo/YAVNEEL.parquet", "dataFormat": "parquet"},
            },
        },
    },
}


class FakeRepositories:
    """What the fake dataToolkit was asked to do, and the rows it answers."""

    def __init__(self):
        self.rows = []
        # Toolkit names FakeToolkitHome refuses to build.
        self.missing = set()
        # True when a load should add no document, as hera does when it skips.
        self.addsNothing = False
        # One entry per loaded tree: (projectName, tree, basedir, overwrite).
        self.jsonLoads = []
        # One entry per "every registered repository" load: (projectName, overwrite).
        self.allLoads = []


def _resolve(tree, basedir):
    """Enough of hera's own resolver for these tests: join the relative resources."""
    resolved = copy.deepcopy(tree)
    for toolkitValue in resolved.values():
        for section in toolkitValue.values():
            for item in section.values():
                if not isinstance(item, dict):
                    continue
                inner = item.get("item", item)
                if "resource" in inner and item.get("isRelativePath") == "True":
                    inner["resource"] = os.path.join(basedir, inner["resource"])
    return resolved


@pytest.fixture
def fake_hera(monkeypatch):
    """Install a fake hera: the database test_documents builds, and a dataToolkit."""
    db = FakeDatabase()
    recorder = FakeRepositories()

    class FakeToolkitHome:
        """Builds a toolkit, or refuses the names the test put in ``missing``."""

        def __init__(self, projectName=None, filesDirectory=None):
            self.projectName = projectName

        def getToolkit(self, toolkitName, **kwargs):
            if toolkitName in recorder.missing:
                raise ImportError(f"Cannot locate class: {toolkitName}")
            return types.SimpleNamespace(toolkitName=toolkitName)

    class FakeDataToolkit:
        def getRepositoryTable(self):
            import pandas as pd
            return pd.DataFrame(recorder.rows)

        def loadAllDatasourcesInRepositoryJSONToProject(self, projectName, repositoryJSON,
                                                        basedir="", overwrite=False):
            recorder.jsonLoads.append((projectName, repositoryJSON, basedir, overwrite))
            if recorder.addsNothing:
                return
            # A real load adds documents, so the project in the answer must change.
            db.add(f"loaded{len(recorder.jsonLoads)}", projectName, type="datasource")

        def loadAllDatasourcesInAllRepositoriesToProject(self, projectName, overwrite=False):
            recorder.allLoads.append((projectName, overwrite))
            db.add(f"all{len(recorder.allLoads)}", projectName, type="datasource")

        @staticmethod
        def resolveDataSourcePaths(tree, basedir=""):
            return _resolve(tree, basedir)

    toolkitModule = types.ModuleType("hera.toolkit")
    toolkitModule.ToolkitHome = FakeToolkitHome

    datalayer = types.ModuleType("hera.datalayer")
    datalayer.All = types.SimpleNamespace(getDocumentsAsDict=db.getDocumentsAsDict)
    toolkit = types.ModuleType("hera.utils.data.toolkit")
    toolkit.dataToolkit = FakeDataToolkit

    hera = types.ModuleType("hera")
    hera.datalayer = datalayer
    hera.toolkit = toolkitModule
    monkeypatch.setitem(sys.modules, "hera", hera)
    monkeypatch.setitem(sys.modules, "hera.toolkit", toolkitModule)
    monkeypatch.setitem(sys.modules, "hera.datalayer", datalayer)
    monkeypatch.setitem(sys.modules, "hera.utils", types.ModuleType("hera.utils"))
    monkeypatch.setitem(sys.modules, "hera.utils.data", types.ModuleType("hera.utils.data"))
    monkeypatch.setitem(sys.modules, "hera.utils.data.toolkit", toolkit)
    recorder.db = db
    return recorder


@pytest.fixture
def warmed(monkeypatch):
    monkeypatch.setattr(dependencies.warmup, "_ready", True)


def _write(path, content):
    """Write a repository file, making its folder on the way."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(content, str):
        path.write_text(content)
    else:
        path.write_text(json.dumps(content))
    return str(path)


def _payload(**fields):
    payload = {
        "projectName": "Alpha",
        "path": "",
        "repositoryJSON": None,
        "baseDir": "",
        "overwrite": False,
    }
    payload.update(fields)
    return types.SimpleNamespace(**payload)


class TestRepositoryJson:
    def test_a_tree_of_toolkits_and_sections_is_a_repository(self):
        assert RepositoryJson.is_repository(REPO_TREE)

    def test_an_empty_tree_is_not_a_repository(self):
        assert not RepositoryJson.is_repository({})

    def test_a_list_is_not_a_repository(self):
        assert not RepositoryJson.is_repository([1, 2, 3])

    def test_a_toolkit_that_is_not_a_dict_is_not_a_repository(self):
        assert not RepositoryJson.is_repository({"MyToolkit": "not a dict"})

    def test_an_unknown_section_is_not_a_repository(self):
        assert not RepositoryJson.is_repository({"MyToolkit": {"config": {}, "nowhere": {}}})

    def test_says_why_a_tree_is_not_a_repository(self):
        why = RepositoryJson.why_not_a_repository
        assert why([1, 2]) == "the top level is list, not an object of toolkit names"
        assert why({}) == "it is empty"
        assert why({"T": "x"}) == "toolkit T is str, not an object"
        assert "T.nowhere is not a section" in why({"T": {"nowhere": {}}})
        assert why(REPO_TREE) == ""

    def test_the_section_names_may_be_in_any_case(self):
        assert RepositoryJson.is_repository({"MyToolkit": {"Config": {}, "Measurements": {}}})

    def test_loads_the_tree_against_its_base_directory(self, fake_hera):
        RepositoryJson(REPO_TREE, "/repos/meteo").load_into("Alpha", overwrite=True)

        assert fake_hera.jsonLoads == [("Alpha", REPO_TREE, "/repos/meteo", True)]

    def test_a_relative_base_directory_is_made_absolute(self, fake_hera):
        RepositoryJson(REPO_TREE, "repos").load_into("Alpha")

        assert fake_hera.jsonLoads[0][2] == os.path.join(os.getcwd(), "repos")


class TestRepoJsonFile:
    def test_reads_the_tree_of_a_repository_file(self, tmp_path):
        path = _write(tmp_path / "repo.json", REPO_TREE)

        assert RepoJsonFile(path).tree() == REPO_TREE

    def test_a_path_that_is_not_there_says_what_the_system_said(self, tmp_path):
        with pytest.raises(NoSuchRepository) as raised:
            RepoJsonFile(str(tmp_path / "nowhere.json")).tree()

        assert "No such file or directory" in str(raised.value)

    def test_a_file_that_is_not_json_says_where_the_json_broke(self, tmp_path):
        path = _write(tmp_path / "broken.json", "this is not json {{{")

        with pytest.raises(NoSuchRepository) as raised:
            RepoJsonFile(path).tree()

        assert "is not JSON" in str(raised.value)
        assert "line 1" in str(raised.value)

    def test_a_json_file_that_is_not_a_repository_says_which_key_is_wrong(self, tmp_path):
        path = _write(tmp_path / "other.json", {"MyToolkit": {"nowhere": {}}})

        with pytest.raises(NoSuchRepository) as raised:
            RepoJsonFile(path).tree()

        assert "MyToolkit.nowhere is not a section" in str(raised.value)

    def test_the_refusal_keeps_the_error_it_came_from(self, tmp_path):
        with pytest.raises(NoSuchRepository) as raised:
            RepoJsonFile(str(tmp_path / "nowhere.json")).tree()

        assert isinstance(raised.value.__cause__, OSError)

    def test_loading_a_file_resolves_against_its_own_folder(self, fake_hera, tmp_path):
        path = _write(tmp_path / "meteo" / "repo.json", REPO_TREE)

        RepoJsonFile(path).load_into("Alpha")

        assert fake_hera.jsonLoads[0][2] == str(tmp_path / "meteo")


class TestRepoFolder:
    def test_finds_the_repository_files_and_skips_the_rest(self, tmp_path):
        _write(tmp_path / "valid.json", REPO_TREE)
        _write(tmp_path / "sub" / "nested.json", REPO_TREE)
        _write(tmp_path / "broken.json", "not json {{{")
        _write(tmp_path / "badsection.json", {"MyToolkit": {"nowhere": {}}})
        _write(tmp_path / "readme.txt", "not a json file")

        files = RepoFolder(str(tmp_path)).files()

        assert files == [str(tmp_path / "sub" / "nested.json"), str(tmp_path / "valid.json")]

    def test_skips_a_case_configuration(self, tmp_path):
        _write(tmp_path / "caseConfiguration.json", REPO_TREE)

        assert RepoFolder(str(tmp_path)).files() == []

    def test_a_folder_that_is_not_there_is_empty(self, tmp_path):
        assert RepoFolder(str(tmp_path / "nowhere")).files() == []
        assert RepoFolder(str(tmp_path / "nowhere")).repositories() == {}

    def test_skips_a_file_it_cannot_read_and_keeps_the_rest(self, tmp_path):
        good = _write(tmp_path / "valid.json", REPO_TREE)
        _write(tmp_path / "broken.json", "not json {{{")

        assert RepoFolder(str(tmp_path)).files() == [good]

    def test_the_contents_come_with_the_files_when_asked(self, tmp_path):
        path = _write(tmp_path / "valid.json", REPO_TREE)

        assert RepoFolder(str(tmp_path)).repositories() == {path: REPO_TREE}

    def test_expands_a_home_folder(self, tmp_path, monkeypatch):
        monkeypatch.setenv("HOME", str(tmp_path))
        path = _write(tmp_path / "repos" / "valid.json", REPO_TREE)

        assert RepoFolder("~/repos").files() == [path]


class TestRepositoryTable:
    def test_answers_the_registered_rows(self, fake_hera):
        fake_hera.rows = [{"datasourceName": "YAVNEEL", "resource": "/repos/meteo.json"}]

        assert RepositoryTable().rows() == fake_hera.rows

    def test_an_empty_table_is_no_rows(self, fake_hera):
        assert RepositoryTable().rows() == []

    def test_loads_every_registered_repository(self, fake_hera):
        RepositoryTable().load_all_into("Alpha", overwrite=True)

        assert fake_hera.allLoads == [("Alpha", True)]


class TestRegisteredRoute:
    def test_answers_the_rows(self, fake_hera, warmed):
        fake_hera.rows = [{"datasourceName": "YAVNEEL", "resource": "/repos/meteo.json"}]

        assert routes.registered() == {"repositories": fake_hera.rows}


class TestGetRoute:
    def test_answers_the_file_and_its_resolved_copy(self, fake_hera, warmed, tmp_path):
        path = _write(tmp_path / "meteo" / "repo.json", REPO_TREE)

        answer = routes.get(path)

        assert answer["json"] == REPO_TREE
        resolved = answer["resolved"]["MeteoLowFreq"]["Datasource"]["YAVNEEL"]["item"]
        assert resolved["resource"] == str(tmp_path / "meteo" / "meteo" / "YAVNEEL.parquet")

    def test_a_path_that_is_not_there_is_not_found(self, fake_hera, warmed, tmp_path):
        with pytest.raises(HTTPException) as raised:
            routes.get(str(tmp_path / "nowhere.json"))

        assert raised.value.status_code == 404

    def test_a_file_that_is_not_json_is_not_found(self, fake_hera, warmed, tmp_path):
        path = _write(tmp_path / "broken.json", "not json {{{")

        with pytest.raises(HTTPException) as raised:
            routes.get(path)

        assert raised.value.status_code == 404

    def test_a_json_file_that_is_not_a_repository_is_not_found(self, fake_hera, warmed, tmp_path):
        path = _write(tmp_path / "other.json", {"MyToolkit": {"nowhere": {}}})

        with pytest.raises(HTTPException) as raised:
            routes.get(path)

        assert raised.value.status_code == 404


class TestJsonsInFolderRoute:
    def test_answers_the_names_only_by_default(self, fake_hera, warmed, tmp_path):
        path = _write(tmp_path / "valid.json", REPO_TREE)

        answer = routes.jsonsInFolder(str(tmp_path))

        assert answer == {"files": [path]}

    def test_answers_the_contents_when_asked(self, fake_hera, warmed, tmp_path):
        path = _write(tmp_path / "valid.json", REPO_TREE)

        answer = routes.jsonsInFolder(str(tmp_path), withContents=True)

        assert answer == {"files": [path], "repoJsons": {path: REPO_TREE}}

    def test_a_folder_that_is_not_there_is_empty_and_not_an_error(self, fake_hera, warmed, tmp_path):
        answer = routes.jsonsInFolder(str(tmp_path / "nowhere"), withContents=True)

        assert answer == {"files": [], "repoJsons": {}}


class TestLoadIntoProjectRoute:
    def test_a_path_loads_the_file_against_its_own_folder(self, fake_hera, warmed, tmp_path):
        path = _write(tmp_path / "meteo" / "repo.json", REPO_TREE)

        answer = routes.loadIntoProject(_payload(path=path))

        assert fake_hera.jsonLoads == [("Alpha", REPO_TREE, str(tmp_path / "meteo"), False)]
        assert answer["project"]["name"] == "Alpha"
        assert len(answer["project"]["documents"]) == 1

    def test_a_tree_loads_against_the_base_directory_it_is_given(self, fake_hera, warmed):
        routes.loadIntoProject(_payload(
            repositoryJSON=REPO_TREE, baseDir="/projects/alpha", overwrite=True))

        assert fake_hera.jsonLoads == [("Alpha", REPO_TREE, "/projects/alpha", True)]

    def test_no_source_loads_every_registered_repository(self, fake_hera, warmed):
        routes.loadIntoProject(_payload())

        assert fake_hera.allLoads == [("Alpha", False)]
        assert fake_hera.jsonLoads == []

    def test_a_path_and_a_tree_at_once_is_refused(self, fake_hera, warmed, tmp_path):
        path = _write(tmp_path / "repo.json", REPO_TREE)

        with pytest.raises(HTTPException) as raised:
            routes.loadIntoProject(_payload(path=path, repositoryJSON=REPO_TREE))

        assert raised.value.status_code == 400
        assert fake_hera.jsonLoads == []

    def test_a_path_that_is_not_there_is_not_found(self, fake_hera, warmed, tmp_path):
        with pytest.raises(HTTPException) as raised:
            routes.loadIntoProject(_payload(path=str(tmp_path / "nowhere.json")))

        assert raised.value.status_code == 404

    def test_a_blank_project_name_is_refused(self, fake_hera, warmed):
        with pytest.raises(HTTPException) as raised:
            routes.loadIntoProject(_payload(projectName="  "))

        assert raised.value.status_code == 400


class TestLoadLog:
    def test_one_complaint_is_reported_once_not_once_per_logger(self):
        # The handler sits on hera, hera.utils, hera.utils.data and so on.
        logging.getLogger("hera.utils.data.toolkit")

        with LoadLog() as log:
            logging.getLogger("hera.utils.data.toolkit").error("said once")

        assert log.problems == ["hera.utils.data.toolkit: said once"]

    def test_keeps_the_errors_hera_logs_instead_of_raising(self):
        with LoadLog() as log:
            logging.getLogger("hera.utils.data.toolkit").error("adding *Datasource* failed")

        assert log.problems == ["hera.utils.data.toolkit: adding *Datasource* failed"]

    def test_keeps_the_line_hera_skips_a_toolkit_with(self):
        with LoadLog() as log:
            logging.getLogger("hera.utils.data.toolkit").info("Skipping key 'GIS_Tiles'")

        assert log.problems == ["hera.utils.data.toolkit: Skipping key 'GIS_Tiles'"]

    def test_keeps_no_ordinary_chatter(self):
        with LoadLog() as log:
            logging.getLogger("hera.utils.data.toolkit").info("Loading repository")

        assert log.problems == []

    def test_hears_a_family_hera_mutes_at_critical_with_no_propagation(self):
        muted = logging.getLogger("hera.measurements")
        muted.setLevel(logging.CRITICAL)
        muted.propagate = False

        with LoadLog() as log:
            logging.getLogger("hera.measurements.GIS.toolkit").error("addDataSource failed")

        assert log.problems == ["hera.measurements.GIS.toolkit: addDataSource failed"]
        assert muted.level == logging.CRITICAL

    def test_keeps_the_traceback_when_the_log_line_carries_one(self):
        with LoadLog() as log:
            try:
                raise ValueError("the real reason")
            except ValueError:
                logging.getLogger("hera.utils.data.toolkit").error("failed", exc_info=True)

        assert "the real reason" in "\n".join(log.problems)
        assert "ValueError" in "\n".join(log.problems)

    def test_names_the_logger_the_complaint_came_from(self):
        with LoadLog() as log:
            logging.getLogger("hera.riskassessment.riskToolkit").error("no agent")

        assert log.problems == ["hera.riskassessment.riskToolkit: no agent"]

    def test_lets_go_of_the_logger_afterwards(self):
        with LoadLog() as log:
            pass
        logging.getLogger("hera.utils.data.toolkit").error("after the load")

        assert log.problems == []
        assert log not in logging.getLogger("hera").handlers


class TestToolkitCheck:
    def test_no_problem_when_every_toolkit_builds(self, fake_hera):
        assert ToolkitCheck(REPO_TREE, "Alpha").problems() == []

    def test_names_the_toolkit_that_cannot_be_built_and_why(self, fake_hera):
        fake_hera.missing = {"MeteoLowFreq"}

        problems = ToolkitCheck(REPO_TREE, "Alpha").problems()

        assert problems == ["MeteoLowFreq: ImportError: Cannot locate class: MeteoLowFreq"]


class TestLoadReportsWhatItDid:
    def test_says_how_many_documents_were_added(self, fake_hera, warmed):
        answer = routes.loadIntoProject(_payload(repositoryJSON=REPO_TREE))

        assert answer["added"] == 1
        assert answer["problems"] == []

    def test_a_load_that_added_nothing_says_which_toolkit_was_skipped(self, fake_hera, warmed):
        fake_hera.addsNothing = True
        fake_hera.missing = {"MeteoLowFreq"}

        answer = routes.loadIntoProject(_payload(repositoryJSON=REPO_TREE))

        assert answer["added"] == 0
        assert answer["problems"] == ["MeteoLowFreq: ImportError: Cannot locate class: MeteoLowFreq"]

    def test_a_load_that_added_nothing_with_every_toolkit_fine_says_so_with_no_problem(
            self, fake_hera, warmed):
        fake_hera.addsNothing = True

        answer = routes.loadIntoProject(_payload(repositoryJSON=REPO_TREE))

        assert answer["added"] == 0
        assert answer["problems"] == []

    def test_keeps_the_errors_hera_logged_while_loading(self, fake_hera, warmed, monkeypatch):
        def complain_and_add_nothing(self, projectName, repositoryJSON, basedir="", overwrite=False):
            logging.getLogger("hera.utils.data.toolkit").error(
                "The error X occured while adding *Datasource* to toolkit MeteoLowFreq... skipping!!!")

        toolkit = sys.modules["hera.utils.data.toolkit"].dataToolkit
        monkeypatch.setattr(toolkit, "loadAllDatasourcesInRepositoryJSONToProject",
                            complain_and_add_nothing, raising=True)

        answer = routes.loadIntoProject(_payload(repositoryJSON=REPO_TREE))

        assert answer["added"] == 0
        assert "skipping!!!" in answer["problems"][0]

    def test_the_registered_repositories_are_many_trees_so_no_toolkit_is_checked(
            self, fake_hera, warmed):
        fake_hera.missing = {"MeteoLowFreq"}

        answer = routes.loadIntoProject(_payload())

        assert answer["added"] == 1
        assert answer["problems"] == []


# Every route reaches hera, so each one refuses while hera is still warming.
WARMING_CALLS = {
    "registered": lambda: routes.registered(),
    "get": lambda: routes.get("/repos/meteo.json"),
    "jsonsInFolder": lambda: routes.jsonsInFolder("/repos"),
    "loadIntoProject": lambda: routes.loadIntoProject(_payload()),
}


@pytest.mark.parametrize("route", sorted(WARMING_CALLS))
def test_every_route_says_warming_up(route, fake_hera, monkeypatch):
    monkeypatch.setattr(dependencies.warmup, "_ready", False)

    with pytest.raises(HTTPException) as raised:
        WARMING_CALLS[route]()

    assert raised.value.status_code == 503
    assert raised.value.detail == "WARMING_UP"
