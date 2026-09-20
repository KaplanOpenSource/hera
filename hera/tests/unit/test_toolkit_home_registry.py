"""``ToolkitHome``'s registry surface: the multi-source listing and the two
JSON importers.

``getToolkitDocuments`` is documented as the *single source of truth* for
"what toolkits exist", merging three independent sources into one normalised
shape::

    {"toolkit": <name>,
     "desc": {"classpath", "type", "source", "repositoryName", "version"}}

The three sources are the static in-memory ``_toolkits`` registry (tagged
``source="internal"``), DB-backed ``ToolkitDataSource`` documents written by
``registerToolkit`` (``source="db"``), and the experiments exposed by
``experimentTK`` (``source="experiment"``).  ``getToolkitTable`` is a thin
projection of that list onto six columns, de-duplicated on
``(toolkit, source)``.  Those two properties -- the source tag on every row,
and the fact that all three sources land in one list -- are what the tests
below assert, rather than the exact registry contents.

Isolation notes, because ``ToolkitHome`` is normally a process-wide
singleton: every test here builds its *own* ``ToolkitHome`` bound to the
mongomock-backed unit project, so nothing is patched on ``hera.toolkitHome``
and the shared ``_toolkits`` dict is never mutated.  ``registerToolkit``
writes through ``dataToolkit()``, i.e. into ``defaultProject``, which the
conftest's ``_reset_unit_database`` drops after every test.

Two of the six targeted methods can never succeed, and both fail for the
same reason -- they call ``registerToolkit`` with keyword names that do not
exist in its signature.  Pinned below as B259 (``auto_register_and_get``)
and B260 (``import_toolkits_from_json``).

Deliberately not covered here:

* ``getToolkit`` and ``registerToolkit`` themselves, and the
  ``setDefaultRepository``/``getDefaultRepository`` pair (with its already
  reported B21/B22) -- see ``test_toolkit_home.py``.
* The dynamic-toolkit *loading* branch of ``getToolkit`` (sys.path
  manipulation and ``pydoc.locate`` of an on-disk package): it needs a real
  importable toolkit package on disk, which is an integration concern.
* ``getRiskAreas``-style end-to-end toolkit instantiation: the importers are
  exercised up to the point where they hand off, which is where they break.
"""
import json

import pytest

from hera.tests.unit.conftest import UNIT_PROJECT_NAME


@pytest.fixture()
def home(unit_files_directory):
    """A private ToolkitHome, so the module-level singleton is untouched."""
    from hera.toolkit import ToolkitHome

    return ToolkitHome(projectName=UNIT_PROJECT_NAME, filesDirectory=unit_files_directory)


class _FakeExperiments:
    """Stands in for experimentHome, which only needs getExperimentsMap()."""

    def __init__(self, *names):
        self._names = names

    def getExperimentsMap(self):
        return {name: {} for name in self._names}


def _by_name(docs):
    return {d["toolkit"]: d["desc"] for d in docs}


# ---------------------------------------------------------------------------
# getToolkitDocuments: the static source
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestStaticToolkitDocuments:
    def test_every_publicly_named_toolkit_is_listed(self, home):
        """The class constants are the documented public names, so each one
        has to resolve to a row."""
        names = _by_name(home.getToolkitDocuments())
        for constant in ("GIS_BUILDINGS", "GIS_TILES", "GIS_LANDCOVER",
                         "GIS_VECTOR_TOPOGRAPHY", "GIS_RASTER_TOPOGRAPHY",
                         "GIS_DEMOGRAPHY", "RISKASSESSMENT", "LSM",
                         "SIMULATIONS_WORKFLOWS", "SIMULATIONS_OPENFOAM",
                         "METEOROLOGY_HIGHFREQ", "METEOROLOGY_LOWFREQ",
                         "EXPERIMENT", "WINDPROFILE", "GAUSSIANDISPERSION",
                         "MACHINELEARNING_DEEPLEARNING"):
            assert getattr(type(home), constant) in names

    def test_a_static_row_is_tagged_as_internal(self, home):
        assert _by_name(home.getToolkitDocuments())[home.LSM]["source"] == "internal"

    def test_a_static_row_carries_an_importable_classpath(self, home):
        """The classpath is what getToolkit feeds to pydoc.locate, so it has
        to resolve -- checked here against the import system, not against a
        copied string."""
        import pydoc

        desc = _by_name(home.getToolkitDocuments())[home.RISKASSESSMENT]
        assert pydoc.locate(desc["classpath"]) is not None

    def test_a_static_row_declares_its_domain(self, home):
        docs = _by_name(home.getToolkitDocuments())
        assert docs[home.LSM]["type"] == "simulations"
        assert docs[home.GIS_TILES]["type"] == "measurements"

    def test_static_rows_have_no_repository_and_no_version(self, home):
        desc = _by_name(home.getToolkitDocuments())[home.LSM]
        assert (desc["repositoryName"], desc["version"]) == ("", "")

    def test_every_row_has_the_documented_shape(self, home):
        for doc in home.getToolkitDocuments():
            assert set(doc) == {"toolkit", "desc"}
            assert set(doc["desc"]) == {"classpath", "type", "source",
                                        "repositoryName", "version"}

    def test_the_name_filter_returns_only_that_toolkit(self, home):
        docs = home.getToolkitDocuments(name=home.LSM)
        assert [d["toolkit"] for d in docs] == [home.LSM]

    def test_a_name_that_matches_nothing_gives_an_empty_list(self, home):
        assert home.getToolkitDocuments(name="NO_SUCH_TOOLKIT") == []


# ---------------------------------------------------------------------------
# getToolkitDocuments: the DB-backed source
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestDynamicToolkitDocuments:
    def test_a_registered_toolkit_appears_in_the_listing(self, home, tmp_path):
        directory = tmp_path / "mytk"
        directory.mkdir()
        home.registerToolkit(toolkit_name="MyTK", toolkit_path=str(directory))

        assert "MyTK" in _by_name(home.getToolkitDocuments())

    def test_a_registered_toolkit_is_tagged_as_coming_from_the_db(self, home, tmp_path):
        directory = tmp_path / "mytk"
        directory.mkdir()
        home.registerToolkit(toolkit_name="MyTK", toolkit_path=str(directory))

        assert _by_name(home.getToolkitDocuments())["MyTK"]["source"] == "db"

    def test_a_registered_toolkit_reports_the_version_it_was_given(self, home, tmp_path):
        directory = tmp_path / "mytk"
        directory.mkdir()
        home.registerToolkit(toolkit_name="MyTK", toolkit_path=str(directory),
                             version=(2, 3, 4))

        assert _by_name(home.getToolkitDocuments())["MyTK"]["version"] == (2, 3, 4)

    def test_a_registered_toolkit_defaults_to_the_measurements_domain(self, home, tmp_path):
        """registerToolkit stores no 'type', so the listing has to supply the
        documented fallback."""
        directory = tmp_path / "mytk"
        directory.mkdir()
        home.registerToolkit(toolkit_name="MyTK", toolkit_path=str(directory))

        assert _by_name(home.getToolkitDocuments())["MyTK"]["type"] == "measurements"

    def test_the_classpath_of_a_registered_toolkit_is_its_directory(self, home, tmp_path):
        """Characterisation, not an endorsement: the column is called
        'classpath' but registerToolkit stores a filesystem directory, which
        is what getToolkit later puts on sys.path."""
        directory = tmp_path / "mytk"
        directory.mkdir()
        home.registerToolkit(toolkit_name="MyTK", toolkit_path=str(directory))

        assert _by_name(home.getToolkitDocuments())["MyTK"]["classpath"] == str(directory)

    def test_the_name_filter_also_narrows_the_db_rows(self, home, tmp_path):
        for name in ("TkOne", "TkTwo"):
            directory = tmp_path / name
            directory.mkdir()
            home.registerToolkit(toolkit_name=name, toolkit_path=str(directory))

        assert [d["toolkit"] for d in home.getToolkitDocuments(name="TkTwo")] == ["TkTwo"]

    def test_a_non_tuple_version_is_exploded_into_characters(self, home, tmp_path):
        """Robustness note (not pinned as a bug): the listing coerces the
        stored version with ``tuple(...)``, which is right for a list and
        wrong for anything string-like -- '1.0' becomes ('1', '.', '0').
        Every in-tree writer stores a 3-tuple, so nothing reaches this today.
        """
        from hera.utils.data.toolkit import dataToolkit

        dtk = dataToolkit()
        dtk._allowWritingToDefaultProject = True
        dtk.addDataSource(dataSourceName="StrVer", resource=str(tmp_path),
                          dataFormat="string", version="1.0")

        assert _by_name(home.getToolkitDocuments())["StrVer"]["version"] == ("1", ".", "0")


# ---------------------------------------------------------------------------
# getExperimentToolkitDocuments
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestExperimentToolkitDocuments:
    def test_without_an_experiment_toolkit_there_are_no_experiment_rows(self, home):
        assert home.experimentTK is None
        assert home.getExperimentToolkitDocuments() == []

    def test_each_experiment_becomes_one_row(self, home, monkeypatch):
        monkeypatch.setattr(home, "experimentTK", _FakeExperiments("EXP_A", "EXP_B"))

        docs = home.getExperimentToolkitDocuments()
        assert sorted(d["toolkit"] for d in docs) == ["EXP_A", "EXP_B"]

    def test_an_experiment_row_is_tagged_as_an_experiment(self, home, monkeypatch):
        monkeypatch.setattr(home, "experimentTK", _FakeExperiments("EXP_A"))

        desc = home.getExperimentToolkitDocuments()[0]["desc"]
        assert desc["source"] == "experiment"
        assert desc["type"] == "experiment"

    def test_an_experiment_row_carries_no_classpath(self, home, monkeypatch):
        """Documented: experiments are built by experimentHome.getExperiment,
        never imported by classpath."""
        monkeypatch.setattr(home, "experimentTK", _FakeExperiments("EXP_A"))

        assert home.getExperimentToolkitDocuments()[0]["desc"]["classpath"] == ""

    def test_the_name_filter_selects_a_single_experiment(self, home, monkeypatch):
        monkeypatch.setattr(home, "experimentTK", _FakeExperiments("EXP_A", "EXP_B"))

        docs = home.getExperimentToolkitDocuments(name="EXP_B")
        assert [d["toolkit"] for d in docs] == ["EXP_B"]

    def test_a_broken_experiment_toolkit_does_not_break_the_listing(self, home, monkeypatch):
        """Documented intent: a failure while querying experiments must not
        take down the unified view."""
        class Boom:
            def getExperimentsMap(self):
                raise RuntimeError("the experiment store is down")

        monkeypatch.setattr(home, "experimentTK", Boom())

        assert home.getExperimentToolkitDocuments() == []
        assert home.getToolkitDocuments()  # the static rows still come back

    def test_experiments_are_merged_into_the_unified_listing(self, home, monkeypatch):
        monkeypatch.setattr(home, "experimentTK", _FakeExperiments("EXP_A"))

        docs = _by_name(home.getToolkitDocuments())
        assert docs["EXP_A"]["source"] == "experiment"
        assert docs[home.LSM]["source"] == "internal"


# ---------------------------------------------------------------------------
# getToolkitTable
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetToolkitTable:
    def test_the_documented_columns_are_present_in_order(self, home):
        table = home.getToolkitTable()
        assert list(table.columns) == ["toolkit", "cls", "source", "type",
                                       "repositoryName", "version"]

    def test_there_is_one_row_per_listed_document(self, home):
        assert len(home.getToolkitTable()) == len(home.getToolkitDocuments())

    def test_the_classpath_lands_in_the_cls_column(self, home):
        table = home.getToolkitTable().set_index("toolkit")
        docs = _by_name(home.getToolkitDocuments())
        assert table.loc[home.LSM, "cls"] == docs[home.LSM]["classpath"]

    def test_all_three_sources_reach_the_table(self, home, tmp_path, monkeypatch):
        directory = tmp_path / "mytk"
        directory.mkdir()
        home.registerToolkit(toolkit_name="MyTK", toolkit_path=str(directory))
        monkeypatch.setattr(home, "experimentTK", _FakeExperiments("EXP_A"))

        table = home.getToolkitTable()
        assert set(table["source"]) == {"internal", "db", "experiment"}

    def test_two_versions_of_one_dynamic_toolkit_collapse_to_one_row(self, home, tmp_path):
        """The documented de-duplication key is (toolkit, source)."""
        for suffix, version in (("a", (0, 0, 1)), ("b", (0, 0, 2))):
            directory = tmp_path / suffix
            directory.mkdir()
            home.registerToolkit(toolkit_name="MyTK", toolkit_path=str(directory),
                                 version=version)

        assert len(home.getToolkitDocuments(name="MyTK")) == 2

        table = home.getToolkitTable()
        assert len(table[table["toolkit"] == "MyTK"]) == 1

    def test_the_surviving_duplicate_is_the_first_one(self, home, tmp_path):
        first = tmp_path / "first"
        first.mkdir()
        second = tmp_path / "second"
        second.mkdir()
        home.registerToolkit(toolkit_name="MyTK", toolkit_path=str(first), version=(0, 0, 1))
        home.registerToolkit(toolkit_name="MyTK", toolkit_path=str(second), version=(0, 0, 2))

        table = home.getToolkitTable()
        row = table[table["toolkit"] == "MyTK"].iloc[0]
        assert row["cls"] == str(first)

    def test_an_empty_listing_still_yields_the_documented_columns(self, home, monkeypatch):
        """The no-rows branch builds its own empty frame, so the schema has
        to survive an empty registry."""
        monkeypatch.setattr(type(home), "getToolkitDocuments",
                            lambda self, name=None: [])

        table = home.getToolkitTable()
        assert table.empty
        assert list(table.columns) == ["toolkit", "cls", "source", "type",
                                       "repositoryName", "version"]

    def test_the_project_name_argument_is_accepted_and_ignored(self, home):
        """Documented as 'if None, uses the current project'; the parameter is
        in fact never read, so both calls agree."""
        withProject = home.getToolkitTable(projectName="SOME_OTHER_PROJECT")
        withoutProject = home.getToolkitTable()
        assert list(withProject["toolkit"]) == list(withoutProject["toolkit"])


# ---------------------------------------------------------------------------
# import_experiments_from_json
# ---------------------------------------------------------------------------

def _experiments_json(tmp_path, payload, name="experiments.json"):
    path = tmp_path / name
    path.write_text(json.dumps(payload))
    return str(path)


@pytest.mark.unit
class TestImportExperimentsFromJSON:
    def test_a_project_name_is_required(self, home, tmp_path):
        path = _experiments_json(tmp_path, {"experiments": []})
        with pytest.raises(ValueError, match="projectName is required"):
            home.import_experiments_from_json(projectName="", json_path=path)

    def test_a_mapping_instead_of_a_list_is_rejected(self, home, tmp_path):
        path = _experiments_json(tmp_path, {"experiments": {"EXP": {}}})
        with pytest.raises(ValueError, match="'experiments' must be a list"):
            home.import_experiments_from_json(projectName=UNIT_PROJECT_NAME, json_path=path)

    def test_a_payload_without_experiments_loads_nothing(self, home, tmp_path):
        path = _experiments_json(tmp_path, {})
        assert home.import_experiments_from_json(
            projectName=UNIT_PROJECT_NAME, json_path=path) == []

    def test_the_experiment_names_are_returned(self, home, tmp_path):
        path = _experiments_json(tmp_path, {"experiments": [
            {"name": "EXP_ONE", "data": []},
            {"name": "EXP_TWO", "data": []},
        ]})
        assert home.import_experiments_from_json(
            projectName=UNIT_PROJECT_NAME, json_path=path) == ["EXP_ONE", "EXP_TWO"]

    def test_a_repeated_name_is_reported_once(self, home, tmp_path):
        path = _experiments_json(tmp_path, {"experiments": [
            {"name": "EXP", "data": []},
            {"name": "EXP", "data": []},
        ]})
        assert home.import_experiments_from_json(
            projectName=UNIT_PROJECT_NAME, json_path=path) == ["EXP"]

    def test_an_unnamed_experiment_is_not_reported(self, home, tmp_path):
        path = _experiments_json(tmp_path, {"experiments": [{"data": []}]})
        assert home.import_experiments_from_json(
            projectName=UNIT_PROJECT_NAME, json_path=path) == []

    def test_a_data_item_becomes_a_measurements_document(self, home, tmp_path, unit_project):
        path = _experiments_json(tmp_path, {"experiments": [
            {"name": "EXP", "data": [
                {"type": "Sonic", "resource": "/data/sonic.parquet",
                 "dataFormat": "parquet", "desc": {"station": "A"}},
            ]},
        ]})
        home.import_experiments_from_json(projectName=UNIT_PROJECT_NAME, json_path=path)

        docs = unit_project.getMeasurementsDocuments(type="Sonic")
        assert len(docs) == 1
        assert docs[0]["resource"] == "/data/sonic.parquet"
        assert docs[0]["dataFormat"] == "parquet"
        assert docs[0].desc == {"station": "A"}

    def test_a_relative_resource_is_resolved_against_the_json_file(self, home, tmp_path, unit_project):
        """CLAUDE.md's isRelativePath contract: repository JSON never carries
        absolute paths, so the importer has to anchor them itself."""
        nested = tmp_path / "nested"
        nested.mkdir()
        path = _experiments_json(nested, {"experiments": [
            {"name": "EXP", "data": [
                {"type": "Sonic", "resource": "measurements/sonic.parquet",
                 "dataFormat": "parquet", "isRelativePath": True},
            ]},
        ]})
        home.import_experiments_from_json(projectName=UNIT_PROJECT_NAME, json_path=path)

        stored = unit_project.getMeasurementsDocuments(type="Sonic")[0]["resource"]
        assert stored == str(nested / "measurements" / "sonic.parquet")

    def test_an_absolute_resource_is_stored_verbatim(self, home, tmp_path, unit_project):
        path = _experiments_json(tmp_path, {"experiments": [
            {"name": "EXP", "data": [
                {"type": "Sonic", "resource": "/absolute/sonic.parquet",
                 "dataFormat": "parquet", "isRelativePath": False},
            ]},
        ]})
        home.import_experiments_from_json(projectName=UNIT_PROJECT_NAME, json_path=path)

        assert unit_project.getMeasurementsDocuments(
            type="Sonic")[0]["resource"] == "/absolute/sonic.parquet"

    @pytest.mark.parametrize("missing", ["type", "resource", "dataFormat"])
    def test_an_item_missing_a_mandatory_field_is_skipped(self, home, tmp_path,
                                                          unit_project, missing):
        item = {"type": "Sonic", "resource": "/x", "dataFormat": "parquet"}
        item.pop(missing)
        path = _experiments_json(tmp_path, {"experiments": [{"name": "EXP", "data": [item]}]})

        home.import_experiments_from_json(projectName=UNIT_PROJECT_NAME, json_path=path)

        assert list(unit_project.getMeasurementsDocuments(type="Sonic")) == []

    def test_a_skipped_item_does_not_skip_its_experiment_name(self, home, tmp_path):
        path = _experiments_json(tmp_path, {"experiments": [
            {"name": "EXP", "data": [{"resource": "/x", "dataFormat": "parquet"}]},
        ]})
        assert home.import_experiments_from_json(
            projectName=UNIT_PROJECT_NAME, json_path=path) == ["EXP"]

    def test_several_data_items_all_land(self, home, tmp_path, unit_project):
        path = _experiments_json(tmp_path, {"experiments": [
            {"name": "EXP", "data": [
                {"type": "Sonic", "resource": "/a", "dataFormat": "parquet"},
                {"type": "Sonic", "resource": "/b", "dataFormat": "parquet"},
                {"type": "Tower", "resource": "/c", "dataFormat": "parquet"},
            ]},
        ]})
        home.import_experiments_from_json(projectName=UNIT_PROJECT_NAME, json_path=path)

        assert len(unit_project.getMeasurementsDocuments(type="Sonic")) == 2
        assert len(unit_project.getMeasurementsDocuments(type="Tower")) == 1

    def test_a_missing_desc_becomes_an_empty_description(self, home, tmp_path, unit_project):
        path = _experiments_json(tmp_path, {"experiments": [
            {"name": "EXP", "data": [
                {"type": "Sonic", "resource": "/a", "dataFormat": "parquet"}]},
        ]})
        home.import_experiments_from_json(projectName=UNIT_PROJECT_NAME, json_path=path)

        assert unit_project.getMeasurementsDocuments(type="Sonic")[0].desc == {}


# ---------------------------------------------------------------------------
# import_toolkits_from_json
# ---------------------------------------------------------------------------

RISK_CLASSPATH = "hera.riskassessment.riskToolkit.RiskToolkit"


def _toolkits_json(tmp_path, payload, name="toolkits.json"):
    path = tmp_path / name
    path.write_text(json.dumps(payload))
    return str(path)


@pytest.mark.unit
class TestImportToolkitsFromJSONValidation:
    def test_a_project_name_is_required(self, home, tmp_path):
        path = _toolkits_json(tmp_path, {"repository": "r", "toolkits": []})
        with pytest.raises(ValueError, match="projectName is required"):
            home.import_toolkits_from_json(projectName="", json_path=path)

    def test_without_any_repository_it_refuses_to_guess(self, home, tmp_path):
        """No 'repository' in the JSON, no default_repository argument, and no
        stored default -- there is nothing to register against."""
        path = _toolkits_json(tmp_path, {"toolkits": []})
        with pytest.raises(ValueError, match="No repository provided in JSON"):
            home.import_toolkits_from_json(projectName=UNIT_PROJECT_NAME, json_path=path)

    def test_the_default_repository_argument_satisfies_the_requirement(self, home, tmp_path):
        path = _toolkits_json(tmp_path, {"toolkits": []})
        assert home.import_toolkits_from_json(
            projectName=UNIT_PROJECT_NAME, json_path=path,
            default_repository="repoA") == []

    def test_a_blank_repository_in_the_json_falls_through_to_the_argument(self, home, tmp_path):
        path = _toolkits_json(tmp_path, {"repository": "   ", "toolkits": []})
        assert home.import_toolkits_from_json(
            projectName=UNIT_PROJECT_NAME, json_path=path,
            default_repository="repoA") == []

    def test_a_mapping_instead_of_a_list_is_rejected(self, home, tmp_path):
        path = _toolkits_json(tmp_path, {"repository": "r", "toolkits": {"a": 1}})
        with pytest.raises(ValueError, match="'toolkits' must be a list"):
            home.import_toolkits_from_json(projectName=UNIT_PROJECT_NAME, json_path=path)

    def test_an_empty_list_registers_nothing(self, home, tmp_path):
        path = _toolkits_json(tmp_path, {"repository": "r", "toolkits": []})
        assert home.import_toolkits_from_json(
            projectName=UNIT_PROJECT_NAME, json_path=path) == []

    @pytest.mark.parametrize("entry", [
        {"classpath": RISK_CLASSPATH},
        {"name": "MyTK"},
        {},
    ])
    def test_an_entry_missing_a_name_or_a_classpath_is_rejected(self, home, tmp_path, entry):
        path = _toolkits_json(tmp_path, {"repository": "r", "toolkits": [entry]})
        with pytest.raises(ValueError, match="missing name/classpath"):
            home.import_toolkits_from_json(projectName=UNIT_PROJECT_NAME, json_path=path)

    def test_a_classpath_that_does_not_resolve_is_reported(self, home, tmp_path):
        path = _toolkits_json(tmp_path, {
            "repository": "r",
            "toolkits": [{"name": "MyTK", "classpath": "no.such.module.Class"}]})
        with pytest.raises(ImportError, match="Cannot locate class by classpath"):
            home.import_toolkits_from_json(projectName=UNIT_PROJECT_NAME, json_path=path)


@pytest.mark.unit
class TestImportToolkitsFromJSONIsBroken:
    """B260: see the module docstring."""

    @staticmethod
    def _one_toolkit(tmp_path):
        return _toolkits_json(tmp_path, {
            "repository": "repoA",
            "toolkits": [{"name": "MyTK", "classpath": RISK_CLASSPATH,
                          "version": [1, 2, 3]}]})

    @pytest.mark.xfail(
        strict=True,
        reason="B260: import_toolkits_from_json calls "
               "self.registerToolkit(toolkitclass=..., repositoryName=..., "
               "datasource_name=..., params=..., version=..., overwrite=...), "
               "but registerToolkit's signature is "
               "(toolkit_name, toolkit_path, params, version, overwrite, "
               "**kwargs).  None of the three leading keywords match, so "
               "every entry raises TypeError: registerToolkit() missing 2 "
               "required positional arguments: 'toolkit_name' and "
               "'toolkit_path'.  The method can never register anything. "
               "See the consolidated findings issue.",
    )
    def test_a_declared_toolkit_should_be_registered(self, home, tmp_path):
        assert home.import_toolkits_from_json(
            projectName=UNIT_PROJECT_NAME,
            json_path=self._one_toolkit(tmp_path)) == ["MyTK"]

    def test_it_currently_raises_a_type_error(self, home, tmp_path):
        """Characterisation of B260."""
        with pytest.raises(TypeError, match="missing 2 required positional arguments"):
            home.import_toolkits_from_json(projectName=UNIT_PROJECT_NAME,
                                           json_path=self._one_toolkit(tmp_path))

    def test_nothing_reaches_the_registry(self, home, tmp_path):
        """Characterisation of B260: the failure is total, not partial."""
        with pytest.raises(TypeError):
            home.import_toolkits_from_json(projectName=UNIT_PROJECT_NAME,
                                           json_path=self._one_toolkit(tmp_path))

        assert home.getToolkitDocuments(name="MyTK") == []

    def test_the_validation_ahead_of_it_still_runs_first(self, home, tmp_path):
        """Characterisation of B260's blast radius: the guards before the
        registerToolkit call are reachable, which is why they are asserted
        above rather than xfailed."""
        path = _toolkits_json(tmp_path, {
            "repository": "r",
            "toolkits": [{"name": "MyTK", "classpath": "no.such.Class"}]})
        with pytest.raises(ImportError):
            home.import_toolkits_from_json(projectName=UNIT_PROJECT_NAME, json_path=path)


# ---------------------------------------------------------------------------
# auto_register_and_get
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestAutoRegisterAndGetValidation:
    def test_without_a_repository_json_there_is_no_classpath_hint(self, home):
        with pytest.raises(ValueError, match="no classpath hint found"):
            home.auto_register_and_get("MyTK")

    def test_a_repository_json_without_a_registry_block_is_rejected(self, home):
        with pytest.raises(ValueError, match="no classpath hint found"):
            home.auto_register_and_get("MyTK", repositoryJSON={"MyTK": {"Config": {}}})

    def test_a_hint_for_a_different_toolkit_does_not_count(self, home):
        with pytest.raises(ValueError, match="no classpath hint found"):
            home.auto_register_and_get(
                "MyTK", repositoryJSON={"OtherTK": {"Registry": {"classpath": RISK_CLASSPATH}}})

    def test_a_hint_without_a_dot_cannot_be_split_into_module_and_class(self, home):
        with pytest.raises(ValueError, match="Invalid classpath hint"):
            home.auto_register_and_get(
                "MyTK", repositoryJSON={"MyTK": {"Registry": {"classpath": "NoDots"}}})

    def test_an_unimportable_hint_is_reported_as_an_import_error(self, home):
        with pytest.raises(ImportError, match="Failed to import 'no.such.Class'"):
            home.auto_register_and_get(
                "MyTK", repositoryJSON={"MyTK": {"Registry": {"cls": "no.such.Class"}}})

    def test_the_cls_key_is_accepted_as_an_alias_for_classpath(self, home):
        """Both spellings are documented; reaching the repository check proves
        the alias resolved."""
        with pytest.raises(ValueError, match="no target repository"):
            home.auto_register_and_get(
                "MyTK", repositoryJSON={"MyTK": {"Registry": {"cls": RISK_CLASSPATH}}})

    def test_without_a_project_and_without_a_repository_it_refuses(self, home, monkeypatch):
        monkeypatch.setattr(home, "_projectName", None)
        with pytest.raises(ValueError, match="projectName is None"):
            home.auto_register_and_get(
                "MyTK", repositoryJSON={"MyTK": {"Registry": {"classpath": RISK_CLASSPATH}}})

    def test_an_unset_default_repository_is_reported(self, home):
        """No default repository is stored for the unit project, so the
        lookup comes back empty and the call stops there."""
        with pytest.raises(ValueError, match="no target repository"):
            home.auto_register_and_get(
                "MyTK", repositoryJSON={"MyTK": {"Registry": {"classpath": RISK_CLASSPATH}}})

    def test_a_whitespace_only_repository_name_counts_as_unset(self, home):
        with pytest.raises(ValueError, match="no target repository"):
            home.auto_register_and_get(
                "MyTK",
                repositoryJSON={"MyTK": {"Registry": {"classpath": RISK_CLASSPATH}}},
                repositoryName="   ")


@pytest.mark.unit
class TestAutoRegisterAndGetIsBroken:
    """B259: see the module docstring."""

    HINT = {"MyTK": {"Registry": {"classpath": RISK_CLASSPATH}}}

    @pytest.mark.xfail(
        strict=True,
        reason="B259: auto_register_and_get calls "
               "self.registerToolkit(toolkitclass=..., datasource_name=..., "
               "params=..., version=..., overwrite=True, repositoryName=...), "
               "but registerToolkit's signature is "
               "(toolkit_name, toolkit_path, params, version, overwrite, "
               "**kwargs).  toolkitclass/datasource_name are absorbed by "
               "**kwargs while the two required positional parameters go "
               "unfilled, so the call raises TypeError: registerToolkit() "
               "missing 2 required positional arguments: 'toolkit_name' and "
               "'toolkit_path'.  The method can never return a toolkit. "
               "See the consolidated findings issue.",
    )
    def test_it_should_return_an_instance_of_the_hinted_class(self, home):
        from hera.riskassessment.riskToolkit import RiskToolkit

        toolkit = home.auto_register_and_get(
            "MyTK", repositoryJSON=self.HINT, repositoryName="repoA")
        assert isinstance(toolkit, RiskToolkit)

    def test_it_currently_raises_a_type_error(self, home):
        """Characterisation of B259."""
        with pytest.raises(TypeError, match="missing 2 required positional arguments"):
            home.auto_register_and_get("MyTK", repositoryJSON=self.HINT,
                                       repositoryName="repoA")

    def test_nothing_is_registered_before_it_fails(self, home):
        """Characterisation of B259: registration and instantiation both
        sit behind the broken call, so the DB is untouched."""
        with pytest.raises(TypeError):
            home.auto_register_and_get("MyTK", repositoryJSON=self.HINT,
                                       repositoryName="repoA")

        assert home.getToolkitDocuments(name="MyTK") == []

    def test_the_same_wall_is_hit_with_explicit_params_and_version(self, home):
        """Characterisation of B259: none of the caller-supplied arguments
        change the outcome."""
        with pytest.raises(TypeError, match="missing 2 required positional arguments"):
            home.auto_register_and_get("MyTK", repositoryJSON=self.HINT,
                                       repositoryName="repoA",
                                       params={"a": 1}, version=(9, 9, 9))
