"""LSM/template.py: the four members test_lsm_template.py left out --
``run``, ``_toNetcdf``, ``_getSimulationsList`` and ``getSimulationsTable``.

``run`` normally shells out to a compiled Fortran binary (``./a.out``).
Here ``subprocess.run`` is replaced with a recording stand-in and
``InputForModelsCreator`` (imported into this module's namespace) with a
fake that records the parameter map and writes a placeholder INPUT file.
What is left under test is exactly what LSMTemplate owns: which flags it
derives from its arguments, what it copies where, which save modes it
refuses, which auxiliary input files it writes, and whether it restores
the working directory.  Nothing runs a solver, touches MongoDB outside
mongomock, or writes outside tmp_path.

``_toNetcdf`` is driven from hand-written OUTD files, so the parsing, the
time ordering, the leading zero frame and the unit scaling are all
checkable without a simulation.

The template document is built by hand for the same reason
test_lsm_template.py does: ``LSMToolkit.loadData`` never persists the
parsed template JSON's own desc (B92), so ``params``/``modelFolder``
would KeyError on a template registered the normal way.

Deliberately not covered: the ``to_xarray=True`` tail of ``run`` (the
netcdf assembly), which needs real OUTD output from a real solver run --
``_toNetcdf`` is tested directly instead.

Bugs pinned here, each with an xfail(strict) for the intended behaviour
and a passing characterisation of what actually happens:

  * B275: ``_toNetcdf`` says it converts the dosage "to s/m**3 instead
    of min/m**3" but multiplies by ``(1 s).m_as(min)`` = 1/60, which is
    the seconds-to-minutes factor -- the reverse of the stated conversion.
  * B276: ``_getSimulationsList`` joins its desc frames to its params
    frames with ``itertools.product`` instead of pairing them, so n
    simulations produce n^2 rows in which one simulation's metadata is
    stapled to another's parameters.
  * B277: ``run(topography=..., canopy=...)`` without ``stations``
    reaches ``stations.columns`` on None and dies with AttributeError.
  * B278: ``run`` chdir's into the simulation directory with no
    try/finally, so any failure past that point leaves the whole process
    in the simulation directory.
  * B279: ``run(stations=...)`` computes the station identifier column
    and then loses it in the xarray resample round-trip, so it only
    works if the caller happens to have supplied a ``station`` column.
  * B280: the same code path joins a never-deduplicated ``onlyStations``
    back onto the resampled series, so every station's wind record is
    written once per original observation.
"""
import os
import subprocess

import pandas
import pytest

import hera.simulations.LSM.template as template_mod
from hera import toolkit as toolkit_mod
from hera.simulations.LSM.template import LSMTemplate

BASE_PARAMS = dict(
    TopoXmin=0, TopoXmax=100, TopoYmin=0, TopoYmax=200,
    sourceRatioX=0.5, sourceRatioY=0.25, particles3D=".true.",
)
TEMPLATE_VERSION = "v4-general"
TEMPLATE_NAME = "myTemplate"


class FakeInputForModelsCreator:
    """Records the parameter map; writes a placeholder INPUT file."""

    instances = []

    def __init__(self, templatesDir):
        self.templatesDir = templatesDir
        self.paramsMap = None
        self.templateName = None
        self.renderedTo = []
        FakeInputForModelsCreator.instances.append(self)

    def setParamsMap(self, paramsMap):
        self.paramsMap = dict(paramsMap)

    def setTemplate(self, templateName):
        self.templateName = templateName

    def render(self, savePath=None):
        self.renderedTo.append(savePath)
        with open(savePath, "w") as inputFile:
            inputFile.write("INPUT\n")


@pytest.fixture()
def lsm(unit_toolkit_factory):
    from hera import toolkitHome

    toolkit = unit_toolkit_factory(toolkitHome.LSM)
    toolkit.to_xarray = False
    return toolkit


@pytest.fixture()
def modelFolder(tmp_path):
    folder = tmp_path / "model"
    (folder / "tozaot" / "machsan").mkdir(parents=True)
    (folder / "tozaot" / "Meteorology").mkdir(parents=True)
    (folder / "modelData.txt").write_text("shipped with the template\n")
    return folder


@pytest.fixture()
def template(lsm, tmp_path, modelFolder):
    def _build(**paramOverrides):
        params = dict(BASE_PARAMS)
        params.update(paramOverrides)
        document = {
            "resource": str(tmp_path / "templateDir"),
            "desc": dict(
                params=params,
                version=TEMPLATE_VERSION,
                datasourceName=TEMPLATE_NAME,
                modelFolder=str(modelFolder),
            ),
        }
        return LSMTemplate(document, lsm)

    return _build


@pytest.fixture()
def ifmc(monkeypatch):
    FakeInputForModelsCreator.instances.clear()
    monkeypatch.setattr(
        template_mod, "InputForModelsCreator", FakeInputForModelsCreator
    )
    return FakeInputForModelsCreator


@pytest.fixture()
def solver(monkeypatch):
    """Replace the ./a.out subprocess with a recorder."""
    record = dict(calls=[], returncode=0)

    def _run(command, *args, **kwargs):
        record["calls"].append((list(command), os.getcwd()))

        class _Completed:
            returncode = record["returncode"]

        return _Completed()

    monkeypatch.setattr(subprocess, "run", _run)
    return record


@pytest.fixture()
def cwd(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return tmp_path


def _stationsFrame(withStationColumn=True, samples=3):
    """Two stations, `samples` observations each, ten minutes apart.

    The wind speed rises 1 m/s per observation so that individual records
    are distinguishable in the meteorology files that get written.
    """
    times = pandas.date_range("2020-01-01 00:00", periods=samples, freq="10Min")
    rows = []
    for x, y in ((0.0, 0.0), (50.0, 60.0)):
        for index, moment in enumerate(times):
            row = dict(
                x=x, y=y, datetime=moment,
                u=float(index + 1), direction=270.0, h=2.0,
            )
            if withStationColumn:
                row["station"] = f"S{int(x)}"
            rows.append(row)
    return pandas.DataFrame(rows)


# ---------------------------------------------------------------------------
# run: the parameter map it hands to the input-file writer
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestRunDerivedFlags:
    def _run(self, template, **kwargs):
        arguments = dict(
            saveMode=toolkit_mod.TOOLKIT_SAVEMODE_NOSAVE, simulationName="sim"
        )
        arguments.update(kwargs)
        return template().run(**arguments)

    def test_no_topography_and_no_stations_means_homogeneous_wind(
        self, template, ifmc, solver, cwd
    ):
        self._run(template)
        assert ifmc.instances[0].paramsMap["homogeneousWind"] == ".TRUE."

    def test_topography_switches_on_the_topography_file(
        self, template, ifmc, solver, cwd
    ):
        self._run(template, topography="TOPOGRAPHY DATA")
        paramsMap = ifmc.instances[0].paramsMap
        assert paramsMap["TopoFile"] == "'TOPO'"
        assert paramsMap["flat"] == ".FALSE."
        assert "homogeneousWind" not in paramsMap

    def test_topography_is_written_into_the_simulation_directory(
        self, template, ifmc, solver, cwd, lsm
    ):
        self._run(template, topography="TOPOGRAPHY DATA")
        written = os.path.join(os.path.abspath(lsm.filesDirectory), "sim", "TOPO")
        with open(written) as topoFile:
            assert topoFile.read() == "TOPOGRAPHY DATA"

    def test_no_canopy_switches_the_canopy_off(self, template, ifmc, solver, cwd):
        self._run(template)
        assert ifmc.instances[0].paramsMap["canopy"] == ".FALSE."

    def test_a_canopy_switches_the_canopy_on(self, template, ifmc, solver, cwd):
        self._run(template, canopy="CANOPY DATA")
        assert ifmc.instances[0].paramsMap["canopy"] == ".TRUE."

    def test_the_canopy_properties_are_written(self, template, ifmc, solver, cwd, lsm):
        self._run(template, canopy="CANOPY DATA")
        written = os.path.join(
            os.path.abspath(lsm.filesDirectory), "sim", "canopy_properties.txt"
        )
        with open(written) as canopyFile:
            assert canopyFile.read() == "CANOPY DATA"

    def test_a_single_deposition_rate_is_duplicated_into_two(
        self, template, ifmc, solver, cwd
    ):
        self._run(template, depositionRates=0.01)
        assert ifmc.instances[0].paramsMap["n_vdep"] == 2

    def test_a_list_of_deposition_rates_is_counted_as_given(
        self, template, ifmc, solver, cwd
    ):
        self._run(template, depositionRates=[0.01, 0.02, 0.03])
        assert ifmc.instances[0].paramsMap["n_vdep"] == 3

    def test_the_deposition_rates_are_written_one_per_line(
        self, template, ifmc, solver, cwd, lsm
    ):
        self._run(template, depositionRates=[0.01, 0.02])
        written = os.path.join(
            os.path.abspath(lsm.filesDirectory), "sim", "INPUT_VDEP"
        )
        with open(written) as depositionFile:
            assert depositionFile.read() == "0.01\n0.02\n"

    def test_the_callers_parameters_override_the_templates_own(
        self, template, ifmc, solver, cwd
    ):
        self._run(template, params=dict(TopoXmax=999))
        assert ifmc.instances[0].paramsMap["TopoXmax"] == 999

    def test_extra_descriptors_reach_the_parameter_map(
        self, template, ifmc, solver, cwd
    ):
        self._run(template, TopoYmax=777)
        assert ifmc.instances[0].paramsMap["TopoYmax"] == 777

    def test_the_templates_own_defaults_are_kept(self, template, ifmc, solver, cwd):
        self._run(template)
        for key, value in BASE_PARAMS.items():
            assert ifmc.instances[0].paramsMap[key] == value


@pytest.mark.unit
class TestRunFileLayout:
    def _run(self, template, **kwargs):
        arguments = dict(
            saveMode=toolkit_mod.TOOLKIT_SAVEMODE_NOSAVE, simulationName="sim"
        )
        arguments.update(kwargs)
        return template().run(**arguments)

    def _saveDir(self, lsm, name="sim"):
        return os.path.join(os.path.abspath(lsm.filesDirectory), name)

    def test_the_input_writer_is_pointed_at_the_templates_directory(
        self, template, ifmc, solver, cwd, tmp_path
    ):
        self._run(template)
        assert ifmc.instances[0].templatesDir == str(tmp_path / "templateDir")

    def test_the_input_template_is_chosen_by_the_template_version(
        self, template, ifmc, solver, cwd
    ):
        self._run(template)
        assert ifmc.instances[0].templateName == f"LSM_{TEMPLATE_VERSION}"

    def test_the_input_file_is_rendered_into_the_simulation_directory(
        self, template, ifmc, solver, cwd, lsm
    ):
        self._run(template)
        assert ifmc.instances[0].renderedTo == [
            os.path.join(self._saveDir(lsm), "INPUT")
        ]

    def test_the_model_folder_contents_are_copied(
        self, template, ifmc, solver, cwd, lsm
    ):
        self._run(template)
        saveDir = self._saveDir(lsm)
        assert os.path.isfile(os.path.join(saveDir, "modelData.txt"))
        assert os.path.isdir(os.path.join(saveDir, "tozaot", "machsan"))

    def test_the_copied_files_keep_their_contents(
        self, template, ifmc, solver, cwd, lsm
    ):
        self._run(template)
        with open(os.path.join(self._saveDir(lsm), "modelData.txt")) as copied:
            assert copied.read() == "shipped with the template\n"

    def test_the_model_folder_listing_is_printed(
        self, template, ifmc, solver, cwd, capsys
    ):
        self._run(template)
        assert "modelData.txt" in capsys.readouterr().out

    def test_the_solver_is_launched_from_the_simulation_directory(
        self, template, ifmc, solver, cwd, lsm
    ):
        self._run(template)
        command, workingDirectory = solver["calls"][0]
        assert command == ["./a.out"]
        assert workingDirectory == self._saveDir(lsm)

    def test_the_working_directory_is_restored_afterwards(
        self, template, ifmc, solver, cwd
    ):
        self._run(template)
        assert os.getcwd() == str(cwd)

    def test_the_deprecated_saveDir_argument_is_ignored(
        self, template, ifmc, solver, cwd, lsm, tmp_path
    ):
        """Documented as deprecated: the toolkit's filesDirectory always wins."""
        elsewhere = tmp_path / "ignoredLocation"
        elsewhere.mkdir()
        self._run(template, saveDir=str(elsewhere))
        assert list(elsewhere.iterdir()) == []
        assert os.path.isfile(os.path.join(self._saveDir(lsm), "INPUT"))

    def test_without_a_simulation_name_the_files_directory_is_used_directly(
        self, template, ifmc, solver, cwd, lsm
    ):
        template().run(saveMode=toolkit_mod.TOOLKIT_SAVEMODE_NOSAVE)
        assert ifmc.instances[0].renderedTo == [
            os.path.join(os.path.abspath(lsm.filesDirectory), "INPUT")
        ]

    def test_it_returns_nothing_when_xarray_output_is_switched_off(
        self, template, ifmc, solver, cwd
    ):
        assert self._run(template) is None


@pytest.mark.unit
class TestRunSolverFailure:
    def test_a_non_zero_exit_code_abandons_the_run(
        self, template, ifmc, solver, cwd
    ):
        solver["returncode"] = 3
        result = template().run(
            saveMode=toolkit_mod.TOOLKIT_SAVEMODE_NOSAVE, simulationName="sim"
        )
        assert result is None

    def test_a_non_zero_exit_code_still_restores_the_working_directory(
        self, template, ifmc, solver, cwd
    ):
        solver["returncode"] = 3
        template().run(
            saveMode=toolkit_mod.TOOLKIT_SAVEMODE_NOSAVE, simulationName="sim"
        )
        assert os.getcwd() == str(cwd)

    def test_a_non_zero_exit_code_produces_no_netcdf(
        self, template, ifmc, solver, cwd, lsm
    ):
        solver["returncode"] = 3
        lsm.to_xarray = True
        template().run(
            saveMode=toolkit_mod.TOOLKIT_SAVEMODE_NOSAVE, simulationName="sim"
        )
        saveDir = os.path.join(os.path.abspath(lsm.filesDirectory), "sim")
        assert not os.path.exists(os.path.join(saveDir, "netcdf"))


@pytest.mark.unit
class TestRunSaveModes:
    def test_file_and_db_registers_a_simulation_document(
        self, template, ifmc, solver, cwd, lsm
    ):
        template().run(saveMode=toolkit_mod.TOOLKIT_SAVEMODE_FILEANDDB)
        documents = lsm.getSimulationsDocuments(type="LSM_run")
        assert len(documents) == 1
        assert documents[0].desc["templateName"] == TEMPLATE_NAME
        assert documents[0].desc["version"] == TEMPLATE_VERSION

    def test_the_simulation_is_named_from_the_project_counter(
        self, template, ifmc, solver, cwd, lsm
    ):
        template().run(saveMode=toolkit_mod.TOOLKIT_SAVEMODE_FILEANDDB)
        documents = lsm.getSimulationsDocuments(type="LSM_run")
        assert documents[0].desc["simulationName"] == "LSM_Simulation_0"

    def test_an_explicit_simulation_name_wins_over_the_counter(
        self, template, ifmc, solver, cwd, lsm
    ):
        template().run(
            saveMode=toolkit_mod.TOOLKIT_SAVEMODE_FILEANDDB, simulationName="myRun"
        )
        documents = lsm.getSimulationsDocuments(type="LSM_run")
        assert documents[0].desc["simulationName"] == "myRun"

    def test_the_document_records_the_directory_when_xarray_is_off(
        self, template, ifmc, solver, cwd, lsm
    ):
        template().run(
            saveMode=toolkit_mod.TOOLKIT_SAVEMODE_FILEANDDB, simulationName="myRun"
        )
        document = lsm.getSimulationsDocuments(type="LSM_run")[0]
        assert document.resource == os.path.join(
            os.path.abspath(lsm.filesDirectory), "myRun"
        )
        assert document.dataFormat == "string"

    def test_repeating_an_identical_run_is_refused(
        self, template, ifmc, solver, cwd, lsm
    ):
        template().run(saveMode=toolkit_mod.TOOLKIT_SAVEMODE_FILEANDDB)
        with pytest.raises(ValueError, match="already exists in the database"):
            template().run(saveMode=toolkit_mod.TOOLKIT_SAVEMODE_FILEANDDB)

    def test_the_refusal_names_the_replace_mode_to_use_instead(
        self, template, ifmc, solver, cwd, lsm
    ):
        template().run(saveMode=toolkit_mod.TOOLKIT_SAVEMODE_FILEANDDB)
        with pytest.raises(
            ValueError, match=toolkit_mod.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE
        ):
            template().run(saveMode=toolkit_mod.TOOLKIT_SAVEMODE_FILEANDDB)

    def test_the_replace_mode_leaves_a_single_document(
        self, template, ifmc, solver, cwd, lsm
    ):
        template().run(saveMode=toolkit_mod.TOOLKIT_SAVEMODE_FILEANDDB)
        template().run(saveMode=toolkit_mod.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE)
        assert len(lsm.getSimulationsDocuments(type="LSM_run")) == 1

    def test_no_save_writes_no_document(self, template, ifmc, solver, cwd, lsm):
        template().run(
            saveMode=toolkit_mod.TOOLKIT_SAVEMODE_NOSAVE, simulationName="sim"
        )
        assert len(lsm.getSimulationsDocuments(type="LSM_run")) == 0

    def test_only_file_refuses_to_overwrite_an_existing_netcdf_directory(
        self, template, ifmc, solver, cwd, lsm
    ):
        saveDir = os.path.join(os.path.abspath(lsm.filesDirectory), "sim")
        os.makedirs(os.path.join(saveDir, "netcdf"))
        with pytest.raises(ValueError, match="Either remove it"):
            template().run(
                saveMode=toolkit_mod.TOOLKIT_SAVEMODE_ONLYFILE, simulationName="sim"
            )

    def test_only_file_replace_clears_the_existing_netcdf_directory(
        self, template, ifmc, solver, cwd, lsm
    ):
        saveDir = os.path.join(os.path.abspath(lsm.filesDirectory), "sim")
        os.makedirs(os.path.join(saveDir, "netcdf"))
        with open(os.path.join(saveDir, "netcdf", "stale.nc"), "w") as stale:
            stale.write("stale\n")

        template().run(
            saveMode=toolkit_mod.TOOLKIT_SAVEMODE_ONLYFILE_REPLACE,
            simulationName="sim",
        )
        assert not os.path.exists(os.path.join(saveDir, "netcdf", "stale.nc"))


@pytest.mark.unit
class TestRunWithStations:
    def _run(self, template, stations, **kwargs):
        arguments = dict(
            saveMode=toolkit_mod.TOOLKIT_SAVEMODE_NOSAVE, simulationName="sim"
        )
        arguments.update(kwargs)
        return template().run(stations=stations, **arguments)

    def _saveDir(self, lsm):
        return os.path.join(os.path.abspath(lsm.filesDirectory), "sim")

    def test_stations_switch_off_the_homogeneous_wind(
        self, template, ifmc, solver, cwd
    ):
        self._run(template, _stationsFrame())
        paramsMap = ifmc.instances[0].paramsMap
        assert paramsMap["homogeneousWind"] == ".FALSE."
        assert paramsMap["StationsFile"] == "'STATIONS'"

    def test_the_station_index_file_counts_the_stations(
        self, template, ifmc, solver, cwd, lsm
    ):
        self._run(template, _stationsFrame())
        with open(os.path.join(self._saveDir(lsm), "STATIONS")) as stationsFile:
            assert stationsFile.read().splitlines()[0] == "2"

    def test_the_station_index_file_lists_each_stations_coordinates(
        self, template, ifmc, solver, cwd, lsm
    ):
        self._run(template, _stationsFrame())
        with open(os.path.join(self._saveDir(lsm), "STATIONS")) as stationsFile:
            lines = stationsFile.read().splitlines()
        assert lines[1].split() == ["AA", "0.0", "0.0"]
        assert lines[2].split() == ["AB", "50.0", "60.0"]

    def test_one_meteorology_file_is_written_per_station(
        self, template, ifmc, solver, cwd, lsm
    ):
        self._run(template, _stationsFrame())
        metDir = os.path.join(self._saveDir(lsm), "tozaot", "Meteorology")
        assert sorted(os.listdir(metDir)) == ["AA_st.txt", "AB_st.txt"]

    def test_each_meteorology_file_is_padded_before_the_real_record(
        self, template, ifmc, solver, cwd, lsm
    ):
        """2901 zero rows precede the interpolated wind, per the writer."""
        self._run(template, _stationsFrame())
        with open(self._metFile(lsm)) as metFile:
            lines = metFile.read().splitlines()
        assert lines[:2901] == ["0 0"] * 2901
        assert lines[2901].split() == ["1.0", "270.0"]

    def _metFile(self, lsm, name="AA"):
        return os.path.join(
            self._saveDir(lsm), "tozaot", "Meteorology", f"{name}_st.txt"
        )

    def _windRecords(self, lsm, name="AA"):
        with open(self._metFile(lsm, name)) as metFile:
            lines = metFile.read().splitlines()
        assert lines[:2901] == ["0 0"] * 2901
        return [float(line.split()[0]) for line in lines[2901:]]

    @pytest.mark.xfail(
        strict=True,
        reason="B280: after resampling the station series to five-minute "
               "steps, run() joins `onlyStations` back on (x, y) -- but "
               "onlyStations was never deduplicated, so it still holds one "
               "row per original observation per station. The join is "
               "therefore many-to-many and every resampled record is "
               "repeated once per original observation: three 10-minute "
               "observations become 5 interpolated steps written 3 times "
               "each, 15 records instead of 5. The station wind series the "
               "solver reads is the interpolated one stuttered N times, so "
               "it runs N times too long and at the wrong times. "
               "drop_duplicates on the (x, y) key is what the join needs. "
               "See the consolidated findings issue.",
    )
    def test_the_wind_is_resampled_to_five_minute_steps(
        self, template, ifmc, solver, cwd, lsm
    ):
        """Three observations 10 minutes apart span 20 minutes, i.e. five
        five-minute steps, interpolating 1 -> 3 m/s."""
        self._run(template, _stationsFrame(samples=3))
        assert self._windRecords(lsm) == pytest.approx([1.0, 1.5, 2.0, 2.5, 3.0])

    def test_every_resampled_record_is_currently_repeated_once_per_observation(
        self, template, ifmc, solver, cwd, lsm
    ):
        """Characterisation of B280."""
        self._run(template, _stationsFrame(samples=3))
        records = self._windRecords(lsm)
        assert len(records) == 5 * 3
        assert records == pytest.approx(
            [value for value in (1.0, 1.5, 2.0, 2.5, 3.0) for _ in range(3)]
        )

    def test_the_repetition_count_currently_follows_the_observation_count(
        self, template, ifmc, solver, cwd, lsm
    ):
        """Characterisation of B280: two observations, each step twice."""
        self._run(template, _stationsFrame(samples=2))
        records = self._windRecords(lsm)
        assert len(records) == 3 * 2
        assert records == pytest.approx(
            [value for value in (1.0, 1.5, 2.0) for _ in range(2)]
        )

    @pytest.mark.xfail(
        strict=True,
        reason="B279: run() computes the station identifier itself "
               "(stations[stationColumn] = x.astype(str) + y.astype(str)), "
               "but then rebuilds the frame from a resampled xarray "
               "(set_index(...).drop(columns=stationColumn).to_xarray() -> "
               "to_dataframe().reset_index()), which does not carry that "
               "column, and joins back only `onlyStations`, computed before "
               "the identifier existed. The very next line reads "
               "stations['station'] and raises KeyError. The path works only "
               "if the caller happens to have supplied a 'station' column of "
               "their own -- which the docstring never mentions and the "
               "computed identifier is meant to make unnecessary. See the "
               "consolidated findings issue.",
    )
    def test_stations_without_a_station_column_should_still_work(
        self, template, ifmc, solver, cwd, lsm
    ):
        self._run(template, _stationsFrame(withStationColumn=False))
        assert os.path.isfile(os.path.join(self._saveDir(lsm), "STATIONS"))

    def test_stations_without_a_station_column_currently_raise_keyerror(
        self, template, ifmc, solver, cwd
    ):
        """Characterisation of B279."""
        with pytest.raises(KeyError, match="station"):
            self._run(template, _stationsFrame(withStationColumn=False))

    def test_the_computed_identifier_is_currently_discarded(
        self, template, ifmc, solver, cwd, lsm
    ):
        """Characterisation of B279: with a caller-supplied 'station'
        column, it is the caller's values that survive, not the x+y
        identifier the code builds."""
        stations = _stationsFrame()
        self._run(template, stations)
        metDir = os.path.join(self._saveDir(lsm), "tozaot", "Meteorology")
        # two distinct caller stations -> two files, named AA/AB from the
        # loop counter, not from the discarded "0.00.0"/"50.060.0" ids
        assert sorted(os.listdir(metDir)) == ["AA_st.txt", "AB_st.txt"]

    def test_a_canopy_needs_the_station_heights(
        self, template, ifmc, solver, cwd, lsm
    ):
        self._run(template, _stationsFrame(), canopy="CANOPY DATA")
        with open(os.path.join(self._saveDir(lsm), "h_stations.txt")) as heights:
            assert heights.read().split() == ["2.0", "2.0"]

    def test_a_canopy_without_station_heights_is_refused(
        self, template, ifmc, solver, cwd
    ):
        stations = _stationsFrame().drop(columns=["h"])
        with pytest.raises(KeyError, match="Height column is missing"):
            self._run(template, stations, canopy="CANOPY DATA")


@pytest.mark.unit
class TestRunCanopyWithoutStations:
    def _run(self, template, **kwargs):
        arguments = dict(
            saveMode=toolkit_mod.TOOLKIT_SAVEMODE_NOSAVE, simulationName="sim"
        )
        arguments.update(kwargs)
        return template().run(**arguments)

    def test_a_canopy_alone_is_accepted(self, template, ifmc, solver, cwd, lsm):
        """With neither topography nor stations the height check is skipped."""
        assert self._run(template, canopy="CANOPY DATA") is None

    @pytest.mark.xfail(
        strict=True,
        reason="B277: with a canopy and topography but no stations, run() "
               "enters `if (topography is not None or stations is not None)` "
               "and then reads `stations.columns` -- on None. Topography "
               "alone is enough to enter the branch, so a canopy over "
               "terrain without station data always dies with "
               "AttributeError: 'NoneType' object has no attribute "
               "'columns'. Only the stations arm of that `or` justifies the "
               "check. See the consolidated findings issue.",
    )
    def test_a_canopy_over_topography_should_not_need_stations(
        self, template, ifmc, solver, cwd
    ):
        assert self._run(template, topography="TOPO DATA", canopy="CANOPY DATA") is None

    def test_a_canopy_over_topography_currently_raises_attributeerror(
        self, template, ifmc, solver, cwd
    ):
        """Characterisation of B277."""
        with pytest.raises(AttributeError, match="'NoneType' object has no attribute 'columns'"):
            self._run(template, topography="TOPO DATA", canopy="CANOPY DATA")

    @pytest.mark.xfail(
        strict=True,
        reason="B278: run() does os.chdir(saveDir) with no try/finally, "
               "so the restoring os.chdir(cur_dir) is skipped whenever "
               "anything between them raises. The failure of one simulation "
               "therefore leaves the whole interpreter parked inside that "
               "simulation's directory, and every later relative path in "
               "the session resolves somewhere unexpected. See the "
               "consolidated findings issue.",
    )
    def test_a_failure_inside_the_run_should_still_restore_the_directory(
        self, template, ifmc, solver, cwd
    ):
        with pytest.raises(AttributeError):
            self._run(template, topography="TOPO DATA", canopy="CANOPY DATA")
        assert os.getcwd() == str(cwd)

    def test_a_failure_inside_the_run_currently_leaves_the_directory_changed(
        self, template, ifmc, solver, cwd, lsm
    ):
        """Characterisation of B278."""
        try:
            with pytest.raises(AttributeError):
                self._run(template, topography="TOPO DATA", canopy="CANOPY DATA")
            assert os.getcwd() == os.path.join(
                os.path.abspath(lsm.filesDirectory), "sim"
            )
        finally:
            os.chdir(str(cwd))


# ---------------------------------------------------------------------------
# _toNetcdf
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestToNetcdf:
    @pytest.fixture()
    def outputFiles(self, tmp_path):
        """Two OUTD files, written newest-first so ordering is observable."""
        directory = tmp_path / "machsan"
        directory.mkdir()
        for moment in (120.0, 60.0):
            (directory / f"OUTD3d03_3_{moment}").write_text(
                "0 0 0 6.0\n0 0 1 3.0\n1 0 0 12.0\n1 0 1 0.0\n"
            )
        return str(directory / "OUTD3d03_3_")

    def test_it_yields_one_frame_per_file_plus_a_leading_zero_frame(
        self, template, outputFiles
    ):
        frames = list(
            template()._toNetcdf(basefiles=outputFiles, datetimeFormat="seconds")
        )
        assert len(frames) == 3

    def test_the_leading_frame_is_all_zeros(self, template, outputFiles):
        frames = list(
            template()._toNetcdf(basefiles=outputFiles, datetimeFormat="seconds")
        )
        assert float(frames[0]["Dosage"].max()) == 0.0
        assert float(frames[0]["Dosage"].min()) == 0.0

    def test_the_leading_frame_sits_at_time_zero(self, template, outputFiles):
        frames = list(
            template()._toNetcdf(basefiles=outputFiles, datetimeFormat="seconds")
        )
        assert list(frames[0].datetime.values) == [0]

    def test_the_leading_frame_can_be_switched_off(self, template, outputFiles):
        frames = list(
            template()._toNetcdf(
                basefiles=outputFiles, addzero=False, datetimeFormat="seconds"
            )
        )
        assert len(frames) == 2
        assert float(frames[0]["Dosage"].max()) > 0

    def test_the_frames_are_ordered_by_the_time_in_the_filename(
        self, template, outputFiles
    ):
        frames = list(
            template()._toNetcdf(
                basefiles=outputFiles, addzero=False, datetimeFormat="seconds"
            )
        )
        assert [float(frame.datetime.values[0]) for frame in frames] == [60.0, 120.0]

    def test_every_frame_is_indexed_by_time_and_the_three_axes(
        self, template, outputFiles
    ):
        frames = list(
            template()._toNetcdf(basefiles=outputFiles, datetimeFormat="seconds")
        )
        assert frames[0]["Dosage"].dims == ("datetime", "x", "y", "z")

    def test_the_grid_comes_from_the_file_contents(self, template, outputFiles):
        frames = list(
            template()._toNetcdf(basefiles=outputFiles, datetimeFormat="seconds")
        )
        frame = frames[1]
        assert list(frame.x.values) == [0, 1]
        assert list(frame.y.values) == [0]
        assert list(frame.z.values) == [0, 1]

    def test_the_timestamp_format_anchors_the_series_at_the_documented_epoch(
        self, template, outputFiles
    ):
        frames = list(
            template()._toNetcdf(basefiles=outputFiles, datetimeFormat="timestamp")
        )
        assert pandas.Timestamp(frames[0].datetime.values[0]) == pandas.Timestamp(
            "2016-01-01 12:00"
        )
        assert pandas.Timestamp(frames[1].datetime.values[0]) == pandas.Timestamp(
            "2016-01-01 12:01"
        )

    @pytest.mark.xfail(
        strict=True,
        reason="B275: the docstring says 'The dosage are converted to "
               "s/m**3 instead of min/m**3', but the code multiplies by "
               "(1*ureg.s / ureg.m**3).m_as(ureg.minute / ureg.m**3), which "
               "is 1/60 -- the seconds-to-minutes factor. Whatever the file "
               "unit is, dividing by 60 cannot be the minutes-to-seconds "
               "conversion the docstring promises; a factor of 60 is. Code "
               "and docstring disagree by a factor of 3600. See the "
               "consolidated findings issue.",
    )
    def test_the_dosage_should_be_scaled_to_seconds_per_cubic_metre(
        self, template, outputFiles
    ):
        frames = list(
            template()._toNetcdf(
                basefiles=outputFiles, addzero=False, datetimeFormat="seconds"
            )
        )
        assert float(frames[0]["Dosage"].sel(x=0, y=0, z=0)) == pytest.approx(
            6.0 * 60
        )

    def test_the_dosage_is_currently_divided_by_sixty(self, template, outputFiles):
        """Characterisation of B275."""
        frames = list(
            template()._toNetcdf(
                basefiles=outputFiles, addzero=False, datetimeFormat="seconds"
            )
        )
        assert float(frames[0]["Dosage"].sel(x=0, y=0, z=0)) == pytest.approx(
            6.0 / 60
        )
        assert float(frames[0]["Dosage"].sel(x=1, y=0, z=0)) == pytest.approx(
            12.0 / 60
        )

    def test_no_files_means_no_frames(self, template, tmp_path):
        empty = str(tmp_path / "nothingHere")
        assert list(template()._toNetcdf(basefiles=empty)) == []


# ---------------------------------------------------------------------------
# _getSimulationsList / getSimulationsTable
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestSimulationsListing:
    def _register(self, lsm, simulationName, **extraParams):
        params = dict(BASE_PARAMS)
        params.update(extraParams)
        return lsm.addSimulationsDocument(
            type="LSM_run",
            resource=f"resource_{simulationName}",
            dataFormat="string",
            desc=dict(
                version=TEMPLATE_VERSION,
                templateName=TEMPLATE_NAME,
                simulationName=simulationName,
                params=params,
            ),
        )

    def test_one_simulation_gives_one_frame(self, template, lsm):
        self._register(lsm, "sim1")
        assert len(template()._getSimulationsList()) == 1

    def test_the_frame_carries_the_document_id(self, template, lsm):
        document = self._register(lsm, "sim1")
        frame = template()._getSimulationsList()[0]
        assert frame["id"][0] == document.id

    def test_the_frame_carries_the_simulation_metadata(self, template, lsm):
        self._register(lsm, "sim1")
        frame = template()._getSimulationsList()[0]
        assert frame["simulationName"][0] == "sim1"
        assert frame["templateName"][0] == TEMPLATE_NAME
        assert frame["version"][0] == TEMPLATE_VERSION

    def test_the_parameters_are_prefixed_so_they_cannot_clash(self, template, lsm):
        self._register(lsm, "sim1")
        frame = template()._getSimulationsList()[0]
        assert frame["params__TopoXmax"][0] == 100
        assert "TopoXmax" not in frame.columns

    def test_a_simulation_with_other_parameters_is_not_listed(self, template, lsm):
        self._register(lsm, "sim1", TopoXmax=100)
        self._register(lsm, "other", TopoXmax=999)
        frames = template()._getSimulationsList()
        assert len(frames) == 1
        assert frames[0]["simulationName"][0] == "sim1"

    def test_a_query_narrows_the_listing(self, template, lsm):
        self._register(lsm, "sim1", TopoXmax=100)
        self._register(lsm, "other", TopoXmax=999)
        frames = template()._getSimulationsList(TopoXmax=999)
        assert len(frames) == 1
        assert frames[0]["simulationName"][0] == "other"

    def test_an_empty_project_gives_no_frames(self, template, lsm):
        assert template()._getSimulationsList() == []

    @pytest.mark.xfail(
        strict=True,
        reason="B276: _getSimulationsList builds one desc frame and one "
               "params frame per document, then combines them with "
               "`product(desc_df_list, params_df_list)` instead of pairing "
               "them. n simulations therefore yield n^2 rows, and every "
               "off-diagonal row staples one simulation's id, name and "
               "version onto a different simulation's parameters. "
               "getSimulationsTable inherits the whole defect, and "
               "LSMToolkit.getSimulationsList / getTemplatesTable repeat the "
               "same construction. See the consolidated findings issue.",
    )
    def test_two_simulations_should_give_two_frames(self, template, lsm):
        self._register(lsm, "sim1", extra=1)
        self._register(lsm, "sim2", extra=2)
        assert len(template()._getSimulationsList()) == 2

    def test_two_simulations_currently_give_four_frames(self, template, lsm):
        """Characterisation of B276."""
        self._register(lsm, "sim1", extra=1)
        self._register(lsm, "sim2", extra=2)
        assert len(template()._getSimulationsList()) == 2**2

    def test_the_cross_join_currently_mixes_one_run_with_anothers_parameters(
        self, template, lsm
    ):
        """Characterisation of B276."""
        self._register(lsm, "sim1", extra=1)
        self._register(lsm, "sim2", extra=2)
        pairs = {
            (frame["simulationName"][0], frame["params__extra"][0])
            for frame in template()._getSimulationsList()
        }
        assert pairs == {("sim1", 1), ("sim1", 2), ("sim2", 1), ("sim2", 2)}

    def test_three_simulations_currently_give_nine_frames(self, template, lsm):
        """Characterisation of B276: the growth is quadratic."""
        for index in range(3):
            self._register(lsm, f"sim{index}", extra=index)
        assert len(template()._getSimulationsList()) == 3**2


@pytest.mark.unit
class TestSimulationsTable:
    def _register(self, lsm, simulationName, **extraParams):
        params = dict(BASE_PARAMS)
        params.update(extraParams)
        return lsm.addSimulationsDocument(
            type="LSM_run",
            resource=f"resource_{simulationName}",
            dataFormat="string",
            desc=dict(
                version=TEMPLATE_VERSION,
                templateName=TEMPLATE_NAME,
                simulationName=simulationName,
                params=params,
            ),
        )

    def test_an_empty_project_is_reported_as_a_missing_file(self, template, lsm):
        """pandas.concat of nothing raises ValueError, remapped by the method."""
        with pytest.raises(FileNotFoundError, match="No simulations found"):
            template().getSimulationsTable()

    def test_one_simulation_gives_one_row_indexed_by_its_id(self, template, lsm):
        document = self._register(lsm, "sim1")
        table = template().getSimulationsTable()
        assert len(table) == 1
        assert list(table.index) == [document.id]

    def test_the_row_holds_the_metadata_and_the_parameters(self, template, lsm):
        self._register(lsm, "sim1")
        table = template().getSimulationsTable()
        assert table["simulationName"].iloc[0] == "sim1"
        assert table["params__TopoXmax"].iloc[0] == 100

    def test_the_id_column_is_moved_into_the_index(self, template, lsm):
        self._register(lsm, "sim1")
        assert "id" not in template().getSimulationsTable().columns

    def test_a_query_narrows_the_table(self, template, lsm):
        self._register(lsm, "sim1", TopoXmax=100)
        self._register(lsm, "other", TopoXmax=999)
        table = template().getSimulationsTable(TopoXmax=999)
        assert len(table) == 1
        assert table["simulationName"].iloc[0] == "other"

    @pytest.mark.xfail(
        strict=True,
        reason="B276: getSimulationsTable inherits _getSimulationsList's "
               "product() cross-join, so a project with n matching "
               "simulations is reported as n^2 rows with duplicated index "
               "values and parameters taken from the wrong run. See the "
               "consolidated findings issue.",
    )
    def test_two_simulations_should_give_two_rows(self, template, lsm):
        self._register(lsm, "sim1", extra=1)
        self._register(lsm, "sim2", extra=2)
        assert len(template().getSimulationsTable()) == 2

    def test_two_simulations_currently_give_four_rows(self, template, lsm):
        """Characterisation of B276."""
        self._register(lsm, "sim1", extra=1)
        self._register(lsm, "sim2", extra=2)
        assert len(template().getSimulationsTable()) == 4

    def test_each_id_currently_appears_twice_in_the_index(self, template, lsm):
        """Characterisation of B276."""
        first = self._register(lsm, "sim1", extra=1)
        second = self._register(lsm, "sim2", extra=2)
        index = list(template().getSimulationsTable().index)
        assert index.count(first.id) == 2
        assert index.count(second.id) == 2
