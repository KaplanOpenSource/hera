"""openFoam/lagrangian/LSM/toolkit.py: OFLSMToolkit and its Analysis layer,
beyond the four properties already pinned in
``test_openfoam_lagrangian_lsm_toolkit.py`` (which also pins B103, the
constructor that can never run).

Two construction routes are used here, both deliberate:

* ``lsm`` -- a real, DB-backed toolkit from ``unit_toolkit_factory``,
  obtained by patching the missing constant ``ToolkitHome.GIS_TOPOGRAPHY``
  onto the *class* for the duration of a test, i.e. by simulating the fix
  for B103. That is the only way to exercise ``__init__``, the topography
  wiring and everything that talks to the datalayer; B103 itself stays
  pinned, unpatched, in the sibling file.
* ``OFLSMToolkit.__new__(OFLSMToolkit)`` -- the technique the sibling file
  introduced -- for the parsers and writers that touch nothing but the
  filesystem, so that no unrelated failure can be blamed on them.

Everything OpenFOAM-shaped is a file under ``tmp_path``; evtk is the
conftest stub, so the VTK tests assert which datasets crossed the seam
rather than any value the stub invented.

Deliberately not covered:

* ``makeDistanceFromGround`` and ``makeUstar``'s uncached path -- both need
  ``getCellDataAndGroundData`` to read a real ``polyMesh`` plus a
  ``ground`` patch out of an OpenFOAM case, and then
  ``topography.analysis.addHeight`` on top of it. ``makeUstar``'s cached
  path *is* covered below, which is what exercises the shared
  ``_loadOrComputeCellData`` / ``_parseOFFieldHeader`` /
  ``_writeScalarField`` chain.
* ``_interpolateNearGroundValues`` -- it is a
  ``coordinateHandler.regularizeTimeSteps`` sandwich over a real
  2-D grid; testing it means testing that helper, which belongs with the
  helper.
* ``__init__``'s ``casePath`` argument -- the unit fixtures build a toolkit
  without extra constructor keywords, and the default (the cwd) is pinned
  instead.
* ``calcConcentrationPointWise``'s dask path -- the identical body on the
  abstract Lagrangian solver is covered by
  ``test_openfoam_lagrangian_analysis.py``.

Bugs pinned below (each with an xfail for the intended behaviour and a
passing characterisation of what happens today):

* B201 ``_extractFile`` / ``_readRecord`` use ``self.logger``, which no
  toolkit has.
* B202 ``loadData``'s saveMode guard rejects a documented mode and its
  error message crashes.
* B203 ``loadData``'s non-parallel branch finds no time steps.
* B204 ``loadData``'s DB branch treats a QuerySet as a document.
* B205 ``loadData`` with ONLYFILE_REPLACE returns an unbound name.
* B206 ``makeSource`` writes the wrong object name into the header.
* B207 ``createRootCaseMeshLink`` references an undefined name.
* B208 ``calcConcentrationTimeStepFullMesh`` sets an attribute xarray
  refuses.
* B209 ``calcDocumentConcentrationPointWise`` never returns.
* B210 ``getConcentrationField(returnFirst=False)`` indexes an empty list.
* B211 ``_extractVelocityField`` cannot read a field file with a real
  boundary section.
"""
import os

import pandas
import pytest

from hera import toolkit as toolkitModule
from hera.simulations.openFoam.lagrangian.LSM import toolkit as lsmModule
from hera.simulations.openFoam.lagrangian.LSM.toolkit import Analysis, OFLSMToolkit

_EXTENTS = dict(xmin=0, xmax=2, ymin=0, ymax=2, zmin=0, zmax=2)


@pytest.fixture()
def lsm(unit_toolkit_factory, monkeypatch):
    """A real OF_LSM toolkit, with B103's missing constant supplied.

    The constant is patched onto the ToolkitHome *class* (never onto the
    toolkitHome instance, which monkeypatch would then "restore" as a
    permanent instance attribute) and disappears again after the test.
    """
    from hera import ToolkitHome, toolkitHome

    monkeypatch.setattr(ToolkitHome, "GIS_TOPOGRAPHY", toolkitHome.GIS_RASTER_TOPOGRAPHY,
                        raising=False)
    return unit_toolkit_factory("OF_LSM")


@pytest.fixture()
def bare():
    """The toolkit without any constructor at all -- see the module docstring."""
    toolkit = OFLSMToolkit.__new__(OFLSMToolkit)
    toolkit._cloudName = "kinematicCloud"
    return toolkit


def _lagrangianFile(path, lines):
    """An OpenFOAM lagrangian field file: _extractFile skips 20 header lines
    and 4 footer lines."""
    header = [f"// header line {index}" for index in range(20)]
    footer = [")", ";", "// ***", "//"]
    path.write_text("\n".join(header + list(lines) + footer) + "\n")


def _cloudDirectory(casePath, timeName="100", cloudName="kinematicCloud"):
    directory = casePath / timeName / "lagrangian" / cloudName
    directory.mkdir(parents=True)
    _lagrangianFile(directory / "globalSigmaPositions", ["(1 2 3)"])
    _lagrangianFile(directory / "origId", ["5"])
    _lagrangianFile(directory / "origProcId", ["0"])
    _lagrangianFile(directory / "globalPositions", ["(10 20 30)"])
    return directory


def _particleData():
    return pandas.DataFrame(dict(time=[1.0, 1.0, 2.0],
                                 globalX=[0.5, 1.5, 0.5],
                                 globalY=[0.5, 0.5, 0.5],
                                 globalZ=[0.5, 0.5, 0.5],
                                 mass=[1.0, 2.0, 3.0]))


class _Recorder:
    def __init__(self, returnValue=None):
        self.calls = []
        self.returnValue = returnValue

    def __call__(self, *args, **kwargs):
        self.calls.append((args, kwargs))
        return self.returnValue

    @property
    def lastArgs(self):
        return self.calls[-1][0]

    @property
    def lastKwargs(self):
        return self.calls[-1][1]


# ---------------------------------------------------------------------------
# Construction (with B103's constant supplied)
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestConstruction:
    def test_the_case_path_defaults_to_the_current_directory(self, lsm):
        assert lsm.casePath == os.getcwd()

    def test_the_cloud_name_defaults_to_the_kinematic_cloud(self, lsm):
        assert lsm.cloudName == "kinematicCloud"

    def test_a_case_is_assumed_to_be_single_processor(self, lsm):
        assert lsm.parallelCase is False

    def test_the_analysis_layer_points_back_at_the_toolkit(self, lsm):
        assert isinstance(lsm.analysis, Analysis)
        assert lsm.analysis.datalayer is lsm

    def test_a_sources_factory_is_built(self, lsm):
        from hera.simulations.openFoam.lagrangian.LSM.sourcesFactoryTool import sourcesFactoryTool

        assert isinstance(lsm.sourcesFactory, sourcesFactoryTool)

    def test_the_topography_toolkit_shares_the_projects_files_directory(self, lsm, unit_files_directory):
        assert lsm.topography.projectName == lsm.projectName
        assert lsm.topography.filesDirectory == unit_files_directory

    def test_the_toolkit_name_is_registered_on_the_project(self, lsm):
        assert lsm.toolkitName == "OF_LSM"
        assert lsm.doctype == "LSMRuns"


# ---------------------------------------------------------------------------
# _extractFile / _readRecord
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestExtractFile:
    def test_a_vector_field_becomes_one_column_per_component(self, bare, tmp_path):
        path = tmp_path / "globalPositions"
        _lagrangianFile(path, ["(1 2 3)", "(4 5 6)"])
        frame = bare._extractFile(str(path), ["x", "y", "z"])
        assert list(frame.columns) == ["x", "y", "z"]
        assert frame.to_dict("records") == [dict(x=1.0, y=2.0, z=3.0), dict(x=4.0, y=5.0, z=6.0)]

    def test_a_scalar_field_becomes_a_single_column(self, bare, tmp_path):
        path = tmp_path / "origId"
        _lagrangianFile(path, ["7", "8"])
        frame = bare._extractFile(str(path), ["id"], vector=False)
        assert frame["id"].tolist() == [7.0, 8.0]

    def test_every_value_is_a_float(self, bare, tmp_path):
        path = tmp_path / "origProcId"
        _lagrangianFile(path, ["0", "1"])
        assert bare._extractFile(str(path), ["procId"], vector=False)["procId"].dtype == float

    def test_the_parentheses_of_a_record_are_stripped(self, bare, tmp_path):
        path = tmp_path / "U"
        _lagrangianFile(path, ["(0.5 -1.5 2.5)"])
        assert bare._extractFile(str(path), ["U_x", "U_y", "U_z"]).iloc[0].tolist() == [0.5, -1.5, 2.5]

    def test_the_specialised_parser_is_unreachable(self, bare, tmp_path):
        """Characterisation of B201.

        A file the csv reader rejects is meant to fall through to the
        hand-written OpenFOAM parser, but the fall-through starts with
        self.logger.execute(...) -- and no toolkit has a `logger` attribute,
        so that branch can only ever raise.
        """
        path = tmp_path / "notACsv"
        _lagrangianFile(path, ["3{(1 2 3)}"])
        with pytest.raises(AttributeError, match="logger"):
            bare._extractFile(str(path), ["x", "y", "z"])


@pytest.mark.unit
class TestReadRecord:
    @pytest.mark.xfail(
        strict=True,
        reason="B201: _readRecord ends with self.logger.debug(...) (and "
               "its except branch starts with one), but neither "
               "abstractToolkit nor Project defines a `logger` attribute -- "
               "they all build a local logger with "
               "get_classMethod_logger(self, ...). Every call therefore ends "
               "in AttributeError, whether or not the time step holds data, "
               "which also makes loadData unusable with the real reader. "
               "See the consolidated findings issue.",
    )
    def test_a_time_step_should_be_read_into_a_frame(self, bare, tmp_path):
        _cloudDirectory(tmp_path)
        frame = bare._readRecord("100", casePath=str(tmp_path))
        assert frame["globalX"].tolist() == [10.0]
        assert frame["time"].tolist() == [100.0]

    def test_reading_a_time_step_currently_raises(self, bare, tmp_path):
        """Characterisation of B201."""
        _cloudDirectory(tmp_path)
        with pytest.raises(AttributeError, match="logger"):
            bare._readRecord("100", casePath=str(tmp_path))

    def test_an_empty_time_step_also_raises(self, bare, tmp_path):
        """Characterisation of B201: the except branch logs too."""
        with pytest.raises(AttributeError, match="logger"):
            bare._readRecord("100", casePath=str(tmp_path))


# ---------------------------------------------------------------------------
# loadData
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestLoadData:
    @pytest.fixture()
    def readSeam(self, monkeypatch):
        """Replace _readRecord: the real one cannot run (B201)."""
        recorder = _Recorder(returnValue=pandas.DataFrame(dict(time=[100.0], globalX=[1.0],
                                                               mass=[2.0])))
        monkeypatch.setattr(OFLSMToolkit, "_readRecord",
                            lambda self, timeName, **kwargs: recorder(timeName, **kwargs))
        return recorder

    @pytest.fixture()
    def parallelCase(self, lsm, tmp_path):
        case = tmp_path / "case"
        (case / "processor0" / "100").mkdir(parents=True)
        (case / "processor0" / "200").mkdir(parents=True)
        lsm.casePath = str(case)
        lsm.parallelCase = True
        return case

    def test_the_data_of_a_decomposed_case_is_returned_without_the_db(self, lsm, parallelCase, readSeam):
        from hera.datalayer import nonDBMetadataFrame

        document = lsm.loadData()
        assert isinstance(document, nonDBMetadataFrame)
        assert document.getData().compute()["globalX"].tolist() == [1.0, 1.0]

    def test_every_processor_time_directory_is_read(self, lsm, parallelCase, readSeam):
        lsm.loadData().getData().compute()
        assert {call[0][0] for call in readSeam.calls} == {os.path.join("processor0", "100"),
                                                           os.path.join("processor0", "200")}

    def test_the_reader_options_are_forwarded(self, lsm, parallelCase, readSeam):
        lsm.loadData(withVelocity=True, withMass=True, withReleaseTimes=True).getData().compute()
        assert readSeam.lastKwargs == dict(casePath=str(parallelCase), withVelocity=True,
                                           withReleaseTimes=True, withMass=True)

    def test_an_explicit_case_path_wins_over_the_toolkits(self, lsm, parallelCase, readSeam, tmp_path):
        other = tmp_path / "other"
        (other / "processor0" / "100").mkdir(parents=True)
        lsm.loadData(casePath=str(other)).getData().compute()
        assert readSeam.lastKwargs["casePath"] == str(other)

    def test_a_decomposed_case_without_processors_is_reported(self, lsm, tmp_path, readSeam):
        empty = tmp_path / "empty"
        empty.mkdir()
        lsm.casePath = str(empty)
        lsm.parallelCase = True
        with pytest.raises(ValueError, match="processor"):
            lsm.loadData()

    def test_the_parquet_is_named_after_the_cloud(self, lsm, parallelCase, readSeam):
        lsm.cloudName = "myCloud"
        # The file is written before the method trips over B205.
        with pytest.raises(UnboundLocalError):
            lsm.loadData(saveMode=toolkitModule.TOOLKIT_SAVEMODE_ONLYFILE_REPLACE)
        assert os.path.exists(os.path.join(lsm.FilesDirectory, "myCloudData.parquet"))

    # -- B202: the saveMode guard -----------------------------------------

    @pytest.mark.xfail(
        strict=True,
        reason="B202: an unknown saveMode is meant to raise ValueError "
               "listing the valid modes, but the message is built with "
               "','.join(...) over a list that contains "
               "TOOLKIT_SAVEMODE_NOSAVE, which is None -- so the guard "
               "itself dies with TypeError before the ValueError is raised. "
               "The same list omits TOOLKIT_SAVEMODE_ONLYFILE (and repeats "
               "FILEANDDB twice), so that documented mode is rejected even "
               "though the body below handles it explicitly. "
               "See the consolidated findings issue.",
    )
    def test_an_unknown_savemode_should_be_reported_as_a_valueerror(self, lsm):
        with pytest.raises(ValueError, match="saveMode must be one of"):
            lsm.loadData(saveMode="noSuchMode")

    def test_an_unknown_savemode_currently_raises_typeerror(self, lsm):
        """Characterisation of B202."""
        with pytest.raises(TypeError, match="expected str instance"):
            lsm.loadData(saveMode="noSuchMode")

    def test_the_onlyfile_mode_is_rejected_by_the_guard(self, lsm, parallelCase, readSeam):
        """Characterisation of B202: a documented, handled mode is refused
        because it is missing from the guard's list of valid modes."""
        with pytest.raises(TypeError, match="expected str instance"):
            lsm.loadData(saveMode=toolkitModule.TOOLKIT_SAVEMODE_ONLYFILE)

    # -- B203: the non-parallel branch ------------------------------------

    @pytest.mark.xfail(
        strict=True,
        reason="B203: the single-processor branch filters the case's "
               "entries with os.path.isdir(x) on the bare name instead of "
               "os.path.join(finalCasePath, x), so it tests the cwd rather "
               "than the case. Every time directory fails the test, the "
               "delayed list comes out empty and dask raises 'Must supply at "
               "least one delayed object'. Same root cause as B83 in "
               "openFoam/toolkit.py's getTimeList. "
               "See the consolidated findings issue.",
    )
    def test_a_single_processor_case_should_be_read(self, lsm, tmp_path, readSeam):
        case = tmp_path / "case"
        (case / "100").mkdir(parents=True)
        lsm.casePath = str(case)
        lsm.loadData().getData().compute()
        assert readSeam.calls != []

    def test_a_single_processor_case_currently_finds_no_time_steps(self, lsm, tmp_path, readSeam):
        """Characterisation of B203."""
        case = tmp_path / "case"
        (case / "100").mkdir(parents=True)
        lsm.casePath = str(case)
        with pytest.raises(TypeError, match="at least one delayed object"):
            lsm.loadData()

    # -- B204: the DB branch ----------------------------------------------

    @pytest.mark.xfail(
        strict=True,
        reason="B204: the FILEANDDB branch does `doc = "
               "self.getSimulationsDocuments(...)`, which returns a QuerySet "
               "-- never None -- and then tests `if doc is not None`. So the "
               "'already in the DB' error fires even for the very first "
               "save, the `elif doc is None` branch that would call "
               "addSimulationsDocument is dead code, and no LSMRuns document "
               "is ever registered. FILEANDDB_REPLACE reaches the else "
               "branch instead and dies on QuerySet.resource/save. "
               "See the consolidated findings issue.",
    )
    def test_saving_to_the_db_should_register_a_document(self, lsm, parallelCase, readSeam):
        lsm.loadData(saveMode=toolkitModule.TOOLKIT_SAVEMODE_FILEANDDB)
        assert len(lsm.getSimulationsDocuments(type="LSMRuns")) == 1

    def test_saving_to_the_db_currently_claims_the_data_is_already_there(self, lsm, parallelCase,
                                                                        readSeam):
        """Characterisation of B204."""
        with pytest.raises(FileExistsError, match="already in the DB"):
            lsm.loadData(saveMode=toolkitModule.TOOLKIT_SAVEMODE_FILEANDDB)
        assert len(lsm.getSimulationsDocuments(type="LSMRuns")) == 0
        # The parquet was written before the DB step, so the data is on disk.
        assert os.path.exists(os.path.join(lsm.FilesDirectory, "kinematicCloudData.parquet"))

    def test_replacing_in_the_db_currently_dies_on_the_queryset(self, lsm, parallelCase, readSeam):
        """Characterisation of B204."""
        with pytest.raises(AttributeError, match="QuerySet"):
            lsm.loadData(saveMode=toolkitModule.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE)

    # -- B205: ONLYFILE_REPLACE -------------------------------------------

    @pytest.mark.xfail(
        strict=True,
        reason="B205: with ONLYFILE_REPLACE the parquet is written but no "
               "document is built -- `doc` is only bound in the DB branches "
               "and in the else that handles NOSAVE -- so `return doc` "
               "raises UnboundLocalError after the file has already been "
               "written. See the consolidated findings issue.",
    )
    def test_saving_only_to_a_file_should_return_something(self, lsm, parallelCase, readSeam):
        assert lsm.loadData(saveMode=toolkitModule.TOOLKIT_SAVEMODE_ONLYFILE_REPLACE) is not None

    def test_saving_only_to_a_file_writes_it_and_then_raises(self, lsm, parallelCase, readSeam):
        """Characterisation of B205."""
        with pytest.raises(UnboundLocalError, match="doc"):
            lsm.loadData(saveMode=toolkitModule.TOOLKIT_SAVEMODE_ONLYFILE_REPLACE)
        assert os.path.exists(os.path.join(lsm.FilesDirectory, "kinematicCloudData.parquet"))


# ---------------------------------------------------------------------------
# makeSource
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestMakeSource:
    @pytest.fixture()
    def case(self, bare, tmp_path):
        from hera.simulations.openFoam.lagrangian.LSM.sourcesFactoryTool import sourcesFactoryTool

        bare._casePath = str(tmp_path)
        bare._sourcesFactory = sourcesFactoryTool()
        (tmp_path / "constant").mkdir()
        return tmp_path

    def test_the_positions_are_written_into_the_constant_directory(self, bare, case):
        bare.makeSource(x=1, y=2, z=3, nParticles=2)
        assert (case / "constant" / "kinematicCloudPositions").exists()

    def test_the_file_name_is_honoured(self, bare, case):
        bare.makeSource(x=1, y=2, z=3, nParticles=1, fileName="mySource")
        assert (case / "constant" / "mySource").exists()

    def test_every_particle_of_a_point_source_is_written_at_the_point(self, bare, case):
        bare.makeSource(x=1.0, y=2.0, z=3.0, nParticles=3)
        content = (case / "constant" / "kinematicCloudPositions").read_text()
        assert content.count("(1.0 2.0 3.0)") == 3

    def test_the_particle_count_precedes_the_list(self, bare, case):
        bare.makeSource(x=0, y=0, z=0, nParticles=4)
        lines = [line.strip() for line in
                 (case / "constant" / "kinematicCloudPositions").read_text().splitlines()]
        assert lines[lines.index("(") - 1] == "4"
        assert lines[-1] == ")"

    def test_the_source_geometry_is_delegated_to_the_factory(self, bare, case):
        bare.makeSource(x=0, y=0, z=0, nParticles=50, type="Sphere", radius=1.0)
        content = (case / "constant" / "kinematicCloudPositions").read_text()
        assert content.count("(") == 51  # 50 records plus the opening bracket

    def test_an_unknown_source_type_is_reported(self, bare, case):
        with pytest.raises(ValueError, match="The type must be"):
            bare.makeSource(x=0, y=0, z=0, nParticles=1, type="noSuchGeometry")

    @pytest.mark.xfail(
        strict=True,
        reason="B206: the FoamFile header is a single hard-coded string "
               "whose `object` entry is always kinematicCloudPositions, even "
               "when fileName says otherwise -- so a source written under any "
               "other name declares the wrong object to OpenFOAM. "
               "See the consolidated findings issue.",
    )
    def test_the_header_should_name_the_file_it_is_in(self, bare, case):
        bare.makeSource(x=0, y=0, z=0, nParticles=1, fileName="mySource")
        assert "object      mySource;" in (case / "constant" / "mySource").read_text()

    def test_the_header_currently_always_names_the_kinematic_cloud(self, bare, case):
        """Characterisation of B206."""
        bare.makeSource(x=0, y=0, z=0, nParticles=1, fileName="mySource")
        assert "object      kinematicCloudPositions;" in (case / "constant" / "mySource").read_text()


# ---------------------------------------------------------------------------
# The OpenFOAM field-file helpers
# ---------------------------------------------------------------------------

def _writeVectorField(path, boundary=("    ground { type zeroGradient; }",)):
    path.write_text("\n".join([
        "FoamFile", "{", "    class       volVectorField;", "}",
        "internalField   nonuniform List<vector>",
        "3", "(", "(1 0 0)", "(0 2 0)", "(0 0 3)", ")", ";",
        "boundaryField", "{", *boundary, "}",
    ]) + "\n")


@pytest.mark.unit
class TestFieldFileHelpers:
    def test_the_header_parser_finds_the_values_and_the_count(self, tmp_path):
        path = tmp_path / "U"
        _writeVectorField(path)
        lines, cellStart, nCells = OFLSMToolkit._parseOFFieldHeader(str(path))
        assert nCells == 3
        assert lines[cellStart].strip() == "(1 0 0)"

    def test_a_file_without_an_internal_field_yields_zeroes(self, tmp_path):
        path = tmp_path / "empty"
        path.write_text("FoamFile\n{\n}\n")
        assert OFLSMToolkit._parseOFFieldHeader(str(path))[1:] == (0, 0)

    def test_the_boundary_parser_finds_the_boundary_line(self, tmp_path):
        path = tmp_path / "U"
        _writeVectorField(path)
        lines, boundaryStart = OFLSMToolkit._parseBoundaryStart(str(path))
        assert lines[boundaryStart].strip() == "boundaryField"

    def test_a_file_without_a_boundary_field_points_past_the_end(self, tmp_path):
        path = tmp_path / "noBoundary"
        path.write_text("a\nb\n")
        lines, boundaryStart = OFLSMToolkit._parseBoundaryStart(str(path))
        assert boundaryStart == len(lines)

    def test_the_scalar_writer_rewrites_the_header_and_the_values(self, tmp_path):
        source = tmp_path / "U"
        _writeVectorField(source)
        lines, cellStart, nCells = OFLSMToolkit._parseOFFieldHeader(str(source))
        boundaryLines, boundaryStart = OFLSMToolkit._parseBoundaryStart(str(source))
        output = tmp_path / "ustar"
        OFLSMToolkit._writeScalarField(lines, cellStart, nCells, boundaryLines, boundaryStart,
                                       pandas.Series([0.2, 0.3, 0.4, 0.9]), str(output))
        written = output.read_text()
        assert "volScalarField" in written and "volVectorField" not in written
        assert "List<scalar>" in written
        # exactly nCells values, then the closing of the list, then the
        # boundary section copied from the template.
        assert written.count("\n0.") == 3 and "0.9" not in written
        assert written.strip().endswith("}")
        assert "ground { type zeroGradient; }" in written


@pytest.mark.unit
class TestExtractVelocityField:
    @pytest.fixture()
    def case(self, bare, tmp_path):
        """A case whose boundary section has at most three tokens per line --
        the only shape the reader survives (see B211 below)."""
        case = tmp_path / "case"
        (case / "50").mkdir(parents=True)
        _writeVectorField(case / "50" / "U", boundary=("    ground", "    {", "    }"))
        bare._casePath = str(case)
        return case

    def test_the_velocity_magnitude_is_joined_onto_the_cell_data(self, bare, case):
        _, cellStart, nCells = OFLSMToolkit._parseOFFieldHeader(str(case / "50" / "U"))
        cellData = pandas.DataFrame(dict(x=[0.0, 1.0, 2.0], y=[0.0, 0.0, 0.0],
                                         height=[1.0, 1.5, 2.0]))
        joined = bare._extractVelocityField(cellData, cellStart, nCells, 50)
        assert joined["U"].tolist() == [1.0, 2.0, 3.0]
        assert list(joined.columns) == ["x", "y", "height", "u", "v", "w", "U"]

    def test_only_the_internal_cells_are_kept(self, bare, case):
        _, cellStart, _ = OFLSMToolkit._parseOFFieldHeader(str(case / "50" / "U"))
        cellData = pandas.DataFrame(dict(x=[0.0, 1.0], y=[0.0, 0.0], height=[1.0, 1.5]))
        assert len(bare._extractVelocityField(cellData, cellStart, 2, 50)) == 2

    @pytest.mark.xfail(
        strict=True,
        reason="B211: the U field is read with pandas.read_csv over the "
               "WHOLE file (skipfooter=0, no nrows) into exactly three "
               "columns, and only trimmed to nCells afterwards with "
               ".head(nCells). Any line below the values that has more than "
               "three whitespace-separated tokens -- i.e. an ordinary "
               "boundary entry such as 'value uniform (0 0 0);' -- makes "
               "pandas raise ParserError, so the reader cannot handle a "
               "real OpenFOAM field file. nrows=nCells would be the fix. "
               "See the consolidated findings issue.",
    )
    def test_a_realistic_boundary_section_should_be_tolerated(self, bare, tmp_path):
        case = tmp_path / "realistic"
        (case / "50").mkdir(parents=True)
        _writeVectorField(case / "50" / "U",
                          boundary=("    inlet", "    {", "        type fixedValue;",
                                    "        value uniform (0 0 0);", "    }"))
        bare._casePath = str(case)
        _, cellStart, nCells = OFLSMToolkit._parseOFFieldHeader(str(case / "50" / "U"))
        cellData = pandas.DataFrame(dict(x=[0.0, 1.0, 2.0], y=[0.0, 0.0, 0.0],
                                         height=[1.0, 1.5, 2.0]))
        assert bare._extractVelocityField(cellData, cellStart, nCells, 50)["U"].tolist() == [1.0, 2.0, 3.0]

    def test_a_realistic_boundary_section_currently_raises(self, bare, tmp_path):
        """Characterisation of B211."""
        case = tmp_path / "realistic"
        (case / "50").mkdir(parents=True)
        _writeVectorField(case / "50" / "U",
                          boundary=("    inlet", "    {", "        type fixedValue;",
                                    "        value uniform (0 0 0);", "    }"))
        bare._casePath = str(case)
        _, cellStart, nCells = OFLSMToolkit._parseOFFieldHeader(str(case / "50" / "U"))
        cellData = pandas.DataFrame(dict(x=[0.0, 1.0, 2.0], y=[0.0, 0.0, 0.0],
                                         height=[1.0, 1.5, 2.0]))
        with pytest.raises(pandas.errors.ParserError, match="Expected 3 fields"):
            bare._extractVelocityField(cellData, cellStart, nCells, 50)


# ---------------------------------------------------------------------------
# makeUstar's cached path (and _loadOrComputeCellData's cached branch)
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestMakeUstarFromCache:
    @pytest.fixture()
    def cachedCase(self, lsm, tmp_path):
        """A case whose cell data and u* are both already in the cache, so
        neither the mesh nor the topography has to be read."""
        case = tmp_path / "case"
        (case / "50").mkdir(parents=True)
        _writeVectorField(case / "50" / "U")
        (case / "50" / "Hmix").write_text("\n".join(["boundaryField", "{", "    ground {}", "}"]) + "\n")
        lsm.casePath = str(case)

        cellDataFile = tmp_path / "cellData.parquet"
        pandas.DataFrame(dict(x=[0.0, 1.0, 2.0], y=[0.0, 0.0, 0.0],
                              height=[1.0, 1.5, 2.0])).to_parquet(str(cellDataFile))
        lsm.topography.addCacheDocument(resource=str(cellDataFile), dataFormat="parquet",
                                        type="cellData",
                                        desc=dict(resolution=10, casePath=str(case)))

        ustarFile = tmp_path / "ustar.parquet"
        pandas.DataFrame(dict(ustar=[0.25, 0.35, 0.45])).to_parquet(str(ustarFile))
        lsm.addCacheDocument(resource=str(ustarFile), dataFormat="parquet", type="ustar",
                             desc=dict(casePath=str(case), time=50))
        return case

    def test_the_cached_ustar_is_written_as_a_scalar_field(self, lsm, cachedCase):
        lsm.makeUstar(times=[50])
        written = (cachedCase / "50" / "ustar").read_text()
        assert "volScalarField" in written
        assert written.count("\n0.") == 3

    def test_the_field_file_is_named_after_the_filename_argument(self, lsm, cachedCase):
        lsm.makeUstar(times=[50], fileName="frictionVelocity")
        assert (cachedCase / "50" / "frictionVelocity").exists()

    def test_the_cached_cell_data_is_used_instead_of_the_mesh(self, lsm, cachedCase, monkeypatch):
        # getCellDataAndGroundData would need a real polyMesh; if the cache
        # branch is taken it is never called.
        def fail(*args, **kwargs):
            raise AssertionError("the mesh must not be read when the cache is warm")

        monkeypatch.setattr(lsmModule, "getCellDataAndGroundData", fail)
        lsm.makeUstar(times=[50])
        assert (cachedCase / "50" / "ustar").exists()


# ---------------------------------------------------------------------------
# createRootCaseMeshLink
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestCreateRootCaseMeshLink:
    @pytest.mark.xfail(
        strict=True,
        reason="B207: the body reads `lastTS`, a name that is never "
               "defined anywhere in the module (and hard-codes '3600' as "
               "the destination time step), so the method raises NameError "
               "for any root case that has at least one processor "
               "directory. See the consolidated findings issue.",
    )
    def test_it_should_link_the_mesh_of_every_processor(self, bare, tmp_path):
        (tmp_path / "processor0" / "constant" / "polyMesh").mkdir(parents=True)
        bare.createRootCaseMeshLink(str(tmp_path))

    def test_it_currently_raises_nameerror(self, bare, tmp_path):
        """Characterisation of B207."""
        (tmp_path / "processor0").mkdir()
        with pytest.raises(NameError, match="lastTS"):
            bare.createRootCaseMeshLink(str(tmp_path))

    def test_a_root_case_without_processors_does_nothing(self, bare, tmp_path):
        assert bare.createRootCaseMeshLink(str(tmp_path)) is None


# ---------------------------------------------------------------------------
# The presentation methods
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestToParaviewCSV:
    def test_one_file_per_time_step_holds_the_global_positions(self, bare, tmp_path):
        outputDirectory = tmp_path / "csv"
        outputDirectory.mkdir()
        bare.to_paraview_CSV(_particleData(), str(outputDirectory), "particles")
        assert sorted(os.listdir(outputDirectory)) == ["particles_1.csv", "particles_2.csv"]
        content = (outputDirectory / "particles_2.csv").read_text().split()
        assert content == ["globalX,globalY,globalZ", "0.5,0.5,0.5"]

    def test_the_time_factor_scales_the_file_name(self, bare, tmp_path):
        outputDirectory = tmp_path / "csv"
        outputDirectory.mkdir()
        bare.to_paraview_CSV(_particleData(), str(outputDirectory), "particles", timeFactor=10)
        assert sorted(os.listdir(outputDirectory)) == ["particles_10.csv", "particles_20.csv"]


@pytest.mark.unit
class TestToUnstructuredVTK:
    @pytest.fixture()
    def vtkSeam(self, monkeypatch):
        recorder = _Recorder()
        monkeypatch.setattr(lsmModule, "pointsToVTK", recorder)
        return recorder

    def test_one_dataset_is_written_per_time_step(self, bare, tmp_path, vtkSeam):
        bare.toUnstructuredVTK(_particleData(), str(tmp_path / "out"), "particles")
        assert [call[0][0] for call in vtkSeam.calls] == [str(tmp_path / "out" / "particles_1"),
                                                          str(tmp_path / "out" / "particles_2")]

    def test_the_output_directory_is_created(self, bare, tmp_path, vtkSeam):
        bare.toUnstructuredVTK(_particleData(), str(tmp_path / "out"), "particles")
        assert (tmp_path / "out").is_dir()

    def test_the_running_index_names_the_files_when_asked(self, bare, tmp_path, vtkSeam):
        bare.toUnstructuredVTK(_particleData(), str(tmp_path), "particles", timeNameOutput=False)
        assert [call[0][0] for call in vtkSeam.calls] == [str(tmp_path / "particles_0"),
                                                          str(tmp_path / "particles_1")]

    def test_the_mass_and_the_positions_of_the_time_step_are_handed_over(self, bare, tmp_path, vtkSeam):
        bare.toUnstructuredVTK(_particleData(), str(tmp_path), "particles")
        path, x, y, z, fields = vtkSeam.calls[0][0]
        assert list(x) == [0.5, 1.5]
        assert list(fields["mass"]) == [1.0, 2.0]


@pytest.mark.unit
class TestToStructuredVTK:
    def test_it_currently_cannot_run(self, bare, tmp_path, monkeypatch):
        """Characterisation of B208, from its caller.

        toStructuredVTK is blocked by the same DataArray attribute
        assignment as calcConcentrationTimeStepFullMesh.
        """
        monkeypatch.setattr(lsmModule, "structuredToVTK", _Recorder())
        bare._analysis = Analysis(bare)
        with pytest.raises(AttributeError, match="filterType"):
            bare.toStructuredVTK(_particleData(), str(tmp_path), "concentration",
                                 extents=_EXTENTS, dxdydz=1)


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestAnalysisConcentration:
    @pytest.fixture()
    def analysis(self, unit_project):
        # Analysis only needs the datalayer's cache/filesDirectory API, so a
        # project is enough -- and keeps the broken constructor out of it.
        return Analysis(unit_project)

    def test_the_datalayer_is_the_one_it_was_built_with(self, analysis, unit_project):
        assert analysis.datalayer is unit_project

    def test_the_mass_of_a_cell_is_summed_and_divided_by_the_volume(self, analysis):
        data = pandas.DataFrame(dict(time=[1.0, 1.0], globalX=[0.5, 0.6], globalY=[0.5, 0.5],
                                     globalZ=[0.5, 0.5], mass=[1.0, 2.0]))
        concentration = analysis.calcConcentrationPointWise(data, dxdydz=2.0)
        assert concentration["C"].tolist() == [3.0 / 8.0]

    def test_particles_in_different_cells_stay_apart(self, analysis):
        data = pandas.DataFrame(dict(time=[1.0, 1.0], globalX=[0.5, 5.5], globalY=[0.5, 0.5],
                                     globalZ=[0.5, 0.5], mass=[1.0, 2.0]))
        assert sorted(analysis.calcConcentrationPointWise(data, dxdydz=1.0)["C"].tolist()) == [1.0, 2.0]

    @pytest.mark.xfail(
        strict=True,
        reason="B208: calcConcentrationTimeStepFullMesh does "
               "`fulldata.filterType = 'C'` on an xarray DataArray, which "
               "refuses attribute assignment (it must go through .attrs, as "
               "the units are set two lines later). Every call raises "
               "AttributeError, which also takes down "
               "calcConcentrationFieldFullMesh and both toStructuredVTK "
               "implementations. See the consolidated findings issue.",
    )
    def test_a_time_step_should_be_embedded_in_the_full_mesh(self, analysis):
        data = pandas.DataFrame(dict(time=[1.0], globalX=[0.5], globalY=[0.5], globalZ=[0.5],
                                     mass=[1.0]))
        field = analysis.calcConcentrationTimeStepFullMesh(data, extents=_EXTENTS, dxdydz=1.0)
        assert field.attrs["units"] == "1*kg/m**3"

    def test_a_time_step_currently_raises(self, analysis):
        """Characterisation of B208."""
        data = pandas.DataFrame(dict(time=[1.0], globalX=[0.5], globalY=[0.5], globalZ=[0.5],
                                     mass=[1.0]))
        with pytest.raises(AttributeError, match="filterType"):
            analysis.calcConcentrationTimeStepFullMesh(data, extents=_EXTENTS, dxdydz=1.0)

    def test_the_full_mesh_field_registers_its_cache_document_before_failing(self, analysis,
                                                                            unit_project):
        """Characterisation of B208, from calcConcentrationFieldFullMesh.

        The cache document (and its Concentrations*.nc resource) is written
        first, so a failed run leaves a document behind pointing at files
        that were never produced. Behind this failure sits a second latent
        typo: the transpose asks for a dimension called "time " (trailing
        space).
        """
        import dask.dataframe

        data = pandas.DataFrame(dict(time=[1.0], globalX=[0.5], globalY=[0.5], globalZ=[0.5],
                                     mass=[1.0]))

        class Document:
            id = "deadbeef"

            def getData(self):
                return dask.dataframe.from_pandas(data, npartitions=1)

        with pytest.raises(AttributeError, match="filterType"):
            analysis.calcConcentrationFieldFullMesh(Document(), extents=_EXTENTS, dxdydz=1.0)
        cached = unit_project.getCacheDocuments(type=Analysis.DOCTYPE_CONCENTRATION)
        assert len(cached) == 1
        assert cached[0].desc["dataID"] == "deadbeef"
        assert cached[0].resource.endswith("Concentrations*.nc")


@pytest.mark.unit
class TestAnalysisCalcDocumentConcentrationPointWise:
    @pytest.fixture()
    def analysis(self, unit_project):
        return Analysis(unit_project)

    @pytest.fixture()
    def document(self):
        import dask.dataframe

        data = pandas.DataFrame(dict(time=[1.0, 1.0], globalX=[0.5, 0.6], globalY=[0.5, 0.5],
                                     globalZ=[0.5, 0.5], mass=[1.0, 2.0]))

        class Document:
            id = "deadbeef"

            def getData(self):
                return dask.dataframe.from_pandas(data, npartitions=1)

        return Document()

    def test_the_cache_document_records_the_simulation_and_the_cell_size(self, analysis, document,
                                                                        unit_project):
        analysis.calcDocumentConcentrationPointWise(document, dxdydz=1.0, saveAsDask=True)
        cached = unit_project.getCacheDocuments(type=Analysis.DOCTYPE_CONCENTRATION_POINTWISE)
        assert len(cached) == 1
        assert cached[0].desc["simID"] == "deadbeef"
        assert cached[0].desc["dxdydz"] == 1.0

    def test_the_extra_metadata_reaches_the_cache_document(self, analysis, document, unit_project):
        analysis.calcDocumentConcentrationPointWise(document, dxdydz=1.0, saveAsDask=True,
                                                    experiment="myExperiment")
        cached = unit_project.getCacheDocuments(type=Analysis.DOCTYPE_CONCENTRATION_POINTWISE)
        assert cached[0].desc["experiment"] == "myExperiment"

    def test_an_explicit_simulation_id_overrides_the_documents(self, analysis, document, unit_project):
        analysis.calcDocumentConcentrationPointWise(document, dxdydz=1.0, saveAsDask=True,
                                                    simulationID="fromAnotherDB")
        cached = unit_project.getCacheDocuments(type=Analysis.DOCTYPE_CONCENTRATION_POINTWISE)
        assert cached[0].desc["simID"] == "fromAnotherDB"

    def test_the_concentration_is_written_next_to_the_document(self, analysis, document, unit_project):
        analysis.calcDocumentConcentrationPointWise(document, dxdydz=1.0, saveAsDask=True)
        cached = unit_project.getCacheDocuments(type=Analysis.DOCTYPE_CONCENTRATION_POINTWISE)
        assert os.path.exists(cached[0].resource)
        assert cached[0].resource.endswith("concentration.parquet")

    def test_a_second_call_reuses_the_cache_document(self, analysis, document, unit_project):
        analysis.calcDocumentConcentrationPointWise(document, dxdydz=1.0, saveAsDask=True)
        analysis.calcDocumentConcentrationPointWise(document, dxdydz=1.0, saveAsDask=True,
                                                    overwrite=True)
        assert len(unit_project.getCacheDocuments(type=Analysis.DOCTYPE_CONCENTRATION_POINTWISE)) == 1

    @pytest.mark.xfail(
        strict=True,
        reason="B209: the method documents ':return: Document' and its "
               "last statement is `ret = doc`, but there is no return "
               "statement at all -- so it always answers None and the caller "
               "cannot reach the document (or the parquet) it just built. "
               "See the consolidated findings issue.",
    )
    def test_it_should_return_the_cache_document(self, analysis, document):
        assert analysis.calcDocumentConcentrationPointWise(document, dxdydz=1.0,
                                                           saveAsDask=True) is not None

    def test_it_currently_returns_none(self, analysis, document):
        """Characterisation of B209."""
        assert analysis.calcDocumentConcentrationPointWise(document, dxdydz=1.0,
                                                           saveAsDask=True) is None


@pytest.mark.unit
class TestAnalysisGetConcentrationField:
    @pytest.fixture()
    def analysis(self, unit_project):
        return Analysis(unit_project)

    @pytest.fixture()
    def document(self):
        class Document:
            id = "deadbeef"

        return Document()

    def test_an_unknown_simulation_gives_none(self, analysis, document):
        assert analysis.getConcentrationField(document) is None

    def test_the_first_matching_document_is_returned(self, analysis, document, unit_project):
        unit_project.addCacheDocument(resource="/x", dataFormat="string",
                                      type=Analysis.DOCTYPE_CONCENTRATION,
                                      desc=dict(dataID="deadbeef", dxdydz=1.0))
        assert analysis.getConcentrationField(document).desc["dataID"] == "deadbeef"

    def test_the_extra_metadata_narrows_the_search(self, analysis, document, unit_project):
        unit_project.addCacheDocument(resource="/x", dataFormat="string",
                                      type=Analysis.DOCTYPE_CONCENTRATION,
                                      desc=dict(dataID="deadbeef", dxdydz=1.0))
        assert analysis.getConcentrationField(document, dxdydz=1.0) is not None
        assert analysis.getConcentrationField(document, dxdydz=99.0) is None

    @pytest.mark.xfail(
        strict=True,
        reason="B210: returnFirst=False is documented as 'returns a list "
               "(might be empty)', but the else branch returns "
               "docList[0].getData(usePandas=True) -- the data of the first "
               "document rather than a list, and an IndexError when nothing "
               "matched. See the consolidated findings issue.",
    )
    def test_return_first_false_should_give_a_list(self, analysis, document):
        assert analysis.getConcentrationField(document, returnFirst=False) == []

    def test_return_first_false_currently_raises_on_an_empty_result(self, analysis, document):
        """Characterisation of B210."""
        with pytest.raises(IndexError):
            analysis.getConcentrationField(document, returnFirst=False)


@pytest.mark.unit
class TestAnalysisGetMassFromLog:
    @pytest.fixture()
    def analysis(self, bare):
        return Analysis(bare)

    @staticmethod
    def _log(tmp_path, lines):
        path = tmp_path / "solver.log"
        path.write_text("\n".join(lines) + "\n")
        return str(path)

    def test_a_release_record_is_parsed(self, analysis, tmp_path):
        logFile = self._log(tmp_path, [
            "Exec   : StochasticLagrangianSolver",
            "Build  : whatever",
            "",
            "Time = 100",
            "Cloud: kinematicCloud",
            "    Current number of parcels 5",
            "source1:",
            "    parcels added 10",
            "    mass introduced 2.5",
            "ExecutionTime = 1 s",
            "",
            "End",
            "",
        ])
        table = analysis.getMassFromLog(logFile)
        assert table.to_dict("records") == [dict(time=100.0, cloudName="kinematicCloud",
                                                 name="source1", action="release",
                                                 parcels=10.0, mass=2.5)]

    def test_a_parcel_fate_record_becomes_an_escape_and_a_stick(self, analysis, tmp_path):
        logFile = self._log(tmp_path, [
            "Exec   : StochasticLagrangianSolver",
            "Build  : whatever",
            "",
            "Time = 200",
            "Cloud: kinematicCloud",
            "    Current number of parcels 5",
            "    Parcel fate: patch ground",
            "        escape  3, 1.5",
            "        stick   4, 2.5",
            "ExecutionTime = 1 s",
            "",
            "End",
            "",
        ])
        table = analysis.getMassFromLog(logFile)
        assert table["action"].tolist() == ["escape", "stick"]
        assert table["parcels"].tolist() == [3.0, 4.0]
        assert table["mass"].tolist() == [1.5, 2.5]
        assert table["name"].tolist() == ["ground", "ground"]

    def test_the_solver_name_selects_where_the_parsing_starts(self, analysis, tmp_path):
        logFile = self._log(tmp_path, [
            "Time = 1",
            "this line is before the Exec banner and must be ignored",
            "Exec   : myOwnSolver",
            "Build  : whatever",
            "",
            "End",
            "",
        ])
        assert analysis.getMassFromLog(logFile, solver="myOwnSolver").empty

    def test_the_cloud_name_comes_from_the_datalayer(self, analysis, bare, tmp_path):
        bare.cloudName = "myCloud"
        logFile = self._log(tmp_path, [
            "Exec   : StochasticLagrangianSolver",
            "Build  : whatever",
            "",
            "Time = 100",
            "Cloud: myCloud",
            "    Current number of parcels 5",
            "source1:",
            "    parcels added 10",
            "    mass introduced 2.5",
            "ExecutionTime = 1 s",
            "",
            "End",
            "",
        ])
        assert analysis.getMassFromLog(logFile)["cloudName"].tolist() == ["myCloud"]
