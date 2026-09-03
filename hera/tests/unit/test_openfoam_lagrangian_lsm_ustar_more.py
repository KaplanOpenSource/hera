"""OFLSMToolkit: the two members ``test_openfoam_lagrangian_lsm_toolkit_more.py``
listed as deliberately not covered.

That file's docstring rules out ``makeDistanceFromGround`` ("needs
``getCellDataAndGroundData`` to read a real ``polyMesh`` plus a ``ground``
patch out of an OpenFOAM case, and then ``topography.analysis.addHeight`` on
top of it") and ``_interpolateNearGroundValues`` ("a
``coordinateHandler.regularizeTimeSteps`` sandwich").  Both turn out to be
reachable:

* ``makeDistanceFromGround`` has a *cached* path -- the same shape
  ``makeUstar`` uses -- in which a ``cellData`` cache document on the
  topography toolkit replaces the mesh read entirely.  Its uncached path is
  reachable too, with ``getCellDataAndGroundData`` (a module-level import, so
  monkeypatchable on the module) and ``addHeight`` (patched on the topography
  analysis *class*) used as recording seams.  What is asserted is the file the
  method writes -- byte for byte the header of the input ``cellCenters``, one
  ``(x y height)`` triplet per cell, and the boundary section of ``Hmix`` --
  and which arguments crossed the two seams.
* ``_interpolateNearGroundValues`` does not need the helper tested, because it
  never reaches it: the name it calls does not exist.

The toolkit is built through ``unit_toolkit_factory`` with B103's missing
``ToolkitHome.GIS_TOPOGRAPHY`` constant patched onto the *class* for the
duration of the test (the technique the sibling file introduced; B103 itself
stays pinned there, unpatched).

Bug pinned here (with a strict xfail for the intended behaviour and a passing
characterisation of what happens today):

* B292 ``_interpolateNearGroundValues`` calls
  ``coordinateHandler.regularizeTimeSteps(...)``, but the module-level import
  is ``from hera.simulations.utils import coordinateHandler``, which binds the
  *module* -- and ``regularizeTimeSteps`` is a method of the
  ``coordinateHandler`` *class* inside it, not a module-level function.  Every
  call raises ``AttributeError: module ... has no attribute
  'regularizeTimeSteps'``.  That kills ``makeUstar``'s uncached path with it:
  the only way ``makeUstar`` can produce a field today is from a warm u*
  cache, which is exactly the one path the sibling file covers.
"""
import os

import pandas
import pytest

from hera.simulations.openFoam.lagrangian.LSM import toolkit as lsmModule
from hera.simulations.openFoam.lagrangian.LSM.toolkit import OFLSMToolkit
from hera.simulations.utils import coordinateHandler as coordinateHandlerModule

_CELL_CENTRES = [(0.0, 0.0, 1.0), (10.0, 0.0, 2.0), (20.0, 0.0, 3.0)]
#: Deliberately disjoint from the z of the centres above, so a test can tell
#: the height above ground apart from the altitude of the cell.
_HEIGHTS = [0.25, 0.5, 0.75]


def _writeCellCentresFile(path, centres=_CELL_CENTRES):
    """An OpenFOAM ``C`` field: header, count, values, then a boundary section."""
    lines = ["FoamFile", "{", "    class volVectorField;", "}",
             "dimensions      [0 1 0 0 0 0 0];", "",
             "internalField   nonuniform List<vector>", f"{len(centres)}", "("]
    lines += [f"({x} {y} {z})" for x, y, z in centres]
    lines += [")", ";", "", "boundaryField", "{", "    ground", "    {",
              "        type            calculated;", "    }", "}"]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def _writeHmixFile(path):
    """The template the boundary section of the output is copied from."""
    lines = ["FoamFile", "{", "    class volScalarField;", "}", "",
             "boundaryField", "{", "    ground", "    {",
             "        type            zeroGradient;", "    }", "}"]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def _writeVelocityFile(path, magnitudes=(1.0, 2.0, 3.0)):
    """A ``U`` field with no multi-token line below the values (see B211)."""
    lines = ["FoamFile", "{", "    class volVectorField;", "}", "",
             "internalField   nonuniform List<vector>", f"{len(magnitudes)}", "("]
    lines += [f"({magnitude} 0 0)" for magnitude in magnitudes]
    lines += [")", ";"]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture()
def lsm(unit_toolkit_factory, monkeypatch):
    """A real OF_LSM toolkit, with B103's missing constant supplied.

    The constant is patched onto the ToolkitHome *class*, never onto the
    toolkitHome singleton, which monkeypatch would restore as a permanent
    instance attribute.
    """
    from hera import ToolkitHome, toolkitHome

    monkeypatch.setattr(ToolkitHome, "GIS_TOPOGRAPHY", toolkitHome.GIS_RASTER_TOPOGRAPHY,
                        raising=False)
    return unit_toolkit_factory("OF_LSM")


@pytest.fixture()
def bare():
    """The toolkit without any constructor -- the sibling file's technique."""
    return OFLSMToolkit.__new__(OFLSMToolkit)


@pytest.fixture()
def case(lsm, tmp_path):
    """A case with cellCenters, Hmix and U, and two time directories."""
    directory = tmp_path / "case"
    _writeCellCentresFile(directory / "0" / "cellCenters")
    _writeHmixFile(directory / "0" / "Hmix")
    for time in ("50", "100"):
        (directory / time).mkdir(parents=True, exist_ok=True)
        _writeHmixFile(directory / time / "Hmix")
        _writeVelocityFile(directory / time / "U")
    lsm.casePath = str(directory)
    return directory


@pytest.fixture()
def warmCellDataCache(lsm, case, tmp_path):
    """A cellData cache document, so the mesh never has to be read."""
    cacheFile = tmp_path / "cellData.parquet"
    pandas.DataFrame(dict(x=[centre[0] for centre in _CELL_CENTRES],
                          y=[centre[1] for centre in _CELL_CENTRES],
                          height=list(_HEIGHTS))).to_parquet(str(cacheFile))
    lsm.topography.addCacheDocument(resource=str(cacheFile), dataFormat="parquet",
                                    type="cellData",
                                    desc=dict(resolution=10, casePath=str(case)))
    return cacheFile


@pytest.fixture()
def refuseMeshRead(monkeypatch):
    """Fail loudly if the mesh reader is reached."""
    def refuse(*args, **kwargs):
        raise AssertionError("the mesh must not be read when the cache is warm")

    monkeypatch.setattr(lsmModule, "getCellDataAndGroundData", refuse)


# ---------------------------------------------------------------------------
# makeDistanceFromGround -- the cached path
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestMakeDistanceFromGroundFromCache:
    def test_the_field_file_is_written_into_every_requested_time(self, lsm, case,
                                                                 warmCellDataCache,
                                                                 refuseMeshRead):
        lsm.makeDistanceFromGround(times=["50", "100"])
        assert (case / "50" / "cellHeights").is_file()
        assert (case / "100" / "cellHeights").is_file()

    def test_the_file_is_named_after_the_filename_argument(self, lsm, case,
                                                           warmCellDataCache,
                                                           refuseMeshRead):
        lsm.makeDistanceFromGround(times=["50"], fileName="distanceFromWalls")
        assert (case / "50" / "distanceFromWalls").is_file()

    def test_the_header_of_the_cellcentres_file_is_copied_verbatim(self, lsm, case,
                                                                   warmCellDataCache,
                                                                   refuseMeshRead):
        lsm.makeDistanceFromGround(times=["50"])
        written = (case / "50" / "cellHeights").read_text().split("\n")
        source = (case / "0" / "cellCenters").read_text().split("\n")
        # "internalField" is source line 6, so cellStart is 9: lines 0-8 are copied
        assert written[:9] == source[:9]

    def test_one_triplet_is_written_per_cell_of_the_mesh(self, lsm, case, warmCellDataCache,
                                                         refuseMeshRead):
        lsm.makeDistanceFromGround(times=["50"])
        written = (case / "50" / "cellHeights").read_text()
        triplets = [line for line in written.split("\n") if line.startswith("(") and " " in line]
        assert len(triplets) == len(_CELL_CENTRES)

    def test_each_triplet_is_the_x_the_y_and_the_height_above_ground(self, lsm, case,
                                                                     warmCellDataCache,
                                                                     refuseMeshRead):
        """The z of the mesh is replaced by the height from the cell data --
        that is the whole point of the method."""
        lsm.makeDistanceFromGround(times=["50"])
        written = (case / "50" / "cellHeights").read_text()
        for (x, y, _), height in zip(_CELL_CENTRES, _HEIGHTS):
            assert f"({x} {y} {height})" in written

    def test_the_altitude_of_the_mesh_does_not_survive_into_the_output(self, lsm, case,
                                                                       warmCellDataCache,
                                                                       refuseMeshRead):
        """The third component is the height above ground, not the z of the
        cell centre -- the two are deliberately disjoint in this fixture."""
        lsm.makeDistanceFromGround(times=["50"])
        valueList = (case / "50" / "cellHeights").read_text().split(")\n;")[0]
        for x, y, z in _CELL_CENTRES:
            assert f"({x} {y} {z})" not in valueList

    def test_the_boundary_section_comes_from_the_hmix_template(self, lsm, case,
                                                               warmCellDataCache,
                                                               refuseMeshRead):
        lsm.makeDistanceFromGround(times=["50"])
        written = (case / "50" / "cellHeights").read_text()
        hmix = (case / "0" / "Hmix").read_text()
        assert written.endswith(hmix[hmix.index("boundaryField"):])

    def test_the_value_list_is_closed_the_way_openfoam_expects(self, lsm, case,
                                                               warmCellDataCache,
                                                               refuseMeshRead):
        lsm.makeDistanceFromGround(times=["50"])
        written = (case / "50" / "cellHeights").read_text()
        assert ")\n;\n\nboundaryField" in written

    def test_it_returns_nothing_and_only_writes(self, lsm, case, warmCellDataCache,
                                               refuseMeshRead):
        assert lsm.makeDistanceFromGround(times=["50"]) is None

    def test_a_cache_entry_for_another_resolution_is_not_used(self, lsm, case,
                                                              warmCellDataCache,
                                                              monkeypatch):
        """The cache is keyed by resolution, so asking for a different one has
        to fall through to the mesh."""
        reached = []
        monkeypatch.setattr(lsmModule, "getCellDataAndGroundData",
                            lambda **kwargs: reached.append(kwargs) or
                            (pandas.DataFrame(), pandas.DataFrame()))
        monkeypatch.setattr(type(lsm.topography.analysis), "addHeight",
                            lambda self, **kwargs: pandas.DataFrame(
                                dict(x=[0.0, 10.0, 20.0], y=[0.0, 0.0, 0.0],
                                     height=[9.0, 8.0, 7.0])),
                            raising=False)
        lsm.makeDistanceFromGround(times=["50"], resolution=25)
        assert len(reached) == 1


# ---------------------------------------------------------------------------
# makeDistanceFromGround -- the uncached path
# ---------------------------------------------------------------------------

@pytest.fixture()
def meshSeams(lsm, monkeypatch):
    """Record the mesh read and the topography height computation."""
    recorded = dict(mesh=[], addHeight=[])

    def readMesh(**kwargs):
        recorded["mesh"].append(kwargs)
        return pandas.DataFrame(dict(x=[centre[0] for centre in _CELL_CENTRES],
                                     y=[centre[1] for centre in _CELL_CENTRES],
                                     Cz=[centre[2] for centre in _CELL_CENTRES])), \
            pandas.DataFrame(dict(x=[0.0], y=[0.0], z=[0.0]))

    def addHeight(self, **kwargs):
        recorded["addHeight"].append(kwargs)
        return pandas.DataFrame(dict(x=[centre[0] for centre in _CELL_CENTRES],
                                     y=[centre[1] for centre in _CELL_CENTRES],
                                     height=list(_HEIGHTS)))

    monkeypatch.setattr(lsmModule, "getCellDataAndGroundData", readMesh)
    monkeypatch.setattr(type(lsm.topography.analysis), "addHeight", addHeight, raising=False)
    return recorded


@pytest.mark.unit
class TestMakeDistanceFromGroundFromTheMesh:
    def test_the_mesh_is_read_from_the_case_when_the_cache_is_cold(self, lsm, case,
                                                                   meshSeams):
        lsm.makeDistanceFromGround(times=["50"])
        assert meshSeams["mesh"] == [dict(casePath=str(case), ground="ground")]

    def test_the_ground_patch_name_is_an_argument(self, lsm, case, meshSeams):
        lsm.makeDistanceFromGround(times=["50"], ground="terrain")
        assert meshSeams["mesh"][0]["ground"] == "terrain"

    def test_the_resolution_reaches_the_height_computation(self, lsm, case, meshSeams):
        lsm.makeDistanceFromGround(times=["50"], resolution=25)
        assert meshSeams["addHeight"][0]["resolution"] == 25

    def test_the_save_mode_and_the_fill_value_reach_it_too(self, lsm, case, meshSeams):
        lsm.makeDistanceFromGround(times=["50"], saveMode="Nothing", fillna=-999)
        assert meshSeams["addHeight"][0]["saveMode"] == "Nothing"
        assert meshSeams["addHeight"][0]["fillna"] == -999

    def test_the_parquet_it_is_told_to_write_sits_beside_the_case(self, lsm, case,
                                                                  meshSeams):
        lsm.makeDistanceFromGround(times=["50"], fileName="cellHeights")
        assert meshSeams["addHeight"][0]["file"] == \
            os.path.join(str(case), "cellHeights.parquet")

    def test_the_heights_it_returns_are_the_ones_written_into_the_field(self, lsm, case,
                                                                        meshSeams):
        lsm.makeDistanceFromGround(times=["50"])
        written = (case / "50" / "cellHeights").read_text()
        for (x, y, _), height in zip(_CELL_CENTRES, _HEIGHTS):
            assert f"({x} {y} {height})" in written

    def test_a_cellcentres_file_without_an_internalfield_marker_is_not_reported(
            self, lsm, case, meshSeams):
        """Characterisation, not a pinned bug: the loop that finds the marker
        leaves cellStart/nCells unbound if it is absent, so a malformed case
        surfaces as UnboundLocalError rather than as a message naming the
        file."""
        (case / "0" / "cellCenters").write_text("FoamFile\n{\n}\n")
        with pytest.raises(UnboundLocalError):
            lsm.makeDistanceFromGround(times=["50"])


# ---------------------------------------------------------------------------
# _interpolateNearGroundValues -- B292
# ---------------------------------------------------------------------------

def _nearGroundFrame():
    """Cells straddling the 1-2 m band, on a 3x2 horizontal grid."""
    return pandas.DataFrame(dict(x=[0.0, 10.0, 20.0, 0.0, 10.0, 20.0],
                                 y=[0.0, 0.0, 0.0, 10.0, 10.0, 10.0],
                                 height=[1.5, 1.5, 1.5, 5.0, 5.0, 5.0],
                                 U=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]))


@pytest.mark.unit
class TestInterpolateNearGroundValues:
    @pytest.mark.xfail(
        strict=True,
        reason="B292: the body calls "
               "coordinateHandler.regularizeTimeSteps(...), but the module "
               "imports `from hera.simulations.utils import coordinateHandler`, "
               "which binds the MODULE -- regularizeTimeSteps is a method of "
               "the coordinateHandler CLASS inside it.  Every call raises "
               "AttributeError, which also makes makeUstar's uncached path "
               "unreachable.  See the consolidated findings issue.",
    )
    def test_it_should_add_a_near_ground_velocity_column(self, bare):
        result = bare._interpolateNearGroundValues(_nearGroundFrame(), [1.0, 2.0], 10,
                                                   None, "D")
        assert "UnearGround" in result.columns
        assert "D" in result.columns

    def test_it_currently_raises_attributeerror_on_the_module(self, bare):
        """Characterisation of B292."""
        with pytest.raises(AttributeError, match="regularizeTimeSteps"):
            bare._interpolateNearGroundValues(_nearGroundFrame(), [1.0, 2.0], 10, None, "D")

    def test_the_name_it_wants_lives_on_the_class_not_on_the_module(self):
        """Characterisation of B292: what the call should have been."""
        assert not hasattr(coordinateHandlerModule, "regularizeTimeSteps")
        assert hasattr(coordinateHandlerModule.coordinateHandler, "regularizeTimeSteps")

    def test_a_displacement_field_does_not_change_the_outcome(self, bare):
        """Characterisation of B292: the failure is on the first regularize
        call, before the dField branch is reached."""
        dField = pandas.DataFrame(dict(x=[0.0, 10.0], y=[0.0, 10.0], D=[0.5, 0.5]))
        with pytest.raises(AttributeError, match="regularizeTimeSteps"):
            bare._interpolateNearGroundValues(_nearGroundFrame(), [1.0, 2.0], 10, dField, "D")

    def test_makeustar_cannot_compute_a_cold_cache_because_of_it(self, lsm, case,
                                                                 warmCellDataCache,
                                                                 refuseMeshRead):
        """Characterisation of B292: with the cell data cached but no u*
        cached, makeUstar walks into _interpolateNearGroundValues and dies --
        so a u* field can only ever come out of a warm u* cache."""
        with pytest.raises(AttributeError, match="regularizeTimeSteps"):
            lsm.makeUstar(times=["50"])
        assert not (case / "50" / "ustar").exists()
