"""datasetToOF: the meteorological-dataset-to-OpenFOAM converter (issue #1060).

The module replaces two scratch scripts that never ran -- both
``toberewritten/netcdf2of.py`` and ``toberewritten/xarrayDataset2OF.py`` raise
``SyntaxError`` when compiled -- so there is no legacy behaviour to pin down.
These tests fix the behaviour of the replacement instead:

* the geometry of a case is read from the ``C`` field through foamlib, for a
  reconstructed and for a decomposed case, and a patch that carries a single
  uniform vector (so its per-face centres are unknown) is skipped rather than
  invented;
* interpolation of a field that is linear in x, y and z is exact -- which is
  what pins down the bilinear horizontal / linear vertical scheme, including
  the terrain-following ``heightVariable`` path;
* a point outside the vertical extent takes the value of the nearest level,
  i.e. the profile is not extrapolated;
* the coordinates of the case are projected into those of the dataset (and
  back, for ``setFieldsDict``), which the old scripts did with a proj4 string
  copied into the source;
* both outputs round-trip through foamlib, which is the check that the emitted
  ``setFieldsDict`` and the written field files are valid OpenFOAM syntax.

The last test is the guard for the third checkbox of the issue ("remove any
paths that are currently exposed in the script").
"""
import logging
import os
import re

import numpy
import pytest

pytest.importorskip("foamlib", reason="datasetToOF reads and writes the case with foamlib")

# hera.simulations.openFoam needs hermes at import time, and its
# preprocessOFObjects package needs PyFoam -- neither is pip-installable.  The
# hermetic layer satisfies both through hera/tests/unit/_stubs.py; outside it,
# the real packages have to be present (CI puts them on PYTHONPATH).  Same
# guard, and for the same reason, as hera/tests/test_of_lsm_flow.py.
pytest.importorskip("hermes", reason="hera.simulations.openFoam imports hermes at import time")
pytest.importorskip("PyFoam", reason="preprocessOFObjects imports PyFoam at import time")

import xarray
from foamlib import FoamFieldFile, FoamFile

from hera.measurements.GIS import ITM, WSG84
from hera.simulations.openFoam.preprocessOFObjects import datasetToOF

# The dataset and the case share a CRS in most tests, so that the projection
# is the identity and any error is an interpolation error.
_ITM_X = numpy.array([180000.0, 180100.0, 180200.0])
_ITM_Y = numpy.array([660000.0, 660100.0])
_LEVELS = numpy.array([10.0, 50.0, 100.0])


def _linearField(x, y, z):
    """A field that is linear in every direction, hence interpolated exactly."""
    return 2.0 * x + 3.0 * y + 5.0 * z


def _rectilinearDataset():
    """A (z, y, x) dataset in ITM metres holding two linear fields."""
    zGrid, yGrid, xGrid = numpy.meshgrid(_LEVELS, _ITM_Y, _ITM_X, indexing="ij")
    return xarray.Dataset(
        data_vars=dict(
            temperature=(("z", "y", "x"), _linearField(xGrid, yGrid, zGrid)),
            u=(("z", "y", "x"), 0.5 * xGrid),
            v=(("z", "y", "x"), -0.25 * yGrid),
        ),
        coords=dict(x=_ITM_X, y=_ITM_Y, z=_LEVELS),
    )


def _writeFoamField(path, fieldClass, internalField, boundaryField=None):
    """Write a minimal but valid OpenFOAM field file.

    Parameters
    ----------
    path : str
        Where to write it.
    fieldClass : str
        ``volScalarField`` or ``volVectorField``.
    internalField : numpy.ndarray
        ``(N,)`` for a scalar and ``(N,3)`` for a vector.
    boundaryField : dict, optional
        ``{patchName: (M,) / (M,3) array, or a single vector for a uniform patch}``.
    """
    isVector = fieldClass == "volVectorField"

    def valueText(values):
        values = numpy.asarray(values, dtype=float)
        if isVector and values.ndim == 1:
            return f"uniform ({' '.join(str(component) for component in values)})"
        if isVector:
            entries = " ".join(f"({' '.join(str(component) for component in row)})" for row in values)
            return f"nonuniform List<vector> {values.shape[0]}({entries})"
        entries = " ".join(str(value) for value in values)
        return f"nonuniform List<scalar> {values.shape[0]}({entries})"

    patches = ""
    for patchName, patchValues in (boundaryField or dict()).items():
        patches += f"""    {patchName}
    {{
        type            calculated;
        value           {valueText(patchValues)};
    }}
"""

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as foamFile:
        foamFile.write(f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       {fieldClass};
    location    "0";
    object      {os.path.basename(path)};
}}
dimensions      [0 0 0 0 0 0 0];

internalField   {valueText(internalField)};

boundaryField
{{
{patches}}}
""")


@pytest.fixture
def caseWithGeometry(tmp_path):
    """A reconstructed case with cell centres and one patch of face centres."""
    cellCentres = numpy.array([[180050.0, 660050.0, 20.0],
                               [180150.0, 660050.0, 60.0],
                               [180050.0, 660050.0, 90.0]])
    inletCentres = numpy.array([[180000.0, 660000.0, 30.0],
                                [180000.0, 660100.0, 70.0]])
    casePath = str(tmp_path / "case")
    _writeFoamField(os.path.join(casePath, "0", "C"), "volVectorField", cellCentres,
                    dict(inlet=inletCentres, top=numpy.array([0.0, 0.0, 100.0])))
    return dict(path=casePath, internalField=cellCentres, inlet=inletCentres)


# ---------------------------------------------------------------------------
#  caseGeometry
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_caseGeometry_readsCellCentresAndFaceCentres(caseWithGeometry):
    geometry = datasetToOF.caseGeometry(caseWithGeometry["path"])

    assert list(geometry) == [datasetToOF.SINGLEPROCESSOR]
    regions = geometry[datasetToOF.SINGLEPROCESSOR]
    numpy.testing.assert_allclose(regions[datasetToOF.INTERNALFIELD], caseWithGeometry["internalField"])
    numpy.testing.assert_allclose(regions["inlet"], caseWithGeometry["inlet"])


@pytest.mark.unit
def test_caseGeometry_skipsAUniformPatchInsteadOfInventingFaces(caseWithGeometry, caplog, monkeypatch):
    # hera's logging config sets 'hera.simulations' to CRITICAL with
    # propagate=False, so the warning has to be let through to be observed.
    simulationsLogger = logging.getLogger("hera.simulations")
    monkeypatch.setattr(simulationsLogger, "propagate", True)
    monkeypatch.setattr(simulationsLogger, "level", logging.WARNING)

    with caplog.at_level(logging.WARNING):
        geometry = datasetToOF.caseGeometry(caseWithGeometry["path"])

    assert "top" not in geometry[datasetToOF.SINGLEPROCESSOR], \
        "a patch with a uniform C value has no per-face centres and must be skipped"
    assert any("uniform value" in record.message for record in caplog.records)


@pytest.mark.unit
def test_caseGeometry_patchListSelectsAndValidates(caseWithGeometry):
    geometry = datasetToOF.caseGeometry(caseWithGeometry["path"], patchList=[])
    assert list(geometry[datasetToOF.SINGLEPROCESSOR]) == [datasetToOF.INTERNALFIELD]

    with pytest.raises(KeyError, match="notAPatch"):
        datasetToOF.caseGeometry(caseWithGeometry["path"], patchList=["notAPatch"])


@pytest.mark.unit
def test_caseGeometry_readsEveryProcessorOfADecomposedCase(tmp_path):
    casePath = str(tmp_path / "decomposed")
    first = numpy.array([[180050.0, 660050.0, 20.0]])
    second = numpy.array([[180150.0, 660050.0, 60.0], [180150.0, 660050.0, 80.0]])
    _writeFoamField(os.path.join(casePath, "processor0", "0", "C"), "volVectorField", first)
    _writeFoamField(os.path.join(casePath, "processor1", "0", "C"), "volVectorField", second)

    geometry = datasetToOF.caseGeometry(casePath)
    assert sorted(geometry) == ["processor0", "processor1"]
    numpy.testing.assert_allclose(geometry["processor0"][datasetToOF.INTERNALFIELD], first)
    numpy.testing.assert_allclose(geometry["processor1"][datasetToOF.INTERNALFIELD], second)

    # readParallel=False looks for the reconstructed case, which does not exist here.
    with pytest.raises(FileNotFoundError):
        datasetToOF.caseGeometry(casePath, readParallel=False)


@pytest.mark.unit
def test_caseGeometry_missingCellCentresRaisesWithTheFix(tmp_path):
    with pytest.raises(FileNotFoundError, match="writeCellCentres"):
        datasetToOF.caseGeometry(str(tmp_path / "empty"))


# ---------------------------------------------------------------------------
#  interpolateDatasetToPoints
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_interpolate_isExactForALinearField():
    points = numpy.array([[180050.0, 660050.0, 30.0],
                          [180000.0, 660000.0, 10.0],
                          [180175.0, 660025.0, 75.0]])

    values = datasetToOF.interpolateDatasetToPoints(dataset=_rectilinearDataset(),
                                                    points=points,
                                                    fieldMap=dict(T="temperature"),
                                                    xCoordinate="x",
                                                    yCoordinate="y",
                                                    verticalCoordinate="z",
                                                    datasetCRS=ITM,
                                                    caseCRS=ITM)

    numpy.testing.assert_allclose(values["T"], _linearField(points[:, 0], points[:, 1], points[:, 2]))


@pytest.mark.unit
def test_interpolate_buildsVectorsAndHonoursConstantComponents():
    points = numpy.array([[180050.0, 660050.0, 30.0], [180150.0, 660050.0, 60.0]])

    values = datasetToOF.interpolateDatasetToPoints(dataset=_rectilinearDataset(),
                                                    points=points,
                                                    fieldMap=dict(U=("u", "v", 0)),
                                                    xCoordinate="x",
                                                    yCoordinate="y",
                                                    verticalCoordinate="z",
                                                    datasetCRS=ITM,
                                                    caseCRS=ITM)

    assert values["U"].shape == (2, 3)
    numpy.testing.assert_allclose(values["U"][:, 0], 0.5 * points[:, 0])
    numpy.testing.assert_allclose(values["U"][:, 1], -0.25 * points[:, 1])
    numpy.testing.assert_allclose(values["U"][:, 2], 0.0)


@pytest.mark.unit
def test_interpolate_doesNotExtrapolateBelowOrAboveTheProfile():
    points = numpy.array([[180050.0, 660050.0, -500.0], [180050.0, 660050.0, 5000.0]])

    values = datasetToOF.interpolateDatasetToPoints(dataset=_rectilinearDataset(),
                                                    points=points,
                                                    fieldMap=dict(T="temperature"),
                                                    xCoordinate="x",
                                                    yCoordinate="y",
                                                    verticalCoordinate="z",
                                                    datasetCRS=ITM,
                                                    caseCRS=ITM)

    lowest = _linearField(180050.0, 660050.0, _LEVELS[0])
    highest = _linearField(180050.0, 660050.0, _LEVELS[-1])
    numpy.testing.assert_allclose(values["T"], [lowest, highest])


@pytest.mark.unit
def test_interpolate_handlesLevelsStoredInDescendingOrder():
    dataset = _rectilinearDataset().isel(z=slice(None, None, -1))
    points = numpy.array([[180050.0, 660050.0, 30.0]])

    values = datasetToOF.interpolateDatasetToPoints(dataset=dataset,
                                                    points=points,
                                                    fieldMap=dict(T="temperature"),
                                                    xCoordinate="x",
                                                    yCoordinate="y",
                                                    verticalCoordinate="z",
                                                    datasetCRS=ITM,
                                                    caseCRS=ITM)

    numpy.testing.assert_allclose(values["T"], _linearField(180050.0, 660050.0, 30.0))


@pytest.mark.unit
def test_interpolate_supportsATerrainFollowingVerticalCoordinate():
    # height(level, y, x) rises with x: the same model level is at a different
    # altitude in every column, which is what WRF's gpt_hgt_M looks like.
    levelIndex = numpy.arange(3.0)
    heights = numpy.stack([numpy.add.outer(numpy.zeros_like(_ITM_Y), _ITM_X - _ITM_X[0]) + 20.0 * (level + 1)
                           for level in levelIndex])
    temperature = 0.5 * heights + 7.0

    dataset = xarray.Dataset(
        data_vars=dict(height=(("level", "y", "x"), heights),
                       temperature=(("level", "y", "x"), temperature)),
        coords=dict(x=_ITM_X, y=_ITM_Y, level=levelIndex),
    )

    points = numpy.array([[_ITM_X[0], _ITM_Y[0], 30.0], [_ITM_X[0], _ITM_Y[0], 50.0]])
    values = datasetToOF.interpolateDatasetToPoints(dataset=dataset,
                                                    points=points,
                                                    fieldMap=dict(T="temperature"),
                                                    xCoordinate="x",
                                                    yCoordinate="y",
                                                    verticalCoordinate="level",
                                                    heightVariable="height",
                                                    datasetCRS=ITM,
                                                    caseCRS=ITM)

    numpy.testing.assert_allclose(values["T"], 0.5 * points[:, 2] + 7.0)


@pytest.mark.unit
def test_interpolate_nearestHandlesACurvilinearGridThatLinearRejects():
    xGrid, yGrid = numpy.meshgrid(_ITM_X, _ITM_Y, indexing="xy")
    dataset = xarray.Dataset(
        data_vars=dict(temperature=(("z", "row", "column"),
                                    numpy.stack([_linearField(xGrid, yGrid, level) for level in _LEVELS])),
                       longitudeLike=(("row", "column"), xGrid),
                       latitudeLike=(("row", "column"), yGrid)),
        coords=dict(z=_LEVELS),
    )
    points = numpy.array([[180048.0, 660010.0, 10.0]])

    nearest = datasetToOF.interpolateDatasetToPoints(dataset=dataset,
                                                     points=points,
                                                     fieldMap=dict(T="temperature"),
                                                     xCoordinate="longitudeLike",
                                                     yCoordinate="latitudeLike",
                                                     verticalCoordinate="z",
                                                     datasetCRS=ITM,
                                                     caseCRS=ITM,
                                                     horizontalMethod=datasetToOF.HORIZONTAL_NEAREST)
    # The nearest column is the grid point (180000, 660000).
    numpy.testing.assert_allclose(nearest["T"], _linearField(_ITM_X[0], _ITM_Y[0], 10.0))

    with pytest.raises(ValueError, match="curvilinear"):
        datasetToOF.interpolateDatasetToPoints(dataset=dataset,
                                               points=points,
                                               fieldMap=dict(T="temperature"),
                                               xCoordinate="longitudeLike",
                                               yCoordinate="latitudeLike",
                                               verticalCoordinate="z",
                                               datasetCRS=ITM,
                                               caseCRS=ITM,
                                               horizontalMethod=datasetToOF.HORIZONTAL_LINEAR)


@pytest.mark.unit
def test_interpolate_projectsTheCaseCoordinatesIntoThoseOfTheDataset():
    # The dataset variable *is* the longitude, so the interpolated value at an
    # ITM point must be the longitude of that point -- which only holds if the
    # point was projected instead of being used as a degree pair.
    longitudes = numpy.array([34.2, 34.6, 35.0, 35.4])
    latitudes = numpy.array([31.4, 31.8, 32.2, 32.6])
    _, longitudeGrid = numpy.meshgrid(latitudes, longitudes, indexing="ij")
    dataset = xarray.Dataset(
        data_vars=dict(longitudeValue=(("latitude", "longitude"), longitudeGrid)),
        coords=dict(longitude=longitudes, latitude=latitudes),
    )

    from pyproj import Transformer
    itmPoint = (180000.0, 660000.0)
    expectedLongitude, expectedLatitude = Transformer.from_crs(
        f"EPSG:{ITM}", f"EPSG:{WSG84}", always_xy=True).transform(*itmPoint)
    assert longitudes[0] < expectedLongitude < longitudes[-1]
    assert latitudes[0] < expectedLatitude < latitudes[-1]

    values = datasetToOF.interpolateDatasetToPoints(dataset=dataset,
                                                    points=numpy.array([[itmPoint[0], itmPoint[1], 0.0]]),
                                                    fieldMap=dict(T="longitudeValue"),
                                                    xCoordinate="longitude",
                                                    yCoordinate="latitude",
                                                    datasetCRS=WSG84,
                                                    caseCRS=ITM)

    numpy.testing.assert_allclose(values["T"], [expectedLongitude], rtol=1e-6)


@pytest.mark.unit
def test_interpolate_selectsTheNearestTimeStep():
    dataset = _rectilinearDataset().expand_dims(time=[0.0, 100.0]).copy()
    dataset["temperature"] = dataset["temperature"] + xarray.DataArray([0.0, 1000.0], dims="time")
    points = numpy.array([[180050.0, 660050.0, 30.0]])

    values = datasetToOF.interpolateDatasetToPoints(dataset=dataset,
                                                    points=points,
                                                    fieldMap=dict(T="temperature"),
                                                    xCoordinate="x",
                                                    yCoordinate="y",
                                                    verticalCoordinate="z",
                                                    time=90.0,
                                                    datasetCRS=ITM,
                                                    caseCRS=ITM)

    numpy.testing.assert_allclose(values["T"], _linearField(180050.0, 660050.0, 30.0) + 1000.0)


@pytest.mark.unit
@pytest.mark.parametrize("fieldMap, message", [
    (dict(U=("u", "v")), "components"),
    (dict(T=object()), "neither a dataset variable name nor a number"),
])
def test_interpolate_rejectsAnIllFormedFieldMap(fieldMap, message):
    with pytest.raises(ValueError, match=message):
        datasetToOF.interpolateDatasetToPoints(dataset=_rectilinearDataset(),
                                               points=numpy.zeros((1, 3)),
                                               fieldMap=fieldMap,
                                               xCoordinate="x",
                                               yCoordinate="y",
                                               verticalCoordinate="z")


@pytest.mark.unit
def test_interpolate_rejectsPointsThatAreNotTriplets():
    with pytest.raises(ValueError, match=r"\(N,3\)"):
        datasetToOF.interpolateDatasetToPoints(dataset=_rectilinearDataset(),
                                               points=numpy.zeros((4, 2)),
                                               fieldMap=dict(T="temperature"),
                                               xCoordinate="x",
                                               yCoordinate="y",
                                               verticalCoordinate="z")


@pytest.mark.unit
def test_interpolate_reportsAMissingVariable():
    with pytest.raises(KeyError, match="notThere"):
        datasetToOF.interpolateDatasetToPoints(dataset=_rectilinearDataset(),
                                               points=numpy.zeros((1, 3)),
                                               fieldMap=dict(T="notThere"),
                                               xCoordinate="x",
                                               yCoordinate="y",
                                               verticalCoordinate="z")


# ---------------------------------------------------------------------------
#  datasetToCaseFields
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_datasetToCaseFields_writesTheInternalFieldAndThePatch(caseWithGeometry):
    casePath = caseWithGeometry["path"]
    nCells = caseWithGeometry["internalField"].shape[0]
    nFaces = caseWithGeometry["inlet"].shape[0]
    _writeFoamField(os.path.join(casePath, "0", "T"), "volScalarField",
                    numpy.zeros(nCells), dict(inlet=numpy.zeros(nFaces)))
    _writeFoamField(os.path.join(casePath, "0", "U"), "volVectorField",
                    numpy.zeros((nCells, 3)), dict(inlet=numpy.zeros((nFaces, 3))))

    written = datasetToOF.datasetToCaseFields(caseDirectory=casePath,
                                              dataset=_rectilinearDataset(),
                                              fieldMap=dict(T="temperature", U=("u", "v", 0)),
                                              xCoordinate="x",
                                              yCoordinate="y",
                                              verticalCoordinate="z",
                                              datasetCRS=ITM,
                                              caseCRS=ITM)

    assert sorted(os.path.basename(path) for path in written) == ["T", "U"]

    cells = caseWithGeometry["internalField"]
    faces = caseWithGeometry["inlet"]

    temperature = FoamFieldFile(os.path.join(casePath, "0", "T"))
    numpy.testing.assert_allclose(numpy.asarray(temperature.internal_field),
                                  _linearField(cells[:, 0], cells[:, 1], cells[:, 2]))
    numpy.testing.assert_allclose(numpy.asarray(temperature.boundary_field["inlet"]["value"]),
                                  _linearField(faces[:, 0], faces[:, 1], faces[:, 2]))

    velocity = FoamFieldFile(os.path.join(casePath, "0", "U"))
    internalVelocity = numpy.asarray(velocity.internal_field)
    assert internalVelocity.shape == (cells.shape[0], 3)
    numpy.testing.assert_allclose(internalVelocity[:, 0], 0.5 * cells[:, 0])
    numpy.testing.assert_allclose(internalVelocity[:, 2], 0.0)


@pytest.mark.unit
def test_datasetToCaseFields_writesEveryProcessorOfADecomposedCase(tmp_path):
    casePath = str(tmp_path / "decomposed")
    centres = dict(processor0=numpy.array([[180050.0, 660050.0, 20.0]]),
                   processor1=numpy.array([[180150.0, 660050.0, 60.0], [180150.0, 660050.0, 80.0]]))
    for processorName, processorCentres in centres.items():
        _writeFoamField(os.path.join(casePath, processorName, "0", "C"), "volVectorField", processorCentres)
        _writeFoamField(os.path.join(casePath, processorName, "0", "T"), "volScalarField",
                        numpy.zeros(processorCentres.shape[0]))

    written = datasetToOF.datasetToCaseFields(caseDirectory=casePath,
                                              dataset=_rectilinearDataset(),
                                              fieldMap=dict(T="temperature"),
                                              xCoordinate="x",
                                              yCoordinate="y",
                                              verticalCoordinate="z",
                                              datasetCRS=ITM,
                                              caseCRS=ITM)

    assert len(written) == 2
    for processorName, processorCentres in centres.items():
        temperature = FoamFieldFile(os.path.join(casePath, processorName, "0", "T"))
        numpy.testing.assert_allclose(
            numpy.atleast_1d(numpy.asarray(temperature.internal_field)),
            _linearField(processorCentres[:, 0], processorCentres[:, 1], processorCentres[:, 2]))


@pytest.mark.unit
def test_datasetToCaseFields_saysWhichFieldFileIsMissing(caseWithGeometry):
    with pytest.raises(FileNotFoundError, match="createEmptyCase"):
        datasetToOF.datasetToCaseFields(caseDirectory=caseWithGeometry["path"],
                                        dataset=_rectilinearDataset(),
                                        fieldMap=dict(T="temperature"),
                                        xCoordinate="x",
                                        yCoordinate="y",
                                        verticalCoordinate="z",
                                        datasetCRS=ITM,
                                        caseCRS=ITM)


# ---------------------------------------------------------------------------
#  datasetToSetFieldsDict
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_datasetToSetFieldsDict_boxesReachHalfwayToTheNeighbours(tmp_path):
    dataset = _rectilinearDataset()
    outputFile = str(tmp_path / "system" / "setFieldsDict")

    regions = datasetToOF.datasetToSetFieldsDict(dataset=dataset,
                                                 fieldMap=dict(T="temperature", U=("u", "v", 0)),
                                                 xCoordinate="x",
                                                 yCoordinate="y",
                                                 verticalCoordinate="z",
                                                 datasetCRS=ITM,
                                                 caseCRS=ITM,
                                                 defaultFieldValues=["volScalarFieldValue", "T", 300],
                                                 outputFile=outputFile)

    # One boxToCell and one boxToFace per cell of the dataset.
    assert len(regions) == 2 * dataset["temperature"].size
    assert [regionType for regionType, _ in regions[:2]] == ["boxToCell", "boxToFace"]

    # The cell at (x=180100, y=660000, z=50) is interior along x and z, so its
    # box spans halfway to its neighbours; along y it has one neighbour only,
    # so the spacing is mirrored.
    boxes = {tuple(numpy.round(numpy.concatenate(box["box"]), 3)) for regionType, box in regions
             if regionType == "boxToCell"}
    assert (180050.0, 659950.0, 30.0, 180150.0, 660050.0, 75.0) in boxes

    # Round-tripping through foamlib is the proof the file is valid OF syntax.
    written = FoamFile(outputFile)
    assert written["defaultFieldValues"] == ["volScalarFieldValue", "T", 300]
    assert len(written["regions"]) == len(regions)
    firstBox = written["regions"][0][1]
    assert firstBox["fieldValues"][:2] == ["volScalarFieldValue", "T"]
    assert "volVectorFieldValue" in firstBox["fieldValues"]


@pytest.mark.unit
def test_datasetToSetFieldsDict_extentDropsOutsideCellsAndClipsTheEdges(tmp_path):
    extent = dict(minX=179950.0, maxX=180150.0, minY=659950.0, maxY=660150.0, minZ=0.0, maxZ=60.0)

    regions = datasetToOF.datasetToSetFieldsDict(dataset=_rectilinearDataset(),
                                                 fieldMap=dict(T="temperature"),
                                                 xCoordinate="x",
                                                 yCoordinate="y",
                                                 verticalCoordinate="z",
                                                 datasetCRS=ITM,
                                                 caseCRS=ITM,
                                                 includeFaces=False)
    unclipped = len(regions)

    clipped = datasetToOF.datasetToSetFieldsDict(dataset=_rectilinearDataset(),
                                                 fieldMap=dict(T="temperature"),
                                                 xCoordinate="x",
                                                 yCoordinate="y",
                                                 verticalCoordinate="z",
                                                 datasetCRS=ITM,
                                                 caseCRS=ITM,
                                                 extent=extent,
                                                 includeFaces=False)

    # x=180200 and z=100 are outside: 2 of 3 x-columns and 2 of 3 levels survive.
    assert unclipped == 18
    assert len(clipped) == 2 * 2 * 2

    lowerLefts = numpy.array([box["box"][0] for _, box in clipped])
    upperRights = numpy.array([box["box"][1] for _, box in clipped])
    assert lowerLefts[:, 0].min() == pytest.approx(extent["minX"])
    assert lowerLefts[:, 1].min() == pytest.approx(extent["minY"])
    assert lowerLefts[:, 2].min() == pytest.approx(extent["minZ"])
    assert upperRights[:, 1].max() == pytest.approx(extent["maxY"])


@pytest.mark.unit
def test_datasetToSetFieldsDict_includeFacesControlsTheBoxToFaceRegions():
    cellsOnly = datasetToOF.datasetToSetFieldsDict(dataset=_rectilinearDataset(),
                                                   fieldMap=dict(T="temperature"),
                                                   xCoordinate="x",
                                                   yCoordinate="y",
                                                   verticalCoordinate="z",
                                                   datasetCRS=ITM,
                                                   caseCRS=ITM,
                                                   includeFaces=False)
    assert {regionType for regionType, _ in cellsOnly} == {"boxToCell"}


@pytest.mark.unit
def test_datasetToSetFieldsDict_projectsTheDatasetIntoTheCoordinatesOfTheCase():
    longitudes = numpy.array([34.7, 34.8])
    latitudes = numpy.array([31.9, 32.0])
    dataset = xarray.Dataset(
        data_vars=dict(temperature=(("z", "latitude", "longitude"), numpy.ones((2, 2, 2)) * 300.0)),
        coords=dict(longitude=longitudes, latitude=latitudes, z=numpy.array([10.0, 20.0])),
    )

    regions = datasetToOF.datasetToSetFieldsDict(dataset=dataset,
                                                 fieldMap=dict(T="temperature"),
                                                 xCoordinate="longitude",
                                                 yCoordinate="latitude",
                                                 verticalCoordinate="z",
                                                 datasetCRS=WSG84,
                                                 caseCRS=ITM,
                                                 includeFaces=False)

    boxCoordinates = numpy.array([box["box"][0] for _, box in regions])
    # ITM eastings/northings around Israel are 6-digit metres, not degrees.
    assert boxCoordinates[:, 0].min() > 100000.0
    assert boxCoordinates[:, 1].min() > 500000.0


@pytest.mark.unit
def test_datasetToSetFieldsDict_rejectsAnUnknownExtentLimit():
    with pytest.raises(ValueError, match="minx"):
        datasetToOF.datasetToSetFieldsDict(dataset=_rectilinearDataset(),
                                           fieldMap=dict(T="temperature"),
                                           xCoordinate="x",
                                           yCoordinate="y",
                                           verticalCoordinate="z",
                                           extent=dict(minx=0.0))


@pytest.mark.unit
def test_datasetToSetFieldsDict_aSingleCellAxisNeedsAnExtent():
    dataset = _rectilinearDataset().isel(z=[0])
    with pytest.raises(ValueError, match="single cell"):
        datasetToOF.datasetToSetFieldsDict(dataset=dataset,
                                           fieldMap=dict(T="temperature"),
                                           xCoordinate="x",
                                           yCoordinate="y",
                                           verticalCoordinate="z",
                                           datasetCRS=ITM,
                                           caseCRS=ITM)


# ---------------------------------------------------------------------------
#  The third checkbox of the issue
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_theModuleExposesNoPathsAndNoHardCodedProjection():
    source = open(datasetToOF.__file__).read()
    code = "\n".join(line for line in source.split("\n") if not line.strip().startswith(("#", "*", '"')))

    assert re.search(r"['\"]/(data|home|mnt|opt)\b", code) is None, "an absolute path leaked into the module"
    assert "proj=tmerc" not in code, "a proj4 projection string leaked into the module"
    assert re.search(r"EPSG:\d", code) is None, "an EPSG code is hard-coded instead of taken from an argument"
    assert ".npz" not in code and ".nc'" not in code, "an input file name leaked into the module"


# ---------------------------------------------------------------------------
#  The pre-existing xarrayToSetFieldsDictDomain
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_xarrayToSetFieldsDictDomain_producesTheDocumentedMapping():
    """The old generic converter, which four defects kept from ever running.

    ``docs/toolkits/simulations/openfoam.md`` documents it with
    ``U=("u10","v10",0), T="temperature"``, but as written it raised on every
    input: ``product`` was handed integers instead of ranges, a string field
    name fell into the vector branch, the vector branch overwrote the dataset
    with the list it was filling, and the vector formatter joined floats as if
    they were strings.  This test covers the mapping the docstring promises.

    Only the logger is taken from ``self``, so a stand-in object is enough --
    the toolkit itself needs a project and a database.
    """
    # Guarded here rather than at module level, because only this one test
    # reaches openFoam/toolkit.py -- the other 27 need datasetToOF alone, which
    # imports neither of these.  toolkit.py does `from evtk import hl` at import
    # time and pulls in VTKPipeline, which imports paraview; paraview is not
    # pip-installable (requirements.txt: 'make install-paraview'), so outside
    # the hermetic layer, which stubs it, this test skips.
    pytest.importorskip("evtk", reason="hera.simulations.openFoam.toolkit imports evtk at module level")
    pytest.importorskip("paraview", reason="openFoam.toolkit -> VTKPipeline imports paraview; not pip-installable")
    from hera.simulations.openFoam.toolkit import OFToolkit

    edges = numpy.array([0.0, 1.0, 2.0])
    dataset = xarray.Dataset(
        data_vars=dict(temperature=(("x", "y", "z"), numpy.arange(27.0).reshape(3, 3, 3)),
                       u=(("x", "y", "z"), numpy.ones((3, 3, 3))),
                       v=(("x", "y", "z"), numpy.full((3, 3, 3), 2.0))),
        coords=dict(x=edges, y=edges, z=edges),
    )

    class _LoggerHost:
        pass

    setFieldsText = OFToolkit.xarrayToSetFieldsDictDomain(_LoggerHost(),
                                                          xarrayData=dataset,
                                                          xColumnName="x",
                                                          yColumnName="y",
                                                          zColumnName="z",
                                                          U=("u", "v", 0),
                                                          T="temperature")

    # 2 cells along each of the three axes.
    assert setFieldsText.count("boxToCell") == 8
    assert "volVectorFieldValue U (1.0 2.0 0 )" in setFieldsText
    assert "volScalarFieldValue T 0.0" in setFieldsText
    assert "box (0.0 0.0 0.0 ) (1.0 1.0 1.0 );" in setFieldsText
