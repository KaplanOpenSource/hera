"""Inject a gridded meteorological dataset into an OpenFOAM case.

This module replaces the two scratch scripts that used to live in
``hera/simulations/openFoam/toberewritten`` -- ``netcdf2of.py`` and
``xarrayDataset2OF.py`` (issue #1060).  Both solved the same problem, namely
the meso-scale-to-CFD nesting step: take the output of an atmospheric model
(GFS/GDAS NetCDF, WRF, or anything else that loads as an
:class:`xarray.Dataset`) and use it as the initial and boundary condition of
an OpenFOAM case.

Two routes are supported, matching the two routes the old scripts took:

``datasetToCaseFields``
    Interpolate the dataset onto the cell centres (and, optionally, onto the
    boundary face centres) of an existing case and write the values straight
    into the field files.

``datasetToSetFieldsDict``
    Emit a ``setFieldsDict`` with one ``boxToCell``/``boxToFace`` region per
    dataset cell, to be applied afterwards by OpenFOAM's ``setFields``.

What changed relative to the old scripts
----------------------------------------
* OpenFOAM files are read and written with `foamlib
  <https://github.com/gerlero/foamlib>`_ instead of the hand-rolled
  tokeniser the scripts carried (which located a patch by name and then
  jumped a fixed number of whitespace-separated tokens, and stripped the
  vector parentheses with ``s[1:]``/``s[:-1]``).
* No paths, file names, projection strings or field names are baked in --
  every one of them is an argument.  The coordinate systems default to
  Hera's :data:`~hera.measurements.GIS.WSG84` and
  :data:`~hera.measurements.GIS.ITM` constants.
* The per-cell Python loop (``np.argmin`` over the whole grid for every
  single cell, printing its progress because it took hours) is replaced by
  vectorised horizontal interpolation and a vectorised vertical search.
* Coordinate axes are transformed with ``always_xy=True``, so ``x`` is
  always the easting and ``y`` always the northing.  The old scripts relied
  on the raw pyproj axis order and swapped the ``x``/``y`` names around it,
  which made the emitted ``box`` entries hard to reason about.

Notes
-----
``foamlib`` is imported lazily, so importing Hera does not require it; the
import error is raised only when one of these functions is actually called.
"""
import functools
import glob
import logging
import os

from hera.utils.lazy import _LazyModule
from hera.measurements.GIS import WSG84, ITM

# Heavy libraries are loaded on first attribute access, so that importing the
# openFOAM toolkit does not pay for numpy/xarray before they are needed.
numpy = _LazyModule("numpy")
xarray = _LazyModule("xarray")

#: Name of the field written by ``postProcess -func writeCellCentres``.  Its
#: ``internalField`` holds the cell centres and its ``boundaryField`` holds
#: the face centres of every patch, which is exactly the geometry needed
#: here -- and the reason the old scripts had to parse ``Cx``/``Cy``/``Cz``
#: by hand.
CELLCENTRES_FIELD = "C"

#: Key used for a case that is not decomposed.  Matches the convention of
#: :class:`~hera.simulations.openFoam.preprocessOFObjects.OFObject.OFObject`.
SINGLEPROCESSOR = "singleProcessor"

#: Key of the internal (volume) region in the geometry/values mappings.
INTERNALFIELD = "internalField"

HORIZONTAL_LINEAR = "linear"
HORIZONTAL_NEAREST = "nearest"


def _getFoamlib():
    """Import foamlib on demand.

    Returns
    -------
    tuple
        ``(FoamFile, FoamFieldFile)``.

    Raises
    ------
    ImportError
        If foamlib is not installed, with the install command in the message.
    """
    try:
        from foamlib import FoamFieldFile, FoamFile
    except ImportError as err:  # pragma: no cover - depends on the environment
        raise ImportError(
            "This function needs foamlib to read and write the OpenFOAM case. "
            "Install it with 'pip install foamlib'."
        ) from err
    return FoamFile, FoamFieldFile


@functools.lru_cache(maxsize=None)
def _getTransformer(inputCRS, outputCRS):
    """Build a pyproj transformer whose input and output are (x, y) ordered.

    Cached, because the same pair of coordinate systems is asked for once per
    region and per processor of a case.

    Parameters
    ----------
    inputCRS : int
        EPSG code of the input coordinates.
    outputCRS : int
        EPSG code of the output coordinates.

    Returns
    -------
    pyproj.Transformer
        A transformer with ``always_xy=True``, i.e. it takes and returns
        (easting, northing) / (longitude, latitude) -- never the
        authority-defined axis order.
    """
    from pyproj import Transformer

    return Transformer.from_crs(f"EPSG:{int(inputCRS)}", f"EPSG:{int(outputCRS)}", always_xy=True)


# ---------------------------------------------------------------------------
#  Reading the geometry of a case
# ---------------------------------------------------------------------------

def caseGeometry(caseDirectory, time=0, patchList=None, readParallel=True, cellCentresField=CELLCENTRES_FIELD):
    """Read the cell centres and the boundary face centres of a case.

    The geometry is read from the ``C`` field, which OpenFOAM writes with
    ``postProcess -func writeCellCentres``; see
    :meth:`~hera.simulations.openFoam.toolkit.OFToolkit.getMesh`, which runs
    that utility.

    Parameters
    ----------
    caseDirectory : str
        The case directory.
    time : str or int
        The time directory the cell centres were written to.
    patchList : list of str, optional
        The patches whose face centres are needed.  ``None`` (the default)
        returns every patch that carries per-face values; an empty list
        returns the internal field only.
    readParallel : bool
        If the case is decomposed and this is True, read every
        ``processor*`` directory.  Otherwise read the reconstructed case.
    cellCentresField : str
        Name of the cell-centres field file.

    Returns
    -------
    dict
        ``{processorName: {"internalField": (N,3) array, patchName: (M,3) array}}``.
        ``processorName`` is :data:`SINGLEPROCESSOR` for a case that is not
        decomposed, and ``processor0``, ``processor1``, ... for one that is.

    Notes
    -----
    A patch whose ``C`` value is *uniform* carries a single vector rather
    than one per face, so there is no way to tell how many faces it has or
    where they are.  Such a patch is skipped (with a warning) instead of
    being silently filled with a constant, which is what the old
    ``netcdf2of.py`` did for its ``top`` patch.
    """
    logger = logging.getLogger(f"{__name__}.caseGeometry")
    _, FoamFieldFile = _getFoamlib()

    processorDirectories = {SINGLEPROCESSOR: caseDirectory}
    if readParallel:
        decomposed = sorted(glob.glob(os.path.join(caseDirectory, "processor*")))
        if len(decomposed) > 0:
            logger.debug(f"Found a decomposed case with {len(decomposed)} processors")
            processorDirectories = dict([(os.path.basename(procDir), procDir) for procDir in decomposed])

    ret = dict()
    for processorName, processorDirectory in processorDirectories.items():
        centresPath = os.path.join(processorDirectory, str(time), cellCentresField)
        if not os.path.isfile(centresPath):
            raise FileNotFoundError(
                f"The cell centres file {centresPath} does not exist. Write it with "
                f"'postProcess -func writeCellCentres' (or call OFToolkit.getMesh, which does that for you)."
            )

        centres = FoamFieldFile(centresPath)
        geometry = {INTERNALFIELD: numpy.atleast_2d(numpy.asarray(centres.internal_field, dtype=float))}

        requestedPatches = list(centres.boundary_field) if patchList is None else list(patchList)
        for patchName in requestedPatches:
            if patchName not in centres.boundary_field:
                raise KeyError(f"Patch '{patchName}' is not in {centresPath}. Found: {list(centres.boundary_field)}.")
            patchValue = numpy.asarray(centres.boundary_field[patchName]["value"], dtype=float)
            if patchValue.ndim == 1:
                logger.warning(
                    f"Patch '{patchName}' has a uniform value in {centresPath}, so its per-face centres are "
                    f"unknown; skipping it. Recompute the cell centres for this patch to include it."
                )
                continue
            geometry[patchName] = patchValue

        ret[processorName] = geometry
        logger.debug(f"{processorName}: {[(region, values.shape) for region, values in geometry.items()]}")

    return ret


# ---------------------------------------------------------------------------
#  Interpolating a dataset onto arbitrary points
# ---------------------------------------------------------------------------

def _fieldMapVariables(fieldMap):
    """List the dataset variables a field map refers to.

    Parameters
    ----------
    fieldMap : dict
        See :func:`interpolateDatasetToPoints`.

    Returns
    -------
    list of str
        The variable names, without repetitions and in order of appearance.
    """
    variableNames = []
    for openFOAMField, mapping in fieldMap.items():
        components = mapping if isinstance(mapping, (tuple, list)) else [mapping]
        if len(components) not in [1, 3, 9]:
            raise ValueError(
                f"The mapping of {openFOAMField} has {len(components)} components; it must have 1 (scalar), "
                f"3 (vector) or 9 (tensor)."
            )
        for component in components:
            if isinstance(component, str):
                if component not in variableNames:
                    variableNames.append(component)
            elif not isinstance(component, (int, float)):
                raise ValueError(
                    f"The mapping of {openFOAMField} contains {component!r}, which is neither a dataset "
                    f"variable name nor a number."
                )
    return variableNames


def _horizontalGrid(dataset, xCoordinate, yCoordinate):
    """Describe the horizontal grid of a dataset.

    Parameters
    ----------
    dataset : xarray.Dataset
        The dataset.
    xCoordinate, yCoordinate : str
        Names of the horizontal coordinates.  Either both 1D (a rectilinear
        grid) or both 2D (a curvilinear one).  They may be coordinates or
        plain data variables -- WRF's ``XLONG``/``XLAT`` are usually the
        latter.

    Returns
    -------
    dict
        ``x``/``y``: the coordinate values as 2D grids, in the order of
        ``dimensions``; ``dimensions``: the names of the two horizontal
        dimensions, matching the axes of those grids; ``curvilinear``:
        whether the original coordinates were 2D.
    """
    if xCoordinate not in dataset:
        raise KeyError(f"'{xCoordinate}' is neither a coordinate nor a variable of the dataset.")
    if yCoordinate not in dataset:
        raise KeyError(f"'{yCoordinate}' is neither a coordinate nor a variable of the dataset.")

    xValues = numpy.asarray(dataset[xCoordinate].values, dtype=float)
    yValues = numpy.asarray(dataset[yCoordinate].values, dtype=float)

    if xValues.ndim == 1 and yValues.ndim == 1:
        yGrid, xGrid = numpy.meshgrid(yValues, xValues, indexing="ij")
        dimensions = [str(dataset[yCoordinate].dims[0]), str(dataset[xCoordinate].dims[0])]
        return dict(x=xGrid, y=yGrid, dimensions=dimensions, curvilinear=False)

    if xValues.ndim == 2 and yValues.ndim == 2:
        if xValues.shape != yValues.shape:
            raise ValueError(
                f"'{xCoordinate}' and '{yCoordinate}' have different shapes, {xValues.shape} and {yValues.shape}."
            )
        return dict(x=xValues, y=yValues,
                    dimensions=[str(dimension) for dimension in dataset[xCoordinate].dims],
                    curvilinear=True)

    raise ValueError(
        f"'{xCoordinate}' and '{yCoordinate}' must both be 1D (rectilinear grid) or both 2D (curvilinear); "
        f"got {xValues.ndim}D and {yValues.ndim}D."
    )


def _horizontalInterpolation(dataArray, pointX, pointY, xCoordinate, yCoordinate, grid, verticalDimension, method):
    """Interpolate a data array horizontally onto a list of points.

    Parameters
    ----------
    dataArray : xarray.DataArray
        The array to interpolate.  May or may not have a vertical dimension.
    pointX, pointY : numpy.ndarray
        The (N,) target coordinates, in the coordinate system of the dataset.
    xCoordinate, yCoordinate : str
        Names of the horizontal coordinates of the dataset.
    grid : dict
        The output of :func:`_horizontalGrid`.
    verticalDimension : str or None
        Name of the vertical dimension, if the array has one.
    method : str
        :data:`HORIZONTAL_LINEAR` (bilinear; needs 1D horizontal
        coordinates) or :data:`HORIZONTAL_NEAREST` (nearest column; works
        for a curvilinear grid whose coordinates are 2D as well).

    Returns
    -------
    numpy.ndarray
        Shape ``(nLevels, N)``; ``nLevels`` is 1 for an array with no
        vertical dimension.
    """
    if method == HORIZONTAL_LINEAR:
        if grid["curvilinear"]:
            raise ValueError(
                f"Bilinear interpolation needs 1D horizontal coordinates, but '{xCoordinate}' and "
                f"'{yCoordinate}' are 2D (a curvilinear grid). Either regrid the dataset to a rectilinear "
                f"grid, or pass horizontalMethod='{HORIZONTAL_NEAREST}'."
            )
        interpolated = dataArray.interp({
            xCoordinate: xarray.DataArray(pointX, dims="point"),
            yCoordinate: xarray.DataArray(pointY, dims="point"),
        })
        if verticalDimension is None:
            return numpy.asarray(interpolated.values, dtype=float).reshape(1, -1)
        return numpy.asarray(interpolated.transpose(verticalDimension, "point").values, dtype=float)

    if method != HORIZONTAL_NEAREST:
        raise ValueError(f"horizontalMethod must be '{HORIZONTAL_LINEAR}' or '{HORIZONTAL_NEAREST}', got {method!r}.")

    # Nearest column: works for a rectilinear and for a curvilinear grid alike.
    from scipy.spatial import cKDTree

    tree = cKDTree(numpy.column_stack([grid["x"].ravel(), grid["y"].ravel()]))
    _, flatIndex = tree.query(numpy.column_stack([pointX, pointY]))

    dimensionOrder = ([verticalDimension] if verticalDimension is not None else []) + list(grid["dimensions"])
    values = numpy.asarray(dataArray.transpose(*dimensionOrder).values, dtype=float)
    if verticalDimension is None:
        return values.reshape(1, -1)[:, flatIndex]
    return values.reshape(values.shape[0], -1)[:, flatIndex]


def _verticalInterpolation(columnValues, columnHeights, targetHeight):
    """Linearly interpolate columns of values to a height, per point.

    Parameters
    ----------
    columnValues : numpy.ndarray
        ``(nLevels, N)`` values.
    columnHeights : numpy.ndarray
        ``(nLevels, N)`` heights of those values, monotonic along axis 0.
    targetHeight : numpy.ndarray
        ``(N,)`` requested heights.

    Returns
    -------
    numpy.ndarray
        ``(N,)`` interpolated values.

    Notes
    -----
    A point below the lowest level or above the highest one takes the value
    of that level: the meteorological profile is not extrapolated (below
    ground level it has no physical meaning anyway).
    """
    logger = logging.getLogger(f"{__name__}._verticalInterpolation")
    if columnValues.shape != columnHeights.shape:
        raise ValueError(
            f"The values and their heights have different shapes, {columnValues.shape} and "
            f"{columnHeights.shape}: the variable and the height variable are on different vertical grids."
        )

    nLevels = columnValues.shape[0]
    if nLevels == 1:
        return columnValues[0, :]

    differences = numpy.diff(columnHeights, axis=0)
    if numpy.all(differences > 0):
        heights, values = columnHeights, columnValues
    elif numpy.all(differences < 0):
        # e.g. levels ordered by descending pressure: flip to ascending height.
        heights, values = columnHeights[::-1, :], columnValues[::-1, :]
    else:
        raise ValueError(
            "The heights of the vertical levels are not monotonic along every column, so a level cannot be "
            "located unambiguously. Check the vertical coordinate / heightVariable of the dataset."
        )

    pointIndex = numpy.arange(columnValues.shape[1])
    upper = numpy.clip(numpy.sum(heights < targetHeight[None, :], axis=0), 1, nLevels - 1)
    lower = upper - 1

    lowerHeight = heights[lower, pointIndex]
    upperHeight = heights[upper, pointIndex]
    thickness = upperHeight - lowerHeight
    weight = numpy.where(thickness == 0, 0.0, (targetHeight - lowerHeight) / numpy.where(thickness == 0, 1.0, thickness))

    outside = int(numpy.count_nonzero((weight < 0) | (weight > 1)))
    if outside > 0:
        logger.warning(
            f"{outside} of {columnValues.shape[1]} points are outside the vertical extent of the dataset; "
            f"they take the value of the nearest level (no extrapolation)."
        )
    weight = numpy.clip(weight, 0.0, 1.0)

    return values[lower, pointIndex] * (1 - weight) + values[upper, pointIndex] * weight


def interpolateDatasetToPoints(dataset,
                               points,
                               fieldMap,
                               xCoordinate,
                               yCoordinate,
                               verticalCoordinate=None,
                               heightVariable=None,
                               time=None,
                               timeCoordinate="time",
                               datasetCRS=WSG84,
                               caseCRS=ITM,
                               horizontalMethod=HORIZONTAL_LINEAR):
    """Interpolate dataset variables onto a list of points.

    Parameters
    ----------
    dataset : xarray.Dataset
        The meteorological data.
    points : numpy.ndarray or pandas.DataFrame
        ``(N,3)`` points -- (x, y, z) in the coordinate system of the case.
        z is the height above sea level, in metres.
    fieldMap : dict
        Maps an OpenFOAM field name to the dataset variable(s) that fill it.
        A string is a scalar; a tuple/list of 3 (9) makes a vector (tensor).
        Each component is either the name of a dataset variable or a fixed
        number.  For example ``dict(U=("ugrd", "vgrd", 0), T="Temp")``.
    xCoordinate, yCoordinate : str
        Names of the horizontal coordinates of the dataset (for example
        ``"grid_xt"``/``"grid_yt"``, or ``"XLONG"``/``"XLAT"``).
    verticalCoordinate : str, optional
        Name of the vertical dimension.  ``None`` for a dataset with no
        vertical dimension.
    heightVariable : str, optional
        A dataset variable holding the height above sea level [m] of every
        level, for a terrain-following vertical coordinate (WRF's
        ``gpt_hgt_M``, for instance).  If ``None``, the values of
        ``verticalCoordinate`` are taken to be those heights.
    time : optional
        If given, the nearest time step is selected first.
    timeCoordinate : str
        Name of the time coordinate.
    datasetCRS : int
        EPSG code of the horizontal coordinates of the dataset.
    caseCRS : int
        EPSG code of the coordinates of the case (and hence of ``points``).
    horizontalMethod : str
        :data:`HORIZONTAL_LINEAR` or :data:`HORIZONTAL_NEAREST`.

    Returns
    -------
    dict
        ``{openFOAMFieldName: array}``, ``(N,)`` for a scalar and
        ``(N, components)`` for a vector or a tensor.
    """
    logger = logging.getLogger(f"{__name__}.interpolateDatasetToPoints")

    pointArray = numpy.asarray(points, dtype=float)
    if pointArray.ndim != 2 or pointArray.shape[1] != 3:
        raise ValueError(f"points must be an (N,3) array of (x,y,z); got shape {pointArray.shape}.")

    if time is not None:
        if timeCoordinate not in dataset.coords and timeCoordinate not in dataset.dims:
            raise KeyError(f"The dataset has no '{timeCoordinate}' coordinate to select a time from.")
        dataset = dataset.sel({timeCoordinate: time}, method="nearest")
        logger.debug(f"Selected the time step nearest to {time}")

    # The case is in a projected CRS (metres); the dataset is usually in
    # lon/lat. always_xy=True keeps 'x' the easting/longitude throughout.
    transformer = _getTransformer(caseCRS, datasetCRS)
    pointX, pointY = transformer.transform(pointArray[:, 0], pointArray[:, 1])
    pointZ = pointArray[:, 2]

    variableNames = _fieldMapVariables(fieldMap)
    logger.debug(f"Interpolating {variableNames} onto {pointArray.shape[0]} points ({horizontalMethod})")

    # xarray's interp needs an ascending coordinate; pressure levels and
    # some latitude axes are stored descending.
    if horizontalMethod == HORIZONTAL_LINEAR:
        sortDimensions = [coordinate for coordinate in (xCoordinate, yCoordinate)
                          if coordinate in dataset.dims and dataset[coordinate].values.ndim == 1]
        if len(sortDimensions) > 0:
            dataset = dataset.sortby(sortDimensions)

    grid = _horizontalGrid(dataset, xCoordinate, yCoordinate)

    columns = dict()
    for variableName in variableNames:
        if variableName not in dataset:
            raise KeyError(f"The variable '{variableName}' is not in the dataset. Found: {list(dataset.data_vars)}.")
        variable = dataset[variableName]
        variableVertical = verticalCoordinate if (verticalCoordinate is not None
                                                  and verticalCoordinate in variable.dims) else None
        columns[variableName] = _horizontalInterpolation(variable, pointX, pointY, xCoordinate, yCoordinate,
                                                         grid, variableVertical, horizontalMethod)

    # The heights of the levels, per point.
    if verticalCoordinate is None:
        columnHeights = None
    elif heightVariable is not None:
        if heightVariable not in dataset:
            raise KeyError(f"The height variable '{heightVariable}' is not in the dataset.")
        columnHeights = _horizontalInterpolation(dataset[heightVariable], pointX, pointY, xCoordinate, yCoordinate,
                                                 grid, verticalCoordinate, horizontalMethod)
    else:
        levelHeights = numpy.asarray(dataset[verticalCoordinate].values, dtype=float)
        columnHeights = numpy.repeat(levelHeights[:, None], pointArray.shape[0], axis=1)

    ret = dict()
    for openFOAMField, mapping in fieldMap.items():
        components = list(mapping) if isinstance(mapping, (tuple, list)) else [mapping]
        componentValues = []
        for component in components:
            if isinstance(component, str):
                columnValues = columns[component]
                if columnValues.shape[0] == 1 or columnHeights is None:
                    componentValues.append(columnValues[0, :])
                else:
                    componentValues.append(_verticalInterpolation(columnValues, columnHeights, pointZ))
            else:
                componentValues.append(numpy.full(pointArray.shape[0], float(component)))

        ret[openFOAMField] = componentValues[0] if len(componentValues) == 1 \
            else numpy.column_stack(componentValues)

    return ret


# ---------------------------------------------------------------------------
#  Route 1: write the values straight into the field files of a case
# ---------------------------------------------------------------------------

def datasetToCaseFields(caseDirectory,
                        dataset,
                        fieldMap,
                        xCoordinate,
                        yCoordinate,
                        verticalCoordinate=None,
                        heightVariable=None,
                        time=None,
                        timeCoordinate="time",
                        datasetCRS=WSG84,
                        caseCRS=ITM,
                        horizontalMethod=HORIZONTAL_LINEAR,
                        writeTime="0",
                        geometryTime=0,
                        patchList=None,
                        readParallel=True):
    """Interpolate a dataset onto a case and write it into the field files.

    This is the modern replacement of ``netcdf2of.py``: for every field in
    ``fieldMap`` the ``internalField`` and the requested patches of
    ``<case>/<writeTime>/<field>`` are overwritten with the interpolated
    values.

    Parameters
    ----------
    caseDirectory : str
        The case to write into.  The field files must already exist -- they
        define the type, the dimensions and the boundary conditions.  Create
        them with
        :meth:`~hera.simulations.openFoam.toolkit.OFToolkit.createEmptyCase`
        if needed.
    dataset : xarray.Dataset
        The meteorological data.
    fieldMap : dict
        OpenFOAM field name -> dataset variable(s); see
        :func:`interpolateDatasetToPoints`.
    xCoordinate, yCoordinate, verticalCoordinate, heightVariable, time, timeCoordinate, datasetCRS, caseCRS, horizontalMethod
        See :func:`interpolateDatasetToPoints`.
    writeTime : str
        The time directory of the field files to write into.
    geometryTime : str or int
        The time directory the cell centres were written to.
    patchList : list of str, optional
        Patches to set as well as the internal field.  ``None`` sets every
        patch that has per-face cell centres and appears in the field file;
        pass ``[]`` for the internal field only.
    readParallel : bool
        Write into every ``processor*`` directory of a decomposed case.

    Returns
    -------
    list of str
        The paths of the field files that were written.
    """
    logger = logging.getLogger(f"{__name__}.datasetToCaseFields")
    _, FoamFieldFile = _getFoamlib()

    geometry = caseGeometry(caseDirectory,
                            time=geometryTime,
                            patchList=patchList,
                            readParallel=readParallel)

    writtenFiles = []
    for processorName, processorGeometry in geometry.items():
        processorDirectory = caseDirectory if processorName == SINGLEPROCESSOR \
            else os.path.join(caseDirectory, processorName)

        interpolatedRegions = dict()
        for regionName, regionPoints in processorGeometry.items():
            interpolatedRegions[regionName] = interpolateDatasetToPoints(
                dataset=dataset,
                points=regionPoints,
                fieldMap=fieldMap,
                xCoordinate=xCoordinate,
                yCoordinate=yCoordinate,
                verticalCoordinate=verticalCoordinate,
                heightVariable=heightVariable,
                time=time,
                timeCoordinate=timeCoordinate,
                datasetCRS=datasetCRS,
                caseCRS=caseCRS,
                horizontalMethod=horizontalMethod)

        for fieldName in fieldMap.keys():
            fieldPath = os.path.join(processorDirectory, str(writeTime), fieldName)
            if not os.path.isfile(fieldPath):
                raise FileNotFoundError(
                    f"The field file {fieldPath} does not exist. Create the case fields first (for example with "
                    f"OFToolkit.createEmptyCase / writeEmptyField) and then set their values."
                )

            fieldFile = FoamFieldFile(fieldPath)
            with fieldFile:
                fieldFile.internal_field = interpolatedRegions[INTERNALFIELD][fieldName]
                for regionName, regionValues in interpolatedRegions.items():
                    if regionName == INTERNALFIELD:
                        continue
                    if regionName not in fieldFile.boundary_field:
                        logger.warning(f"Patch '{regionName}' is not in {fieldPath}; not setting it.")
                        continue
                    fieldFile.boundary_field[regionName]["value"] = regionValues[fieldName]

            writtenFiles.append(fieldPath)
            logger.info(f"Wrote {fieldPath}")

    return writtenFiles


# ---------------------------------------------------------------------------
#  Route 2: emit a setFieldsDict
# ---------------------------------------------------------------------------

def _cellEdges(centres, axis, lowerLimit=None, upperLimit=None):
    """Return the lower and the upper edge of every cell along one axis.

    An edge between two cells is placed halfway between their centres --
    which is what the old scripts computed with their ``dxm``/``dxp``
    arithmetic.  The outermost edges are placed at ``lowerLimit`` /
    ``upperLimit`` when those are given, and otherwise mirror the spacing of
    the adjacent cell.

    Parameters
    ----------
    centres : numpy.ndarray
        The cell centres; may have any number of dimensions.
    axis : int
        The axis to compute the edges along.
    lowerLimit, upperLimit : float, optional
        The position of the outermost edges.

    Returns
    -------
    tuple of numpy.ndarray
        ``(lower, upper)``, both shaped like ``centres``.
    """
    centres = numpy.asarray(centres, dtype=float)
    if centres.shape[axis] < 2:
        if lowerLimit is None or upperLimit is None:
            raise ValueError(
                f"Axis {axis} has a single cell, so its extent cannot be derived from the spacing of its "
                f"neighbours. Pass an 'extent' to define it."
            )
        return (numpy.full_like(centres, float(lowerLimit)), numpy.full_like(centres, float(upperLimit)))

    midPoints = 0.5 * (numpy.take(centres, numpy.arange(centres.shape[axis] - 1), axis=axis) +
                       numpy.take(centres, numpy.arange(1, centres.shape[axis]), axis=axis))

    firstCentre = numpy.take(centres, [0], axis=axis)
    firstMid = numpy.take(midPoints, [0], axis=axis)
    lastCentre = numpy.take(centres, [centres.shape[axis] - 1], axis=axis)
    lastMid = numpy.take(midPoints, [midPoints.shape[axis] - 1], axis=axis)

    lowerEdge = firstCentre - (firstMid - firstCentre) if lowerLimit is None \
        else numpy.full_like(firstCentre, float(lowerLimit))
    upperEdge = lastCentre + (lastCentre - lastMid) if upperLimit is None \
        else numpy.full_like(lastCentre, float(upperLimit))

    lower = numpy.concatenate([lowerEdge, midPoints], axis=axis)
    upper = numpy.concatenate([midPoints, upperEdge], axis=axis)
    return lower, upper


def _fieldValueTokens(fieldMap, dataArrays, index):
    """Build the ``fieldValues`` token list of one region.

    Parameters
    ----------
    fieldMap : dict
        OpenFOAM field name -> dataset variable(s) or number(s).
    dataArrays : dict
        variable name -> ``(nz, ny, nx)`` numpy array.
    index : tuple
        The ``(iz, iy, ix)`` index of the cell.

    Returns
    -------
    list
        Tokens in the order OpenFOAM expects, e.g.
        ``["volScalarFieldValue", "T", 300.0, "volVectorFieldValue", "U", [1.0, 2.0, 0.0]]``.
    """
    tokens = []
    for openFOAMField, mapping in fieldMap.items():
        components = list(mapping) if isinstance(mapping, (tuple, list)) else [mapping]
        values = [float(dataArrays[component][index]) if isinstance(component, str) else float(component)
                  for component in components]
        if len(values) == 1:
            tokens += ["volScalarFieldValue", openFOAMField, values[0]]
        elif len(values) == 3:
            tokens += ["volVectorFieldValue", openFOAMField, values]
        else:
            tokens += ["volTensorFieldValue", openFOAMField, values]
    return tokens


def datasetToSetFieldsDict(dataset,
                           fieldMap,
                           xCoordinate,
                           yCoordinate,
                           verticalCoordinate,
                           heightVariable=None,
                           time=None,
                           timeCoordinate="time",
                           datasetCRS=WSG84,
                           caseCRS=ITM,
                           extent=None,
                           defaultFieldValues=None,
                           includeFaces=True,
                           outputFile=None):
    """Build a ``setFieldsDict`` with one region per cell of the dataset.

    This is the modern replacement of the ``setFieldsDict`` half of
    ``xarrayDataset2OF.py``/``wrf2of.py``: every dataset cell becomes a
    ``boxToCell`` (and, by default, a ``boxToFace``) region whose box spans
    halfway to the neighbouring cells, carrying the values of the mapped
    fields.  The dictionary is written by foamlib, so its header and its
    syntax are not assembled by hand.

    Parameters
    ----------
    dataset : xarray.Dataset
        The meteorological data.
    fieldMap : dict
        OpenFOAM field name -> dataset variable(s); see
        :func:`interpolateDatasetToPoints`.
    xCoordinate, yCoordinate : str
        Names of the horizontal coordinates of the dataset.  They may be 1D
        (a rectilinear grid) or 2D (a curvilinear one).
    verticalCoordinate : str
        Name of the vertical dimension.
    heightVariable : str, optional
        Height above sea level [m] of every level; see
        :func:`interpolateDatasetToPoints`.
    time, timeCoordinate, datasetCRS, caseCRS
        See :func:`interpolateDatasetToPoints`.  The horizontal coordinates
        are transformed from ``datasetCRS`` to ``caseCRS``, because the box
        coordinates must be in the coordinate system of the case.
    extent : dict, optional
        ``dict(minX=, maxX=, minY=, maxY=, minZ=, maxZ=)`` -- cells whose
        centre falls outside are dropped, and the outermost boxes are
        clipped to it.  Any subset of the six keys may be given.
    defaultFieldValues : list, optional
        Tokens of the ``defaultFieldValues`` entry, e.g.
        ``["volScalarFieldValue", "T", 300]``.  Omitted from the dictionary
        if ``None``.
    includeFaces : bool
        Emit a ``boxToFace`` region beside every ``boxToCell`` one, so that
        boundary faces are set as well.
    outputFile : str, optional
        Where to write the dictionary (typically
        ``<case>/system/setFieldsDict``).  If ``None``, nothing is written.

    Returns
    -------
    list
        The ``regions`` entry as foamlib represents it: a list of
        ``(regionType, dict)`` pairs.
    """
    logger = logging.getLogger(f"{__name__}.datasetToSetFieldsDict")
    FoamFile, _ = _getFoamlib()

    if time is not None:
        if timeCoordinate not in dataset.coords and timeCoordinate not in dataset.dims:
            raise KeyError(f"The dataset has no '{timeCoordinate}' coordinate to select a time from.")
        dataset = dataset.sel({timeCoordinate: time}, method="nearest")

    extent = dict() if extent is None else extent
    unknownLimits = set(extent) - {"minX", "maxX", "minY", "maxY", "minZ", "maxZ"}
    if len(unknownLimits) > 0:
        raise ValueError(f"Unknown extent limits {sorted(unknownLimits)}; use minX/maxX/minY/maxY/minZ/maxZ.")

    if verticalCoordinate not in dataset.dims:
        raise KeyError(f"The dataset has no '{verticalCoordinate}' dimension. Found: {list(dataset.dims)}.")

    # --- the (x, y) of every column, in the coordinate system of the case.
    grid = _horizontalGrid(dataset, xCoordinate, yCoordinate)
    horizontalDimensions = grid["dimensions"]

    transformer = _getTransformer(datasetCRS, caseCRS)
    caseX, caseY = transformer.transform(grid["x"], grid["y"])
    caseX = numpy.asarray(caseX, dtype=float)
    caseY = numpy.asarray(caseY, dtype=float)

    dimensionOrder = [verticalCoordinate] + horizontalDimensions
    nLevels = int(dataset.sizes[verticalCoordinate])

    # --- the height of every cell centre, shaped (nz, ny, nx).
    if heightVariable is not None:
        if heightVariable not in dataset:
            raise KeyError(f"The height variable '{heightVariable}' is not in the dataset.")
        caseZ = numpy.asarray(dataset[heightVariable].transpose(*dimensionOrder).values, dtype=float)
    else:
        levelHeights = numpy.asarray(dataset[verticalCoordinate].values, dtype=float)
        caseZ = numpy.broadcast_to(levelHeights[:, None, None], (nLevels,) + caseX.shape).copy()

    dataArrays = dict()
    for variableName in _fieldMapVariables(fieldMap):
        if variableName not in dataset:
            raise KeyError(f"The variable '{variableName}' is not in the dataset. Found: {list(dataset.data_vars)}.")
        variable = dataset[variableName]
        if verticalCoordinate in variable.dims:
            values = numpy.asarray(variable.transpose(*dimensionOrder).values, dtype=float)
        else:
            values = numpy.broadcast_to(
                numpy.asarray(variable.transpose(*horizontalDimensions).values, dtype=float),
                (nLevels,) + caseX.shape).copy()
        dataArrays[variableName] = values

    # --- the box of every cell: halfway to its neighbours along each axis.
    xLower, xUpper = _cellEdges(caseX, axis=1, lowerLimit=extent.get("minX"), upperLimit=extent.get("maxX"))
    yLower, yUpper = _cellEdges(caseY, axis=0, lowerLimit=extent.get("minY"), upperLimit=extent.get("maxY"))
    zLower, zUpper = _cellEdges(caseZ, axis=0, lowerLimit=extent.get("minZ"), upperLimit=extent.get("maxZ"))

    inside = numpy.ones(caseZ.shape, dtype=bool)
    for limitName, coordinateValues, comparison in [("minX", caseX, numpy.greater_equal),
                                                    ("maxX", caseX, numpy.less_equal),
                                                    ("minY", caseY, numpy.greater_equal),
                                                    ("maxY", caseY, numpy.less_equal),
                                                    ("minZ", caseZ, numpy.greater_equal),
                                                    ("maxZ", caseZ, numpy.less_equal)]:
        if limitName in extent:
            mask = comparison(coordinateValues, float(extent[limitName]))
            inside &= mask if mask.shape == inside.shape else mask[None, :, :]

    regions = []
    for levelIndex, rowIndex, columnIndex in zip(*numpy.nonzero(inside)):
        index = (levelIndex, rowIndex, columnIndex)
        lowerLeft = [float(xLower[rowIndex, columnIndex]),
                     float(yLower[rowIndex, columnIndex]),
                     float(zLower[index])]
        upperRight = [float(xUpper[rowIndex, columnIndex]),
                      float(yUpper[rowIndex, columnIndex]),
                      float(zUpper[index])]
        box = dict(box=(lowerLeft, upperRight), fieldValues=_fieldValueTokens(fieldMap, dataArrays, index))
        regions.append(("boxToCell", box))
        if includeFaces:
            regions.append(("boxToFace", dict(box)))

    logger.info(f"Built {len(regions)} regions out of {caseZ.size} dataset cells")

    if outputFile is not None:
        outputDirectory = os.path.dirname(os.path.abspath(outputFile))
        if len(outputDirectory) > 0:
            os.makedirs(outputDirectory, exist_ok=True)
        if os.path.exists(outputFile):
            os.remove(outputFile)
        setFieldsDict = FoamFile(outputFile)
        with setFieldsDict:
            setFieldsDict["FoamFile"] = dict(version=2.0, format="ascii", **{"class": "dictionary"},
                                             location='"system"', object="setFieldsDict")
            if defaultFieldValues is not None:
                setFieldsDict["defaultFieldValues"] = list(defaultFieldValues)
            setFieldsDict["regions"] = regions
        logger.info(f"Wrote {outputFile}")

    return regions
