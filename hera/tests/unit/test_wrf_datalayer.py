"""simulations/WRF/wrfDatalayer.py: ``wrfDatalayer.find_i`` and
``wrfDatalayer.getPandas``.

``find_i`` converts a requested coordinate value (a longitude, say) into
a fractional index along the matching WRF dimension, so that xarray's
``interp`` can be pointed at it.  On a WRF grid the projected coordinate
varies along *both* horizontal dimensions, so each index ``i`` of
``west_east`` covers a *range* of longitudes; the function scans for the
index whose mean the request sits at or past and then interpolates
linearly between that column's mean and the next one's:

    find_i(request) = i + (request - mean(XLONG_i)) / (mean(XLONG_i+1) - mean(XLONG_i))

The tests build a WRF-shaped ``xarray.Dataset`` (dims ``Time``,
``bottom_top``, ``south_north``, ``west_east`` and their staggered twins)
with a deliberately tilted coordinate field, so the index the function
must return is arithmetic rather than guesswork.

Four defects surfaced:

* B246, fixed: ``delta`` and ``request_delta`` were assigned only inside
  the matching branch while the ``return`` statement used them
  unconditionally, and the branch tested the request against one column's
  *extent*.  A request falling between two columns' extents -- which on a
  tilted grid is most of the domain, since each column covers only its own
  spread -- exited the loop with the names unbound and raised
  ``UnboundLocalError``.  The bracket is now taken between consecutive
  column means, which do tile the axis, and a request outside the
  dimension altogether raises a ``ValueError`` that says so.
* B247, fixed: the scan read ``i + 1`` without a bounds check, so a
  request inside the *last* column of the dimension raised ``IndexError:
  index 4 is out of bounds for axis 2 with size 4`` -- the last row and
  column of every WRF domain were unaddressable.  The spacing for the last
  index is now read from the pair that ends at it.
* B248, fixed: ``import geopandas`` was commented out at the top of the
  module, but ``getPandas`` calls ``geopandas.GeoDataFrame`` and
  ``geopandas.points_from_xy`` in three places.  Passing an ITM
  coordinate -- the documented input, and the whole reason for the
  ``> 360`` branch -- raised ``NameError: name 'geopandas' is not
  defined``.  The final GeoDataFrame construction at the end of the
  function was dead for the same reason, so ``getPandas`` could not
  return.  The import is back.
* B249: ``import wrf`` is wrapped in ``try/except ImportError`` that
  only ``print``s.  With python-wrf absent the name is never bound, so
  ``getPandas`` dies with ``NameError: name 'wrf' is not defined``
  half-way through instead of the ``ImportError`` the message promises.

Not covered: ``getPandas``'s successful path.  It needs python-wrf
(``wrf.getvar``), which is not installed in this environment and is not
in requirements.txt, so those tests importorskip.
"""
import numpy
import pytest
import xarray

from hera.simulations.WRF.wrfDatalayer import wrfDatalayer

NI, NJ, NT, NZ = 4, 3, 2, 2


def _tilted(nRows, nCols, origin, rowStep, colStep):
    """A coordinate field that varies along both horizontal dimensions, the
    way a projected WRF coordinate does."""
    return numpy.array(
        [
            [origin + rowStep * row + colStep * col for col in range(nCols)]
            for row in range(nRows)
        ]
    )


def _grid():
    """A WRF-shaped dataset.

    XLONG spans 0.1 deg per west_east index and 0.02 deg per south_north
    index, so column i covers [34 + 0.1 i, 34 + 0.1 i + 0.04] with mean
    34 + 0.1 i + 0.02.
    """
    xlong = _tilted(NJ, NI, 34.0, 0.02, 0.1)
    xlongU = _tilted(NJ, NI + 1, 34.0, 0.02, 0.1)
    xlat = _tilted(NJ, NI, 32.0, 0.1, 0.02)
    xlatV = _tilted(NJ + 1, NI, 32.0, 0.1, 0.02)
    return xarray.Dataset(
        {
            "XLONG": (("Time", "south_north", "west_east"),
                      numpy.stack([xlong] * NT)),
            "XLONG_U": (("Time", "south_north", "west_east_stag"),
                        numpy.stack([xlongU] * NT)),
            "XLAT": (("Time", "south_north", "west_east"),
                     numpy.stack([xlat] * NT)),
            "XLAT_V": (("Time", "south_north_stag", "west_east"),
                       numpy.stack([xlatV] * NT)),
            "W": (("Time", "bottom_top", "south_north", "west_east"),
                  numpy.zeros((NT, NZ, NJ, NI))),
            "U": (("Time", "bottom_top", "south_north", "west_east_stag"),
                  numpy.zeros((NT, NZ, NJ, NI + 1))),
            "V": (("Time", "bottom_top", "south_north_stag", "west_east"),
                  numpy.zeros((NT, NZ, NJ + 1, NI))),
        }
    )


@pytest.fixture()
def datalayer():
    return wrfDatalayer()


@pytest.fixture()
def grid():
    return _grid()


@pytest.fixture()
def wrfFile(tmp_path):
    """The same grid on disk, since getPandas takes a path."""
    path = tmp_path / "wrfout_d01.nc"
    _grid().to_netcdf(path)
    return str(path)


@pytest.mark.unit
class TestFindIndex:
    def test_a_request_at_a_column_mean_returns_that_column(self, datalayer, grid):
        """mean(XLONG at west_east=1) is 34.12, so the offset term vanishes."""
        assert datalayer.find_i(34.12, grid, "west_east", "XLONG") == pytest.approx(1.0)

    def test_a_request_off_the_mean_interpolates_linearly(self, datalayer, grid):
        """34.13 is a tenth of the 0.1 deg column spacing past the mean of
        column 1, so the index is 1.1."""
        assert datalayer.find_i(34.13, grid, "west_east", "XLONG") == pytest.approx(1.1)

    def test_the_index_increases_with_the_request(self, datalayer, grid):
        near = datalayer.find_i(34.11, grid, "west_east", "XLONG")
        far = datalayer.find_i(34.21, grid, "west_east", "XLONG")
        assert far > near

    def test_consecutive_columns_are_one_index_apart(self, datalayer, grid):
        """The dimension is uniformly spaced, so equal coordinate steps give
        equal index steps."""
        first = datalayer.find_i(34.12, grid, "west_east", "XLONG")
        second = datalayer.find_i(34.22, grid, "west_east", "XLONG")
        assert second - first == pytest.approx(1.0)

    def test_it_works_on_the_latitude_dimension_too(self, datalayer, grid):
        """The dimension and coordinate names are parameters, so XLAT over
        south_north behaves identically. XLAT spans four west_east columns
        rather than three south_north rows, so row 1's mean is 32.13."""
        assert datalayer.find_i(32.13, grid, "south_north", "XLAT") == pytest.approx(
            1.0
        )

    def test_it_works_on_a_staggered_dimension(self, datalayer, grid):
        """getPandas asks for XLONG_U over west_east_stag for the U wind."""
        assert datalayer.find_i(
            34.12, grid, "west_east_stag", "XLONG_U"
        ) == pytest.approx(1.0)

    def test_the_returned_index_is_usable_for_interpolation(self, datalayer, grid):
        """Which is what the docstring promises it is for."""
        index = datalayer.find_i(34.13, grid, "west_east", "XLONG")
        interpolated = grid.isel(Time=0, bottom_top=0).interp(west_east=index)
        assert float(interpolated.XLONG.mean()) == pytest.approx(34.13)


@pytest.mark.unit
class TestFindIndexBetweenColumns:
    """B246: see the module docstring."""

    def test_a_request_inside_the_domain_resolves_to_an_index(self, datalayer, grid):
        """34.07 lies well inside the domain's longitude span
        [34.00, 34.34], so it has to map to some index."""
        index = datalayer.find_i(34.07, grid, "west_east", "XLONG")
        assert 0.0 <= index <= NI - 1

    def test_a_request_between_columns_interpolates_from_the_lower_mean(
        self, datalayer, grid
    ):
        """34.07 is half a column spacing past column 0's mean of 34.02."""
        assert datalayer.find_i(34.07, grid, "west_east", "XLONG") == pytest.approx(
            0.5
        )

    def test_a_request_outside_the_domain_is_reported_as_such(self, datalayer, grid):
        """B246: an out-of-domain request used to take the same unbound-name
        path as a gap; it now names the coordinate and the dimension."""
        with pytest.raises(ValueError, match="outside the west_east dimension"):
            datalayer.find_i(-120.0, grid, "west_east", "XLONG")

    def test_the_columns_really_do_leave_gaps(self, grid):
        """The premise, read off the data rather than the code: consecutive
        columns are 0.1 deg apart but each spans only 0.04 deg."""
        column0 = grid.isel(Time=0, south_north=slice(None), west_east=0).XLONG
        column1 = grid.isel(Time=0, south_north=slice(None), west_east=1).XLONG
        assert float(column0.max()) < 34.07 < float(column1.min())


@pytest.mark.unit
class TestFindIndexAtTheLastColumn:
    """B247: see the module docstring."""

    def test_the_last_column_is_addressable(self, datalayer, grid):
        """34.32 is the mean of the last column, squarely inside the domain."""
        index = datalayer.find_i(34.32, grid, "west_east", "XLONG")
        assert index == pytest.approx(NI - 1)

    def test_the_last_column_keeps_the_dimension_spacing(self, datalayer, grid):
        """B247: the spacing for the last index comes from the pair before
        it, so a request past its mean still advances by the same step."""
        assert datalayer.find_i(34.33, grid, "west_east", "XLONG") == pytest.approx(
            NI - 1 + 0.1
        )

    def test_the_request_really_is_inside_the_last_column(self, grid):
        """The premise: 34.32 is bracketed by the last column's own extents."""
        last = grid.isel(Time=0, west_east=NI - 1).XLONG
        assert float(last.min()) <= 34.32 <= float(last.max())


@pytest.mark.unit
class TestGetPandasArgumentHandling:
    def test_neither_latitude_nor_longitude_is_rejected(self, datalayer, wrfFile):
        """The docstring says either one of lat or lon must be given."""
        with pytest.raises(KeyError, match="Choose longtitude or latitude"):
            datalayer.getPandas(wrfFile)


@pytest.mark.unit
class TestGetPandasGeopandas:
    """B248: see the module docstring."""

    def test_an_itm_coordinate_no_longer_fails_on_geopandas(self, datalayer, wrfFile):
        """lat > 360 is the documented signal for an ITM coordinate.  The
        conversion itself runs now; whatever stops the call further on, it is
        no longer the missing name."""
        with pytest.raises(Exception) as raised:
            datalayer.getPandas(wrfFile, lat=650000.0)
        assert "geopandas" not in str(raised.value)

    def test_the_module_binds_geopandas(self):
        """B248's mechanism: the import was commented out.  Note the module
        has to be fetched from sys.modules -- the package's __init__ rebinds
        the name `wrfDatalayer` to the class."""
        import sys

        module = sys.modules["hera.simulations.WRF.wrfDatalayer"]
        assert hasattr(module, "geopandas")


@pytest.mark.unit
class TestGetPandasMissingWrf:
    """B249: see the module docstring."""

    @pytest.mark.xfail(
        strict=True,
        reason="B249: the module catches ImportError for `import wrf` and only "
               "prints 'You must install python-wrf to use this package', leaving "
               "the name unbound. getPandas then runs a long way in -- opening the "
               "file, resolving the horizontal indices -- before dying on "
               "`wrf.getvar(data, 'ter')` with NameError: name 'wrf' is not "
               "defined, which names neither the package nor the remedy. "
               "See the consolidated findings issue.",
    )
    def test_a_missing_dependency_is_reported_as_such(self, datalayer, wrfFile):
        with pytest.raises(ImportError, match="python-wrf"):
            datalayer.getPandas(wrfFile, lat=32.12)

    def test_it_currently_raises_a_bare_nameerror(self, datalayer, wrfFile):
        """Characterisation of B249, and proof the horizontal index lookup
        upstream of it succeeded."""
        pytest.importorskip(
            "netCDF4", reason="getPandas opens the file with netCDF4.Dataset"
        )
        if _wrfInstalled():
            pytest.skip("python-wrf is installed here, so the name is bound")
        with pytest.raises(NameError, match="wrf"):
            datalayer.getPandas(wrfFile, lat=32.12)

    def test_the_module_never_binds_wrf(self):
        """Characterisation of B249's mechanism."""
        import hera.simulations.WRF.wrfDatalayer as module

        if _wrfInstalled():
            pytest.skip("python-wrf is installed here, so the name is bound")
        assert not hasattr(module, "wrf")

    def test_the_documented_happy_path_needs_python_wrf(self, datalayer, wrfFile):
        """Skipped unless python-wrf is present; it is not installed here and
        is not declared in requirements.txt."""
        pytest.importorskip("wrf", reason="python-wrf is not installed")
        frame = datalayer.getPandas(wrfFile, lat=32.12)
        assert {"U", "V", "W", "height"} <= set(frame.columns)


def _wrfInstalled():
    import importlib.util

    return importlib.util.find_spec("wrf") is not None
