"""The raster / archive / class ``DataHandler_*`` families in
``datalayer/datahandler.py``.

Every handler in that module is a pair of ``@staticmethod``s -- ``saveData``
(object -> file, returning a dict of *store parameters* that
``MetadataFrame.getData`` later feeds back in as kwargs) and ``getData``
(resource -> object).  This file covers the five handlers whose formats are
either write-only-in-theory or need an optional binary dependency:

* ``geotiff`` -- ``getData`` goes through GDAL, ``saveData`` is a stub.
* ``tif``     -- ``getData`` goes through rasterio, ``saveData`` is a stub.
* ``HDF``     -- ``getData`` goes through dask+PyTables, ``saveData`` is a stub.
* ``zarr_xarray`` -- both halves need the ``zarr`` package.
* ``Class``   -- ``saveData`` is a stub by design (a class is not a file).

Environment, and why some tests skip: ``rasterio`` is installed, so the
``tif`` reader is driven against a real three-by-four GeoTIFF written into
``tmp_path``.  ``osgeo`` (GDAL) is *not* installed, so ``geotiff.getData``
is exercised on its documented "gdal not installed" branch rather than on a
real raster.  ``zarr`` and ``tables`` are *not* installed, so those two
readers are guarded with ``pytest.importorskip`` per the convention in this
suite.

The four ``saveData`` stubs get one contract test each.  ``.coveragerc``
excludes ``raise NotImplementedError`` lines from coverage, so these buy no
coverage -- they are here only so that a future implementation cannot land
silently without a decision about the contract.

B262 is pinned below: ``zarr_xarray.getData`` can never return.  Because
``zarr`` is missing, that pin is driven with ``xarray.open_zarr`` stubbed
out -- the defect is in the line *after* the load, so the stub only removes
the optional dependency, it does not manufacture the failure.

Deliberately not covered here: the ``string``/``time``/``csv_pandas``/
``JSON_*``/``geopandas``/``parquet``/``image``/``pickle``/``dict``/
``numpy_*`` handlers and the ``datatypes`` dispatch tables -- see
``test_datalayer_datahandler.py``, ``_more.py`` and ``_formats.py``.
"""
import numpy
import pytest

from hera.datalayer.datahandler import (DataHandler_Class, DataHandler_geotiff,
                                        DataHandler_HDF, DataHandler_tif,
                                        DataHandler_zarr_xarray, datatypes)


@pytest.fixture()
def geotiff_file(tmp_path):
    """A real single-band float GeoTIFF, written with rasterio."""
    rasterio = pytest.importorskip("rasterio")
    from rasterio.transform import from_origin

    payload = numpy.arange(12, dtype="float32").reshape(3, 4)
    path = str(tmp_path / "small.tif")
    with rasterio.open(path, "w", driver="GTiff", height=3, width=4, count=1,
                       dtype="float32", crs="EPSG:4326",
                       transform=from_origin(10.0, 20.0, 1.0, 1.0)) as dst:
        dst.write(payload, 1)
    return path, payload


# ---------------------------------------------------------------------------
# The saveData stubs
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestUnimplementedSavers:
    """Four handlers are read-only.  Each has to say so, not corrupt a file."""

    @pytest.mark.parametrize("handler, message", [
        (DataHandler_geotiff, "Not implemented yet"),
        (DataHandler_HDF, "HDF saver not implemented yet"),
        (DataHandler_tif, "tif format is not implemented"),
        (DataHandler_Class, "Saving a Class datatype is not supported"),
    ])
    def test_saving_is_refused_with_an_explanation(self, handler, message, tmp_path):
        with pytest.raises(NotImplementedError, match=message):
            handler.saveData("anything", str(tmp_path / "out"))

    @pytest.mark.parametrize("handler", [DataHandler_geotiff, DataHandler_HDF,
                                         DataHandler_tif, DataHandler_Class])
    def test_no_file_is_created_by_the_refusal(self, handler, tmp_path):
        target = tmp_path / "out"
        with pytest.raises(NotImplementedError):
            handler.saveData("anything", str(target))
        assert not target.exists()

    @pytest.mark.parametrize("handler", [DataHandler_geotiff, DataHandler_HDF,
                                         DataHandler_tif, DataHandler_Class])
    def test_extra_store_parameters_do_not_change_the_refusal(self, handler, tmp_path):
        """All four accept **kwargs, so a caller passing store parameters gets
        the same answer rather than a TypeError."""
        with pytest.raises(NotImplementedError):
            handler.saveData("anything", str(tmp_path / "out"), compression="gzip")


# ---------------------------------------------------------------------------
# geotiff
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGeotiffReader:
    def test_gdal_is_not_installed_in_this_environment(self):
        """States the precondition for the two tests below, so that a future
        environment with GDAL fails loudly here instead of silently changing
        what they mean."""
        with pytest.raises(ImportError):
            import osgeo.gdal  # noqa: F401

    def test_a_missing_gdal_is_reported_as_an_import_error(self, geotiff_file):
        path, _ = geotiff_file
        with pytest.raises(ImportError, match="gdal not installed"):
            DataHandler_geotiff.getData(path)

    def test_the_original_module_not_found_error_is_chained(self, geotiff_file):
        """``raise ... from exc`` -- the cause must survive so the operator can
        see which import failed."""
        path, _ = geotiff_file
        with pytest.raises(ImportError) as info:
            DataHandler_geotiff.getData(path)
        assert isinstance(info.value.__cause__, ModuleNotFoundError)

    def test_the_failure_does_not_depend_on_the_resource(self, tmp_path):
        """The import is attempted before the resource is touched, so a path
        that does not exist gives the same ImportError, not an OSError."""
        with pytest.raises(ImportError, match="gdal not installed"):
            DataHandler_geotiff.getData(str(tmp_path / "no-such-file.tif"))

    def test_the_documented_extra_arguments_are_accepted(self, geotiff_file):
        path, _ = geotiff_file
        with pytest.raises(ImportError):
            DataHandler_geotiff.getData(path, rasterBand=2, desc={"a": 1})


# ---------------------------------------------------------------------------
# tif
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestTifReader:
    def test_it_returns_an_open_rasterio_dataset(self, geotiff_file):
        path, _ = geotiff_file
        dataset = DataHandler_tif.getData(path)
        try:
            assert not dataset.closed
        finally:
            dataset.close()

    def test_the_pixels_come_back_unchanged(self, geotiff_file):
        path, payload = geotiff_file
        dataset = DataHandler_tif.getData(path)
        try:
            assert numpy.array_equal(dataset.read(1), payload)
        finally:
            dataset.close()

    def test_the_raster_shape_matches_the_file_that_was_written(self, geotiff_file):
        path, payload = geotiff_file
        dataset = DataHandler_tif.getData(path)
        try:
            assert (dataset.height, dataset.width) == payload.shape
            assert dataset.count == 1
        finally:
            dataset.close()

    def test_the_coordinate_reference_system_is_preserved(self, geotiff_file):
        """CLAUDE.md treats the CRS as load-bearing metadata, so the reader
        must not drop it."""
        path, _ = geotiff_file
        dataset = DataHandler_tif.getData(path)
        try:
            assert dataset.crs.to_epsg() == 4326
        finally:
            dataset.close()

    def test_the_dataset_remembers_the_resource_it_was_opened_from(self, geotiff_file):
        path, _ = geotiff_file
        dataset = DataHandler_tif.getData(path)
        try:
            assert dataset.name == path
        finally:
            dataset.close()

    def test_a_description_argument_is_accepted_and_ignored(self, geotiff_file):
        path, payload = geotiff_file
        dataset = DataHandler_tif.getData(path, desc={"crs": 2039})
        try:
            assert numpy.array_equal(dataset.read(1), payload)
        finally:
            dataset.close()

    def test_a_missing_file_is_reported_by_rasterio(self, tmp_path):
        with pytest.raises(Exception) as info:
            DataHandler_tif.getData(str(tmp_path / "no-such-file.tif"))
        assert not isinstance(info.value, ImportError)


# ---------------------------------------------------------------------------
# HDF
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestHDFReader:
    """Needs PyTables, which is not installed here -- both tests skip."""

    @pytest.fixture()
    def hdf_file(self, tmp_path):
        pytest.importorskip("tables")
        import pandas

        frame = pandas.DataFrame({"a": [1.0, 2.0, 3.0]},
                                 index=pandas.RangeIndex(3, name="idx"))
        path = str(tmp_path / "store.h5")
        frame.to_hdf(path, key="the_key", format="table")
        return path, frame

    def test_it_reads_the_requested_key_as_a_dask_frame(self, hdf_file):
        import dask.dataframe

        path, frame = hdf_file
        result = DataHandler_HDF.getData({"path": path, "key": "the_key"})
        assert isinstance(result, dask.dataframe.DataFrame)
        assert list(result.compute()["a"]) == list(frame["a"])

    def test_use_pandas_computes_eagerly(self, hdf_file):
        import pandas

        path, frame = hdf_file
        result = DataHandler_HDF.getData({"path": path, "key": "the_key"},
                                         usePandas=True)
        assert isinstance(result, pandas.DataFrame)
        assert list(result["a"]) == list(frame["a"])


# ---------------------------------------------------------------------------
# zarr_xarray
# ---------------------------------------------------------------------------

def _dataset():
    import xarray

    data = xarray.Dataset({"a": ("x", [1.0, 2.0, 3.0])}, coords={"x": [0, 1, 2]})
    data.attrs["origin"] = "unit test"
    return data


@pytest.mark.unit
class TestZarrWriter:
    """Needs the zarr package, which is not installed here -- these skip."""

    def test_it_writes_an_archive_that_xarray_can_reopen(self, tmp_path):
        pytest.importorskip("zarr")
        import xarray

        target = str(tmp_path / "store.zarr")
        assert DataHandler_zarr_xarray.saveData(_dataset(), target) == {}
        assert list(xarray.open_zarr(target)["a"].values) == [1.0, 2.0, 3.0]

    def test_writing_twice_overwrites_rather_than_appends(self, tmp_path):
        """Documented: "Rewrites on the file.  If you want to append, do it
        manually." """
        pytest.importorskip("zarr")
        import xarray

        target = str(tmp_path / "store.zarr")
        DataHandler_zarr_xarray.saveData(_dataset(), target)
        DataHandler_zarr_xarray.saveData(_dataset(), target)
        assert xarray.open_zarr(target)["a"].size == 3


@pytest.mark.unit
class TestZarrReaderIsBroken:
    """B262: see the module docstring."""

    @pytest.fixture()
    def stubbed_open_zarr(self, monkeypatch):
        """Removes the missing optional dependency, nothing else.

        The defect under test is on the line *after* the load, so what
        ``open_zarr`` returns is irrelevant: any dataset reaches the same
        statement.
        """
        import xarray

        loaded = _dataset()
        monkeypatch.setattr(xarray, "open_zarr", lambda resource, **kwargs: loaded)
        return loaded

    @pytest.mark.xfail(
        strict=True,
        reason="B262: DataHandler_zarr_xarray.getData loads the archive into "
               "`df` and then writes `resource.attrs = "
               "JSONToConfiguration(resource.attrs)` -- `resource` is the "
               "path string it was handed, not the dataset, so it raises "
               "AttributeError: 'str' object has no attribute 'attrs' before "
               "returning.  The line should read `df.attrs = "
               "JSONToConfiguration(df.attrs)`, mirroring the deserialisation "
               "that DataHandler_netcdf_xarray.getData performs.  No zarr "
               "resource can be read through the datalayer today. "
               "See the consolidated findings issue.",
    )
    def test_it_should_return_the_loaded_dataset(self, stubbed_open_zarr):
        assert DataHandler_zarr_xarray.getData("some/path.zarr") is stubbed_open_zarr

    def test_it_currently_raises_on_the_path_string(self, stubbed_open_zarr):
        """Characterisation of B262."""
        with pytest.raises(AttributeError, match="'str' object has no attribute 'attrs'"):
            DataHandler_zarr_xarray.getData("some/path.zarr")

    def test_the_failure_is_after_the_load_not_before_it(self, monkeypatch):
        """Characterisation of B262: open_zarr is reached and does return,
        which is what makes this a plain typo rather than a missing
        dependency."""
        import xarray

        calls = []

        def _spy(resource, **kwargs):
            calls.append(resource)
            return _dataset()

        monkeypatch.setattr(xarray, "open_zarr", _spy)
        with pytest.raises(AttributeError):
            DataHandler_zarr_xarray.getData("some/path.zarr")
        assert calls == ["some/path.zarr"]

    def test_a_resource_that_does_carry_attrs_gets_past_the_line(self, stubbed_open_zarr):
        """Characterisation of B262: handing getData a *dataset* instead of
        a path makes it return, which localises the defect to the operand.
        No caller in hera does that -- MetadataFrame.getData always passes
        document.resource."""
        import xarray

        resource = _dataset()
        assert isinstance(DataHandler_zarr_xarray.getData(resource), xarray.Dataset)


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestHandlerLookup:
    @pytest.mark.parametrize("constant, handler", [
        ("GEOTIFF", DataHandler_geotiff),
        ("HDF", DataHandler_HDF),
        ("ZARR_XARRAY", DataHandler_zarr_xarray),
        ("CLASS", DataHandler_Class),
    ])
    def test_each_format_constant_resolves_to_its_handler(self, constant, handler):
        assert datatypes.getHandler(getattr(datatypes, constant)) is handler

    def test_the_tif_handler_has_no_format_constant_of_its_own(self, ):
        """``tif`` is reachable by name but is not in the ``datatypes``
        constants -- GEOTIFF is the advertised raster format, so a document
        has to spell 'tif' literally to get the rasterio reader."""
        assert not any(getattr(datatypes, name) == "tif"
                       for name in dir(datatypes) if name.isupper())
        assert datatypes.getHandler("tif") is DataHandler_tif
