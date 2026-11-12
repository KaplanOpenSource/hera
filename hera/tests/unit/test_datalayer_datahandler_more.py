"""datalayer/datahandler.py: the remaining format handlers that don't need
a heavy geospatial/HDF/zarr dependency -- guessHandler's module-level
wrapper, DataHandler_tif (real GeoTIFF via rasterio, which is installed),
DataHandler_numpy_dict_array (.npz round trip), and DataHandler_Class
(dynamic class loading, both the sys.path and resource-fallback paths).

DataHandler_geotiff/_HDF/_zarr_xarray need osgeo/tables/zarr respectively,
none of which are installed even under the pinned requirements.txt in this
environment -- left uncovered, same as earlier batches.
"""
import sys

import numpy
import pytest

from hera.datalayer import datahandler


@pytest.mark.unit
class TestGuessHandlerModuleFunction:
    def test_it_delegates_to_datatypes_guess_handler(self):
        import pandas

        handler = datahandler.guessHandler(pandas.DataFrame({"a": [1]}))
        assert handler is datahandler.DataHandler_parquet


@pytest.mark.unit
class TestDataHandlerTif:
    @pytest.fixture()
    def tif_path(self, tmp_path):
        rasterio = pytest.importorskip("rasterio")
        from rasterio.transform import from_origin

        path = tmp_path / "test.tif"
        data = numpy.arange(16, dtype="uint8").reshape(4, 4)
        transform = from_origin(0, 4, 1, 1)
        with rasterio.open(
            str(path), "w", driver="GTiff", height=4, width=4, count=1,
            dtype="uint8", crs="EPSG:4326", transform=transform,
        ) as dst:
            dst.write(data, 1)
        return path

    def test_get_data_opens_a_real_rasterio_dataset(self, tif_path):
        obj = datahandler.DataHandler_tif.getData(str(tif_path))
        try:
            assert numpy.array_equal(obj.read(1), numpy.arange(16, dtype="uint8").reshape(4, 4))
        finally:
            obj.close()

    def test_save_data_is_not_implemented(self):
        with pytest.raises(NotImplementedError, match="tif format"):
            datahandler.DataHandler_tif.saveData(None, "x")


@pytest.mark.unit
class TestDataHandlerNumpyDictArray:
    def test_round_trips_a_dict_of_arrays(self, tmp_path):
        target = tmp_path / "arrays.npz"
        resource = {"a": numpy.array([1, 2, 3]), "b": numpy.array([[1.0, 2.0]])}
        meta = datahandler.DataHandler_numpy_dict_array.saveData(resource, str(target))
        assert meta == {}
        restored = datahandler.DataHandler_numpy_dict_array.getData(str(target))
        assert numpy.array_equal(restored["a"], resource["a"])
        assert numpy.array_equal(restored["b"], resource["b"])


@pytest.mark.unit
class TestDataHandlerClass:
    @pytest.fixture()
    def package_dir(self, tmp_path):
        pkg = tmp_path / "myclasspkg"
        pkg.mkdir()
        (pkg / "__init__.py").write_text("")
        (pkg / "mymod.py").write_text(
            "class MyClass:\n"
            "    def __init__(self, a=1, b=2):\n"
            "        self.a = a\n"
            "        self.b = b\n"
        )
        return tmp_path

    def test_it_instantiates_a_class_resolvable_on_sys_path(self, package_dir, monkeypatch):
        monkeypatch.syspath_prepend(str(package_dir))
        obj = datahandler.DataHandler_Class.getData(
            resource=None,
            desc={"classpath": "myclasspkg.mymod.MyClass", "parameters": {"a": 10}},
            b=99,
        )
        assert (obj.a, obj.b) == (10, 99)
        sys.modules.pop("myclasspkg", None)
        sys.modules.pop("myclasspkg.mymod", None)

    def test_desc_parameters_override_kwargs_on_conflict(self, package_dir, monkeypatch):
        monkeypatch.syspath_prepend(str(package_dir))
        obj = datahandler.DataHandler_Class.getData(
            resource=None,
            desc={"classpath": "myclasspkg.mymod.MyClass", "parameters": {"a": 10}},
            a=1,
        )
        assert obj.a == 10
        sys.modules.pop("myclasspkg", None)
        sys.modules.pop("myclasspkg.mymod", None)

    def test_instantiate_false_returns_the_class_not_an_instance(self, package_dir, monkeypatch):
        monkeypatch.syspath_prepend(str(package_dir))
        cls = datahandler.DataHandler_Class.getData(
            resource=None,
            desc={"classpath": "myclasspkg.mymod.MyClass", "instantiate": False},
        )
        assert isinstance(cls, type)
        assert cls.__name__ == "MyClass"
        sys.modules.pop("myclasspkg", None)
        sys.modules.pop("myclasspkg.mymod", None)

    def test_falls_back_to_loading_from_the_resource_directory(self, package_dir):
        """When the module isn't importable via sys.path, it's loaded
        directly from the resource directory instead."""
        obj = datahandler.DataHandler_Class.getData(
            resource=str(package_dir), desc={"classpath": "myclasspkg.mymod.MyClass"},
        )
        assert (obj.a, obj.b) == (1, 2)
        sys.modules.pop("myclasspkg", None)
        sys.modules.pop("myclasspkg.mymod", None)

    def test_a_missing_classpath_raises(self):
        with pytest.raises(ValueError, match="classpath"):
            datahandler.DataHandler_Class.getData(resource=None, desc={})

    def test_save_data_is_not_implemented(self):
        with pytest.raises(NotImplementedError, match="Saving a Class datatype"):
            datahandler.DataHandler_Class.saveData(None, "x")
