import numpy
import pandas
import dask.dataframe
import xarray
import json
import geopandas
from osgeo import gdal
import matplotlib.image as mpimg
import sys
import pickle
import io
import rasterio
from hera.utils import loadJSON
import importlib
import os


version = sys.version_info[0]
if version == 3:
    from json import JSONDecodeError
elif version == 2:
    from simplejson import JSONDecodeError




class datatypes:
    STRING = "string"
    TIME = "time"
    CSV_PANDAS = "csv_pandas"
    HDF = "HDF"
    NETCDF_XARRAY = "netcdf_xarray"
    ZARR_XARRAY = "zarr_xarray"
    JSON_DICT = "JSON_dict"
    JSON_PANDAS = "JSON_pandas"
    JSON_GEOPANDAS = "JSON_geopandas"
    GEOPANDAS = "geopandas"
    GEOTIFF = "geotiff"
    PARQUET = "parquet"
    IMAGE = "image"
    PICKLE = "pickle"
    DICT = "dict"
    NUMPY_ARRAY = "numpy_array"
    NUMPY_DICT_ARRAY = "numpy_dict_array"  # A dict of numpy arrays, no automatic detection.
    CLASS = "Class"

    @staticmethod
    def get_obj_or_instance_fullName(obj):
        """
        Returns the fully qualified name of a class or instance, including its module.

        Examples:
            >>> get_full_name(SomeClass)
            'package.module.SomeClass'

            >>> get_full_name(SomeClass())
            'package.module.SomeClass'
        """
        # If it's a class
        if isinstance(obj, type):
            cls = obj
        else:
            cls = obj.__class__

        module = cls.__module__
        qualname = cls.__qualname__

        if module == "builtins":
            return qualname  # No need to show 'builtins' for int, str, etc.
        return f"{module}.{qualname}"

    typeDatatypeMap = {
        "str": dict(typeName=STRING, ext="txt"),
        "pandas.core.frame.DataFrame": dict(typeName=PARQUET, ext="parquet"),
        'pandas.core.series.Series': dict(typeName=JSON_PANDAS, ext="json"),
        "dask_expr._collection.DataFrame": dict(typeName=PARQUET, ext="parquet"),
        'geopandas.geodataframe.GeoDataFrame': dict(typeName=GEOPANDAS, ext="gpkg"),
        'xarray.core.dataarray.DataArray': dict(typeName=ZARR_XARRAY, ext="zarr"),
        "dict": dict(typeName=PICKLE, ext="pckle"),
        "list": dict(typeName=PICKLE, ext="pckle"),
        "bytes": dict(typeName=PICKLE, ext="pckle"),
        "object": dict(typeName=PICKLE, ext="pckle"),
        "numpy.ndarray": dict(typeName=NUMPY_ARRAY, ext="npy")
    }

    @staticmethod
    def getDataFormatName(obj_or_class):
        """
            Tries to find the datatype name in hera for the object.
            if cannot found, use general object.

        Parameters
        ----------
        obj_or_class : object or type.

        Returns
        -------
            A dict with
                - typeName : the string that identifies the datahandler.
                -ext : the extension of the file name.
        """
        objTypeName = datatypes.get_obj_or_instance_fullName(obj_or_class)


        dataItemName = datatypes.typeDatatypeMap["object"] if objTypeName not in datatypes.typeDatatypeMap else \
        datatypes.typeDatatypeMap[objTypeName]

        return dataItemName["typeName"]

    @staticmethod
    def getDataFormatExtension(obj_or_class):
        """
            Tries to find the datatype name in hera for the object.
            if cannot found, use general object.

        Parameters
        ----------
        obj_or_class : object or type.

        Returns
        -------
            A dict with
                - typeName : the string that identifies the datahandler.
                -ext : the extension of the file name.
        """
        objTypeName = datatypes.get_obj_or_instance_fullName(obj_or_class)


        dataItemName = datatypes.typeDatatypeMap["object"] if objTypeName not in datatypes.typeDatatypeMap else \
        datatypes.typeDatatypeMap[objTypeName]

        return dataItemName["ext"]

    @staticmethod
    def guessHandler(obj_or_class):

        dataTypeName = datatypes.getDataFormatName(obj_or_class)

        return datatypes.getHandler(objectType=dataTypeName)

    @staticmethod
    def getHandler(objectType):
        dataHandlerModule = importlib.import_module("hera.datalayer.datahandler")

        handlerName = f"DataHandler_{objectType}"

        if not hasattr(dataHandlerModule, handlerName):
            raise ValueError(f"The data handler for the type {objectType} is not known")

        return getattr(dataHandlerModule, handlerName)

dataFormatClass = datatypes

def guessHandler(obj_or_class):
    """
        Tries to estimate the type of the object and re
    Parameters
    ----------
    obj

    Returns
    -------

    """
    return datatypes.guessHandler(obj_or_class)


def getHandler(objType):
    """
        Returns the handler.
    Parameters
    ----------
    objType

    Returns
    -------

    """
    return datatypes.getHandler(objType)


class DataHandler_geotiff(object):
    """
        Loads a single key from HDF file or files.

        Returns a pandas or a dask dataframe.

        The structure of the resource is a dictionary with the keys:
         -  path: the path to the HDF file (can be a pattern to represent a list of files).
         -  key : a single key.
    """

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        raise NotImplementedError("Not implemented yet")

    @staticmethod
    def getData(resource, rasterBand=1):
        """
        Reads a geotiff

        Parameters
        ----------
        resource : dict
            A dictionary with path to the HDF file in the 'path' key, and HDF key in the 'key' key.


        Returns
        -------
            The gdal at image coordinates.
        """
        ds = gdal.Open(resource)
        return ds


class DataHandler_string(object):
    """
        The resource is a string.
    """

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        with open(fileName, "w") as outFile:
            outFile.write(resource)
        return dict()

    @staticmethod
    def getData(resource, desc={},**kwargs):
        """
        The data in the record is a string.

        Parameters
        ----------
        resource : str
            String

        Returns
        -------
        str
        """
        return resource


class DataHandler_time(object):
    """
        The resource is a timestamp.
    """

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        """
            No need to save as it is stored in the resource.
        Parameters
        ----------
        resource
        fileName
        kwargs

        Returns
        -------

        """
        return dict()

    @staticmethod
    def getData(resource, desc={},**kwargs):
        """
        The data in the record is a timestamp.

        Parameters
        ----------
        resource : timestamp
            Time

        Returns
        -------
        pandas.Timestamp
        """
        return pandas.Timestamp(resource)


class DataHandler_csv_pandas(object):
    """
        Loads a csv file into pandas dataframe.

        Returns pandas dataframe.
    """

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        resource.to_csv(fileName,**kwargs)
        return dict()

    @staticmethod
    def getData(resource, desc={},**kwargs):
        """
        Loads a csv file into pandas dataframe.

        Parameters
        ----------
        resource : str
            Path to a csv file

        Returns
        -------
        panda.dataframe
        """

        df = pandas.read_csv(resource,**kwargs)

        return df


class DataHandler_HDF(object):
    """
        Loads a single key from HDF file or files.

        Returns a pandas or a dask dataframe.

        The structure of the resource is a dictionary with the keys:
         -  path: the path to the HDF file (can be a pattern to represent a list of files).
         -  key : a single key.
    """

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        raise NotImplementedError("HDF saver not implemented yet")

    @staticmethod
    def getData(resource, usePandas=False, desc=None):
        """
        Loads a key from a HDF file or files.

        Parameters
        ----------
        resource : dict
            A dictionary with path to the HDF file in the 'path' key, and HDF key in the 'key' key.

        usePandas : bool, optional, default True
            if False use dask if True use pandas.

        Returns
        -------
        dask.DataFrame or pandas.DataFrame
        """
        df = dask.dataframe.read_hdf(resource['path'], resource['key'], sorted_index=True)

        if usePandas:
            df = df.compute()

        return df


class DataHandler_netcdf_xarray(object):

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        resource.to_netcdf(fileName,**kwargs)
        return dict()

    @staticmethod
    def getData(resource, desc={}, **kwargs):
        """
        Loads netcdf file into xarray using the open_mfdataset.

        Parameters
        ----------
        resource : str
            Path to the netcdf file.

        kwargs:
            parameters to the xarray.open_mfdataset function

        Returns
        -------
        xarray
        """
        df = xarray.open_mfdataset(resource, combine='by_coords', **kwargs)

        return df


class DataHandler_zarr_xarray(object):

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        """
            Write the zarr. Rewrites on the file.
            If you want to append, do it manuyally.
        Parameters
        ----------
        resource
        fileName

        Returns
        -------

        """
        resource.to_zarr(fileName, mode="w",**kwargs)
        return dict()

    @staticmethod
    def getData(resource, desc={}, **kwargs):
        """
        Loads netcdf file into xarray using the open_mfdataset.

        Parameters
        ----------
        resource : str
            Path to the netcdf file.

        kwargs:
            parameters to the xarray.open_mfdataset function

        Returns
        -------
        xarray
        """
        df = xarray.open_zarr(resource, **kwargs)
        return df


class DataHandler_JSON_dict(object):

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        with open(fileName, "w") as outFile:
            json.dump(resource, outFile,**kwargs)
        return dict()

    @staticmethod
    def getData(resource, desc={},**kwargs):
        """
        Loads JSON to dict

        Parameters
        ----------
        resource : str
            The data in a JSON format.

        Returns
        -------
        dict
        """
        df = loadJSON(resource)
        return df


class DataHandler_JSON_pandas(object):

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        resource.to_json(fileName,**kwargs)

        if isinstance(resource, pandas.Series):
            ret = dict(pandasSeries=True,usePandas=True)
        elif isinstance(resource,pandas.DataFrame):
            ret = dict(usePandas=True)
        else:
            ret = dict()

        return ret

    @staticmethod
    def getData(resource, usePandas=True, desc={},**kwargs):
        """
        Loads JSON to pandas/dask

        Parameters
        ----------
        resource : str
            The data in a JSON Format

        usePandas : bool, optional, default True
            if False use dask if True use pandas.

        Returns
        -------
        pandas.DataFrame or dask.DataFrame
        """
        if usePandas:
            isSeries = kwargs.get("pandasSeries",False)
            if isSeries:
                readParams = dict(typ="series")
            else:
                readParams = dict()
            df = pandas.read_json(resource,**readParams)
        else:
            df = dask.dataframe.read_json(resource)

        return df


class DataHandler_JSON_geopandas(object):

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        resource.to_json(fileName,**kwargs)
        return dict(crs = resource.crs )

    @staticmethod
    def getData(resource, desc={}, **kwargs):
        df = geopandas.GeoDataFrame.from_features(loadJSON(resource)["features"])
        if "crs" in desc:
            df.crs = desc['crs']

        return df


class DataHandler_geopandas(object):

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        resource.to_file(fileName, driver="GPKG",**kwargs)
        return dict(crs=resource.crs)

    @staticmethod
    def getData(resource, desc={}, **kwargs):
        df = geopandas.read_file(resource, **kwargs)
        if "crs" in desc:
            df.crs = desc['crs']

        return df


class DataHandler_parquet(object):

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        if isinstance(resource, pandas.DataFrame):
            # pandas. write as a single file.
            # if any of the columns is integer it breaks the dask
            resource.to_parquet(fileName,**kwargs)
            ret = dict(usePandas=True)
        else:
            # dask - automatically split it
            kwargs.setdefault('partition_size',"100MB")
            resource.to_parquet(fileName,**kwargs)
            ret = dict(usePandas=False)
        return ret

    @staticmethod
    def getData(resource, desc={}, usePandas=False, **kwargs):
        """
        Loads a parquet file to dask/pandas.

        Parameters
        ----------
        resource : str
            The directory of the parquet file.

        usePandas : bool, optional, default False
            if False use dask if True use pandas.

        Returns
        -------
        dask.Dataframe or pandas.DataFrame
        """
        try:
            df = dask.dataframe.read_parquet(resource, **kwargs)
            if usePandas:
                df = df.compute()
        except ValueError:
            # dask cannot read parquet with multi index. so we try to load it with pandas.
            df = pandas.read_parquet(resource, **kwargs)

        return df


class DataHandler_image(object):

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        mpimg.imsave(fileName, resource,**kwargs)
        return dict()

    @staticmethod
    def getData(resource, desc={},**kwargs):
        """
        Loads an image using the resource.

        Parameters
        ----------
        resource : str
            The path of the image.

        Returns
        -------
        img
        """
        img = mpimg.imread(resource)

        return img


class DataHandler_pickle(object):

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        with open(fileName, 'wb') as f:
            pickle.dump(resource, f,**kwargs)
        return dict()

    @staticmethod
    def getData(resource, desc={},**kwargs):
        """
        Loads an pickled object using the resource.

        Parameters
        ----------
        resource : str
            The path to the pickled object

        Returns
        -------
        img
        """
        with open(resource, 'rb') as f:
            obj = pickle.load(f)

        return obj


class DataHandler_dict(object):

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        """
            No need to save as it is stored in the resource.
        Parameters
        ----------
        resource
        fileName
        kwargs

        Returns
        -------

        """
        return dict()

    @staticmethod
    def getData(resource, desc=None):
        """
        The resource is a dict.

        Parameters
        ----------
        resource : dict
            The resrouce

        Returns
        -------
        dict
        """
        return dict(resource)


class DataHandler_tif(object):

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        raise NotImplementedError("tif format is not implemented")

    @staticmethod
    def getData(resource, desc={},**kwargs):
        """
        Loads an pickled object using the resource.

        Parameters
        ----------
        resource : str
            The path to the pickled object

        Returns
        -------
        img
        """
        obj = rasterio.open(resource)

        return obj


class DataHandler_numpy_array:

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        numpy.save(fileName, resource,**kwargs)
        return dict()

    @staticmethod
    def getData(resource, desc={},**kwargs):
        """
        Loads a numpy array

        Parameters
        ----------
        resource : str
            The path to the pickled object

        Returns
        -------
        img
        """
        obj = numpy.load(resource)

        return obj


class DataHandler_numpy_dict_array:
    """
        A dict of numpy arrays.
    """

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        numpy.savez(fileName, **resource,**kwargs)
        return dict()

    @staticmethod
    def getData(resource, desc={},**kwargs):
        """
        Loads a numpy array

        Parameters
        ----------
        resource : str
            The path to the pickled object

        Returns
        -------
        dict of numpy arrays.
        """
        obj = numpy.load(resource)

        return obj


class DataHandler_Class(object):
    """
    DataType: Class

    resource:
        Filesystem directory used to resolve imports:
        - If 'resource' points to the *package directory itself* (contains __init__.py),
          we add ONLY its parent directory to sys.path (so 'import mypkg.mymod' works).
        - If 'resource' points to a *parent directory* that contains the package,
          we add that parent directory to sys.path.
        - If resource is empty/None, we do not modify sys.path.

    desc:
        - 'classpath' (str, REQUIRED): fully qualified class path, e.g. 'mypkg.mymod.MyClass'.
        - 'params' or 'parameters' (dict, optional): kwargs for the class constructor.
        - 'instantiate' (bool, optional; default True):
            If True, return an instance (cls(**...)).
            If False, return the class object (cls).

    Merge rule:
        When both desc.parameters and **kwargs provide the same key,
        desc.parameters take precedence (they override **kwargs).
    """

    @staticmethod
    def saveData(resource, fileName, **kwargs):
        # This handler is for dynamically loading Python classes; saving is not supported.
        raise NotImplementedError("Saving a Class datatype is not supported")

    @staticmethod
    def getData(resource, desc=None, **kwargs):
        """
        Load a Python class dynamically.

        Behavior for 'resource':
        - If 'resource' is a *package directory* (contains __init__.py), we must add
          the *parent directory* to sys.path so that 'import top_pkg.module' resolves.
        - If 'resource' is a directory that *contains* the top-level package, we add
          that directory itself to sys.path.

        Extra robustness:
        - Move the chosen path to the *front* of sys.path (even if it already exists).
        - Purge relevant entries from sys.modules before import (avoid stale cache).
        - Invalidate import caches between attempts.
        """
        import os
        import sys
        import importlib

        # ---------- 1) Resolve classpath ----------
        desc = desc or {}
        classpath = desc.get("classpath") or kwargs.get("classpath")
        if not classpath:
            raise ValueError('For dataFormat=Class you must provide desc["classpath"]')

        module_name, _, class_name = classpath.rpartition(".")
        if not module_name or not class_name:
            raise ValueError(f"Invalid classpath '{classpath}'. Expected 'pkg.mod.Class'.")
        top_pkg = module_name.split(".", 1)[0]

        # ---------- 2) sys.path manipulation helpers ----------
        def move_to_front(path: str):
            """Ensure 'path' is at sys.path[0] (remove previous occurrence if any)."""
            try:
                while path in sys.path:
                    sys.path.remove(path)
            except ValueError:
                pass
            sys.path.insert(0, path)

        def ensure_parent_for_package_dir(pkg_dir: str):
            """resource == package dir: put its *parent* at sys.path[0]."""
            parent = os.path.dirname(pkg_dir)
            move_to_front(parent)

        def ensure_parent_dir_contains_pkg(parent_dir: str):
            """resource == parent dir that contains the package."""
            move_to_front(parent_dir)

        def purge_modules():
            """Drop relevant modules from sys.modules to avoid stale imports."""
            for mod in (top_pkg, module_name):
                if mod in sys.modules:
                    del sys.modules[mod]

        # ---------- 3) Decide which path to put on sys.path ----------
        if resource:
            abs_path = os.path.abspath(resource)
            if os.path.isdir(abs_path):
                pkg_init = os.path.join(abs_path, "__init__.py")
                if os.path.isfile(pkg_init):
                    # Case A: resource is the *package directory*
                    ensure_parent_for_package_dir(abs_path)
                else:
                    # Case B: resource is a parent dir that *may* contain the package
                    cand_pkg = os.path.join(abs_path, top_pkg)
                    if os.path.isdir(cand_pkg) and os.path.isfile(os.path.join(cand_pkg, "__init__.py")):
                        ensure_parent_dir_contains_pkg(abs_path)
                    else:
                        # Fallback: still add the given dir
                        ensure_parent_dir_contains_pkg(abs_path)

        importlib.invalidate_caches()
        purge_modules()  # clear caches BEFORE first import attempt

        # ---------- 4) Import module & get class (with one retry and debug) ----------
        def import_and_get():
            module = importlib.import_module(module_name)
            try:
                return getattr(module, class_name)
            except AttributeError:
                raise ImportError(f"Module '{module_name}' has no attribute '{class_name}'")

        try:
            cls = import_and_get()
        except Exception as e1:
            # Retry once with a stronger path setup: ensure both parent and pkg dir exist at the front
            debug_notes = []
            try:
                if resource:
                    abs_path = os.path.abspath(resource)
                    if os.path.isdir(abs_path):
                        pkg_init = os.path.join(abs_path, "__init__.py")
                        if os.path.isfile(pkg_init):
                            # Make sure BOTH parent and pkg dir are at the very front (parent first)
                            parent = os.path.dirname(abs_path)
                            move_to_front(parent)
                            move_to_front(abs_path)  # harmless; parent is still first
                        else:
                            cand_pkg = os.path.join(abs_path, top_pkg)
                            if os.path.isdir(cand_pkg) and os.path.isfile(os.path.join(cand_pkg, "__init__.py")):
                                move_to_front(abs_path)
                                move_to_front(cand_pkg)
                importlib.invalidate_caches()
                purge_modules()
                cls = import_and_get()
            except Exception as e2:
                # Build rich diagnostics
                try:
                    head_sys_path = list(sys.path[:8])
                except Exception:
                    head_sys_path = []
                try:
                    parent_ls = []
                    pkg_ls = []
                    if resource:
                        ap = os.path.abspath(resource)
                        if os.path.isdir(ap):
                            parent = os.path.dirname(ap)
                            parent_ls = sorted(os.listdir(parent)) if os.path.isdir(parent) else []
                            pkg_ls = sorted(os.listdir(ap)) if os.path.isdir(ap) else []
                    debug_notes.append(f"sys.path(head)={head_sys_path}")
                    debug_notes.append(f"top_pkg={top_pkg} | module_name={module_name}")
                    debug_notes.append(f"resource={resource}")
                    debug_notes.append(f"parent_dir_ls={parent_ls[:20]}")
                    debug_notes.append(f"pkg_dir_ls={pkg_ls[:20]}")
                except Exception:
                    pass
                msg = (
                        f"Cannot import module '{module_name}' for classpath '{classpath}': {e2}\n"
                        + "\n".join(debug_notes)
                )
                raise ImportError(msg) from e2

        # ---------- 5) Build kwargs (desc.parameters override **kwargs) ----------
        params = desc.get("parameters") or desc.get("params") or {}
        instantiate = desc.get("instantiate", True)
        call_kwargs = dict(kwargs)
        call_kwargs.update(params)

        # ---------- 6) Return instance or class ----------
        return cls(**call_kwargs) if instantiate else cls
