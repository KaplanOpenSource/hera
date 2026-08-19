import importlib
import json

import cloudpickle

from hera.utils.jsonutils import ConfigurationToJSON, JSONToConfiguration
from hera.utils.lazy import _LazyModule

numpy  = _LazyModule("numpy")
pandas = _LazyModule("pandas")


class datatypes:
    """
    Registry of supported data format constants and dispatch logic for data handlers.

    Each constant (e.g. ``STRING``, ``PARQUET``, ``HDF``) identifies a data format.
    Use ``getHandler(formatName)`` to retrieve the corresponding ``DataHandler_*`` class,
    or ``getDataFormatName(obj)`` to auto-detect the format from a Python object.
    """
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
        """
        Auto-detect the data format and return the appropriate handler class.

        Parameters
        ----------
        obj_or_class : object or type
            The data object or class to detect the format for.

        Returns
        -------
        DataHandler class
            The handler class for the detected format.
        """
        dataTypeName = datatypes.getDataFormatName(obj_or_class)

        return datatypes.getHandler(objectType=dataTypeName)

    @staticmethod
    def getHandler(objectType):
        """
        Return the DataHandler class for the given data format name.

        Parameters
        ----------
        objectType : str
            A data format name (e.g. ``datatypes.PARQUET``).

        Returns
        -------
        DataHandler class

        Raises
        ------
        ValueError
            If no handler exists for the given type.
        """
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
        """Save geotiff data to file (not implemented)."""
        raise NotImplementedError("Not implemented yet")

    @staticmethod
    def getData(resource, rasterBand=1, desc=None, **kwargs):
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
        try:
            from osgeo import gdal
        except ImportError as exc:
            raise ImportError("gdal not installed, no support for shapefiles") from exc
        ds = gdal.Open(resource)
        return ds


class DataHandler_string(object):
    """
        The resource is a string.
    """

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        """Save a string resource to a text file."""
        with open(fileName, "w") as outFile:
            outFile.write(resource)
        return dict()

    @staticmethod
    def getData(resource, desc=None, **kwargs):
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
    def getData(resource, desc=None, **kwargs):
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
        """Save a pandas DataFrame to a CSV file."""
        resource.to_csv(fileName,**kwargs)
        return dict()

    @staticmethod
    def getData(resource, desc=None, **kwargs):
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
        """Save data to HDF format (not implemented)."""
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
        import dask.dataframe
        df = dask.dataframe.read_hdf(resource['path'], resource['key'], sorted_index=True)

        if usePandas:
            df = df.compute()

        return df


class DataHandler_netcdf_xarray(object):
    """Handler for xarray datasets stored as NetCDF files."""
    SERIALIZED_ATTRS_KEY = "_serialized_attrs"

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        """Save an xarray dataset to NetCDF format."""
        flattened_attrs = ConfigurationToJSON(resource.attrs)
        resource.attrs = {DataHandler_netcdf_xarray.SERIALIZED_ATTRS_KEY:json.dumps(flattened_attrs)}
        resource.to_netcdf(fileName,**kwargs)
        return dict()

    @staticmethod
    def getData(resource, desc=None, **kwargs):
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
        import xarray
        ds = xarray.open_mfdataset(resource, combine='by_coords', **kwargs)
        if DataHandler_netcdf_xarray.SERIALIZED_ATTRS_KEY in ds.attrs:
            ds.attrs = json.loads(ds.attrs[DataHandler_netcdf_xarray.SERIALIZED_ATTRS_KEY])
        ds.attrs = JSONToConfiguration(ds.attrs)

        return ds


class DataHandler_zarr_xarray(object):
    """Handler for xarray datasets stored as Zarr archives."""

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
        resource.attrs = ConfigurationToJSON(resource.attrs)
        resource.to_zarr(fileName, mode="w",**kwargs)
        return dict()

    @staticmethod
    def getData(resource, desc=None, **kwargs):
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
        import xarray
        df = xarray.open_zarr(resource, **kwargs)
        resource.attrs = JSONToConfiguration(resource.attrs)
        return df


class DataHandler_JSON_dict(object):
    """Handler for Python dicts stored as JSON files."""

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        """Save a dict to a JSON file."""
        with open(fileName, "w") as outFile:
            json.dump(resource, outFile,**kwargs)
        return dict()

    @staticmethod
    def getData(resource, desc=None, **kwargs):
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
        from hera.utils.jsonutils import loadJSON
        df = loadJSON(resource)
        return df


class DataHandler_JSON_pandas(object):
    """Handler for pandas DataFrames/Series stored as JSON files."""

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        """Save a pandas DataFrame or Series to JSON format."""
        resource.to_json(fileName,**kwargs)

        if isinstance(resource, pandas.Series):
            ret = dict(pandasSeries=True,usePandas=True)
        elif isinstance(resource,pandas.DataFrame):
            ret = dict(usePandas=True)
        else:
            ret = dict()

        return ret

    @staticmethod
    def getData(resource, usePandas=True, desc=None, **kwargs):
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
            import dask.dataframe
            df = dask.dataframe.read_json(resource)

        return df


class DataHandler_JSON_geopandas(object):
    """Handler for GeoDataFrames stored as GeoJSON files."""

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        """Save a GeoDataFrame to GeoJSON format."""
        resource.to_json(fileName,**kwargs)
        return dict(crs = resource.crs )

    @staticmethod
    def getData(resource, desc=None, **kwargs):
        """Load a GeoDataFrame from a GeoJSON file."""
        import geopandas

        from hera.utils.jsonutils import loadJSON
        desc = desc or {}
        df = geopandas.GeoDataFrame.from_features(loadJSON(resource)["features"])
        if "crs" in desc:
            df.crs = desc['crs']

        return df


class DataHandler_geopandas(object):
    """Handler for GeoDataFrames stored as GeoPackage files."""

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        """Save a GeoDataFrame to GeoPackage (GPKG) format."""
        resource.to_file(fileName, driver="GPKG",**kwargs)
        return dict(crs=resource.crs)

    @staticmethod
    def getData(resource, desc=None, **kwargs):
        """Load a GeoDataFrame from a geospatial file."""
        import geopandas
        desc = desc or {}
        df = geopandas.read_file(resource, **kwargs)
        if "crs" in desc:
            df.crs = desc['crs']

        return df


class DataHandler_parquet(object):
    """Handler for pandas/dask DataFrames stored as Parquet files."""

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        """Save a pandas or dask DataFrame to Parquet format."""
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
    def getData(resource, desc=None, usePandas=False, **kwargs):
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
        import dask.dataframe
        try:
            df = dask.dataframe.read_parquet(resource, **kwargs)
            if usePandas:
                df = df.compute()
        except ValueError:
            df = pandas.read_parquet(resource, **kwargs)

        return df


class DataHandler_image(object):
    """Handler for image data stored via matplotlib."""

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        """Save an image array to a file using ``matplotlib.image.imsave``."""
        import matplotlib.image as mpimg
        mpimg.imsave(fileName, resource,**kwargs)
        return dict()

    @staticmethod
    def getData(resource, desc=None, **kwargs):
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
        import matplotlib.image as mpimg
        img = mpimg.imread(resource)

        return img


class DataHandler_pickle(object):
    """Handler for arbitrary Python objects stored as pickle files."""

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        """Save an object to a pickle file."""
        with open(fileName, 'wb') as f:
            cloudpickle.dump(resource, f,**kwargs)
        return dict()

    @staticmethod
    def getData(resource, desc=None, **kwargs):
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
            obj = cloudpickle.load(f)

        return obj


class DataHandler_dict(object):
    """Handler for dict data stored directly in the document resource field."""

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
    """Handler for GeoTIFF raster files (read-only via rasterio)."""

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        """Not implemented; raises ``NotImplementedError``."""
        raise NotImplementedError("tif format is not implemented")

    @staticmethod
    def getData(resource, desc=None, **kwargs):
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
        try:
            import rasterio
        except ImportError as exc:
            raise ImportError("rasterio not installed, no support for image data types.") from exc
        obj = rasterio.open(resource)

        return obj


class DataHandler_numpy_array:
    """Handler for numpy arrays stored as ``.npy`` files."""

    @staticmethod
    def saveData(resource, fileName,**kwargs):
        """Save a numpy array to ``.npy`` format."""
        numpy.save(fileName, resource,**kwargs)
        return dict()

    @staticmethod
    def getData(resource, desc=None, **kwargs):
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
        """Save a dict of numpy arrays to ``.npz`` format."""
        numpy.savez(fileName, **resource,**kwargs)
        return dict()

    @staticmethod
    def getData(resource, desc=None, **kwargs):
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
        Filesystem directory that contains the code for the class (will be added to sys.path).
        If empty/None, nothing is added to sys.path.
        If resource points directly to a Python package directory (contains __init__.py),
        the parent directory is also added so that `import top_pkg.module` resolves.

    desc:
        - 'classpath' (str): fully qualified import path of the class,
          e.g. 'mypkg.mymodule.MyClass' (REQUIRED).
        - 'params' or 'parameters' (dict, optional): keyword-args for the class constructor.
        - 'instantiate' (bool, optional; default True):
            If True, return an instance (cls(**...)).
            If False, return the class object (cls) itself.

    Merge rule (Option B):
        When both desc.parameters and **kwargs provide the same key,
        desc.parameters should take precedence (override kwargs).
    """

    @staticmethod
    def saveData(resource, fileName, **kwargs):
        """Not implemented; class datatypes cannot be saved to disk."""
        # Storing a "Class" datatype as a file is not supported by this handler.
        raise NotImplementedError("Saving a Class datatype is not supported")

    @staticmethod
    def getData(resource, desc=None, **kwargs):
        """Import and optionally instantiate a class from ``desc['classpath']``.

        Uses importlib.util.spec_from_file_location for resource-based loading
        so that sys.path is never mutated by DB-supplied paths [1.6, 3.3].
        """
        import importlib
        import importlib.util
        import os

        # 1) Resolve metadata
        desc = desc or {}
        classpath = desc.get("classpath") or kwargs.get("classpath")
        if not classpath:
            raise ValueError('For dataFormat=Class you must provide desc["classpath"]')

        params = desc.get("parameters") or desc.get("params") or {}
        instantiate = desc.get("instantiate", True)

        # 2) Import module and get class by name
        module_name, _, class_name = classpath.rpartition(".")
        if not module_name or not class_name:
            raise ValueError(
                f"Invalid classpath '{classpath}'. Expected something like 'pkg.mod.Class'."
            )

        # Try loading via existing sys.path first (safe path)
        try:
            module = importlib.import_module(module_name)
        except ModuleNotFoundError:
            # Fall back to loading from resource directory without mutating sys.path
            if not resource:
                raise ImportError(
                    f"Cannot import module '{module_name}' and no resource path provided."
                )
            abs_resource = os.path.abspath(resource)
            parts = module_name.split(".")
            # resource may point to the top-level package dir or to its parent.
            # When its basename matches the first package component, strip that
            # component so we look inside the package dir rather than for a
            # nested sub-directory with the same name.
            if parts[0] == os.path.basename(abs_resource):
                inner_parts = parts[1:]
            else:
                inner_parts = parts
            module_file = os.path.join(abs_resource, os.sep.join(inner_parts) + ".py")
            if not os.path.isfile(module_file):
                raise ImportError(
                    f"Cannot find module '{module_name}' in sys.path or in {abs_resource!r}"
                )
            spec = importlib.util.spec_from_file_location(module_name, module_file)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
        except Exception as e:
            raise ImportError(
                f"Cannot import module '{module_name}' for classpath '{classpath}': {e}"
            )

        try:
            cls = getattr(module, class_name)
        except AttributeError:
            raise ImportError(f"Module '{module_name}' has no attribute '{class_name}'")

        # 3) Merge constructor kwargs so that desc.parameters override duplicates
        call_kwargs = dict(kwargs)
        call_kwargs.update(params)

        # 4) Return an instance or the class object
        return cls(**call_kwargs) if instantiate else cls
