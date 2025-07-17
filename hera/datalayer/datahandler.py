import numpy
import pandas
import dask.dataframe
import xarray
import json
import geopandas
try:
    from osgeo import gdal
except ImportError:
    print("gdal not installed, no support for shapefiles")

import matplotlib.image as mpimg
import sys
import pickle
import io
import rasterio
from hera.utils import loadJSON
import importlib

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
