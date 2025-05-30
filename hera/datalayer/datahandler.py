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
from ..utils import loadJSON


version = sys.version_info[0]
if version == 3:
    from json import JSONDecodeError
elif version == 2:
    from simplejson import JSONDecodeError

class datatypes:
    STRING = "string"
    TIME   = "time"
    CSV_PANDAS = "csv_pandas"
    HDF    = "HDF"
    NETCDF_XARRAY = "netcdf_xarray"
    JSON_DICT = "JSON_dict"
    JSON_PANDAS = "JSON_pandas"
    JSON_GEOPANDAS = "JSON_geopandas"
    GEOPANDAS  = "geopandas"
    GEOTIFF    = "geotiff"
    PARQUET    =  "parquet"
    IMAGE = "image"
    PICKLE = "pickle"
    DICT   = "dict"

    type_to_datatype_string = {
            str: datatypes.STRING,
            pd.DataFrame: datatypes.CSV_PANDAS,
            gpd.GeoDataFrame: datatypes.GEOPANDAS,
            xr.Dataset: datatypes.NETCDF_XARRAY,
            dict: datatypes.JSON_DICT,
            list: datatypes.JSON_DICT,  # If applicable
            Image.Image: datatypes.IMAGE,
            bytes: datatypes.PICKLE,
            object: datatypes.PICKLE
        }

    datatype_to_extension = {
        datatypes.STRING: ".txt",
        datatypes.TIME: ".txt",  # or .time if you have a custom convention
        datatypes.CSV_PANDAS: ".csv",
        datatypes.HDF: ".h5",
        datatypes.NETCDF_XARRAY: ".nc",
        datatypes.JSON_DICT: ".json",
        datatypes.JSON_PANDAS: ".json",
        datatypes.JSON_GEOPANDAS: ".geojson",
        datatypes.GEOPANDAS: ".gpkg",  # or .shp
        datatypes.GEOTIFF: ".tif",
        datatypes.PARQUET: ".parquet",
        datatypes.IMAGE: ".png",  # could also be .jpg, .tiff depending on use
        datatypes.PICKLE: ".pkl",
        datatypes.DICT: ".json",  # often serialized as JSON
    }


def guessHandler(datatype):
    """
        Tries to estimate the type of the object and re
    Parameters
    ----------
    obj

    Returns
    -------

    """

    return getHandler(type_to_datatype_string[datatype])

def getHandler(objType):
    return globals()[f"DataHandler_{objType}"]

class DataHandler_geotiff(object):
    """
        Loads a single key from HDF file or files.

        Returns a pandas or a dask dataframe.

        The structure of the resource is a dictionary with the keys:
         -  path: the path to the HDF file (can be a pattern to represent a list of files).
         -  key : a single key.
    """

    @staticmethod
    def saveData(resource,fileName):
        raise NotImplementedError("Not implemented yet")

    @staticmethod
    def getData(resource,rasterBand=1):
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
    def saveData(resource,fileName):
        with open(fileName,"w") as outFile:
            outFile.write(resource)


    @staticmethod
    def getData(resource,desc=None):
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
    def saveData(resource, fileName):
        raise NotImplementedError("time is not implemented for saving")

    @staticmethod
    def getData(resource,desc=None):
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
    def saveData(resource,fileName):
        resource.to_csv(fileName)

    @staticmethod
    def getData(resource,desc=None):
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

        df = pandas.read_csv(resource)

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
    def saveData(resource,fileName):
        raise NotImplementedError("HDF saver not implemented yet")

    @staticmethod
    def getData(resource, usePandas=False,desc=None):
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
    def saveData(resource,fileName):
        resource.to_netcdf(fileName)

    @staticmethod
    def getData(resource,desc=None,**kwargs):
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
        df = xarray.open_mfdataset(resource, combine='by_coords',**kwargs)

        return df


class DataHandler_JSON_dict(object):

    @staticmethod
    def saveData(resource,fileName):
        with open(fileName,"w") as outFile:
            json.dump(resource,outFile)


    @staticmethod
    def getData(resource,desc=None):
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
    def saveData(resource, fileName):
        resource.to_json(fileName)

    @staticmethod
    def getData(resource, usePandas=True,desc=None):
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
            df = pandas.read_json(resource)
        else:
            df = dask.dataframe.read_json(resource)

        return df


class DataHandler_JSON_geopandas(object):


    @staticmethod
    def saveData(resource, fileName):
        resource.to_json(fileName)

    @staticmethod
    def getData(resource,desc=None,**kwargs):
        df = geopandas.GeoDataFrame.from_features(loadJSON(resource)["features"])
        if "crs" in desc:
            df.crs = desc['crs']

        return df


class DataHandler_geopandas(object):

    @staticmethod
    def saveData(resource, fileName):
        resource.to_file(fileName, driver="GPKG")


    @staticmethod
    def getData(resource,desc=None,**kwargs):

        df = geopandas.read_file(resource,**kwargs)
        if "crs" in desc:
            df.crs = desc['crs']

        return df


class DataHandler_parquet(object):

    @staticmethod
    def saveData(resource, fileName):
        resource.to_parquet(fileName)

    @staticmethod
    def getData(resource,desc=None, usePandas=False,**kwargs):
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
            df = dask.dataframe.read_parquet(resource,**kwargs)
            if usePandas:
                df = df.compute()
        except ValueError:
            # dask cannot read parquet with multi index. so we try to load it with pandas.
            df = pandas.read_parquet(resource,**kwargs)

        return df


class DataHandler_image(object):

    @staticmethod
    def saveData(resource, fileName):
        mpimg.imsave(fileName, resource)

    @staticmethod
    def getData(resource,desc=None):
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
    def saveData(resource, fileName):
        with open(fileName, 'wb') as f:
            pickle.dump(resource, f)

    @staticmethod
    def getData(resource,desc=None):
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
        obj = pickle.load(resource)

        return obj

class DataHandler_dict(object):

    @staticmethod
    def saveData(resource,fileName):
        pass

    @staticmethod
    def getData(resource,desc=None):
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
    def saveData(resource,fileName):
        raise NotImplementedError("tif format is not implemented")

    @staticmethod
    def getData(resource,desc=None):
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