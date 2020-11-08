import pandas
import dask.dataframe
import xarray
import json
import geopandas
import matplotlib.image as mpimg
import sys
import pickle
version = sys.version_info[0]
if version == 3:
    from json import JSONDecodeError
elif version == 2:
    from simplejson import JSONDecodeError

class datatypes:
    STRING = "string"
    TIME   = "time"
    CSV_PANDAS = "csv_pandas"
    HDF    = "hdf"
    NETCDF_XARRAY = "netcdf_xarray"
    JSON_DICT = "JSON_dict"
    JSON_PANDAS = "JSON_pandas"
    GEOPANDAS  = "geopandas"
    PARQUET    =  "parquet"
    IMAGE = "image"
    PICKLE = "pickle"


def getHandler(type):
    return globals()['DataHandler_%s' % type]


class DataHandler_string(object):
    """
        The resource is a string.
    """

    @staticmethod
    def getData(resource):
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
    def getData(resource):
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
    def getData(resource):
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
    def getData(resource, usePandas=False):
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
    def getData(resource):
        """
        Loads netcdf file into xarray.

        Parameters
        ----------
        resource : str
            Path to the netcdf file.

        Returns
        -------
        xarray
        """
        df = xarray.open_mfdataset(resource, combine='by_coords')

        return df


class DataHandler_JSON_dict(object):

    @staticmethod
    def getData(resource):
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
        try:
            df = json.loads(resource)
        except JSONDecodeError:
            with open(resource, 'r') as myFile:
                df = json.load(myFile)

        return df


class DataHandler_JSON_pandas(object):

    @staticmethod
    def getData(resource, usePandas=True):
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


class DataHandler_geopandas(object):
    @staticmethod
    def getData(resource):
        df = geopandas.read_file(resource)

        return df


class DataHandler_parquet(object):

    @staticmethod
    def getData(resource, usePandas=False):
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

        df = dask.dataframe.read_parquet(resource)
        if usePandas:
            df = df.compute()

        return df


class DataHandler_image(object):

    @staticmethod
    def getData(resource):
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
    def getData(resource):
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