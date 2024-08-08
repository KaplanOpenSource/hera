*************
Data handlers
*************

Each record in the mongodb represents a single item of data. In order to
simplify the loading procedure, we keep track of the data format (which is determined by the user when it is created).

The following data handlers are implemented:

- string          :py:class:`hera.datalayer.datahandler.DataHandler_string`
- time            :py:class:`hera.datalayer.datahandler.DataHandler_time`
- HDF             :py:class:`hera.datalayer.datahandler.HDF`
- netcdf_xarray   :py:class:`hera.datalayer.datahandler.netcdf_xarray`
- JSON_dict       :py:class:`hera.datalayer.datahandler.JSON_dict`
- JSON_pandas     :py:class:`hera.datalayer.datahandler.JSON_pandas`
- geopandas       :py:class:`hera.datalayer.datahandler.geopandas`
- parquet         :py:class:`hera.datalayer.datahandler.parquet`
- image           :py:class:`hera.datalayer.datahandler.image`

string
------
.. autoclass:: hera.datalayer.datahandler.DataHandler_string
    :members:
    :undoc-members:

time
----
.. autoclass:: hera.datalayer.datahandler.DataHandler_time
    :members:
    :undoc-members:

HDF
---
.. autoclass:: hera.datalayer.datahandler.DataHandler_HDF
    :members:
    :undoc-members:

netcdf_xarray
-------------
.. autoclass:: hera.datalayer.datahandler.DataHandler_netcdf_xarray
    :members:
    :undoc-members:

JSON_dict
---------
.. autoclass:: hera.datalayer.datahandler.DataHandler_JSON_dict
    :members:
    :undoc-members:

JSON_pandas
-----------
.. autoclass:: hera.datalayer.datahandler.DataHandler_JSON_pandas
    :members:
    :undoc-members:

geopandas
---------
.. autoclass:: hera.datalayer.datahandler.DataHandler_geopandas
    :members:
    :undoc-members:

parquet
-------
.. autoclass:: hera.datalayer.datahandler.DataHandler_parquet
    :members:
    :undoc-members:

image
-----
.. autoclass:: hera.datalayer.datahandler.DataHandler_image
    :members:
    :undoc-members:

