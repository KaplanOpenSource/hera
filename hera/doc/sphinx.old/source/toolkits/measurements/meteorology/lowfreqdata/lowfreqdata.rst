Parsing the data and getting pandas DataFrame
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Currently, the tool is able to parse the IMS-JSON format. So, first
download file from https://ims.data.gov.il/ims/7

.. code:: ipython3

    from hera.measurements.meteorology.lowfreqdata import lowfreqDataLayer


.. parsed-literal::

    WARNING:root:FreeCAD not Found, cannot convert to STL


.. parsed-literal::

    You must install python-wrf to use this package 


Loading the data from the file. We use the default IMS-JSON foramat.

.. code:: ipython3

    loadedData = lowfreqDataLayer.getStationDataFromFile(pathToData="ims-results-json.json", 
                                                         fileFormat=lowfreqDataLayer.JSONIMS,
                                                         storeDB =False)

The data that was loaded is returned as a metadata document. Therefore,
we need to use getData()

.. code:: ipython3

    print(loadedData.getData().head())


.. parsed-literal::

                          stn_name  stn_num      BP  DiffR  Grad   NIP  RH  Rain  \
    time_obs                                                                       
    2020-09-01 00:00:00  BET DAGAN       54  1001.4      0 -9999 -9999   0     0   
    2020-09-01 00:10:00  BET DAGAN       54  1001.3      0 -9999 -9999   0     0   
    2020-09-01 00:20:00  BET DAGAN       54  1001.2      0 -9999 -9999   0     0   
    2020-09-01 00:30:00  BET DAGAN       54  1001.2      0 -9999 -9999   0     0   
    2020-09-01 00:40:00  BET DAGAN       54  1001.2      0 -9999 -9999   0     0   
    
                         STDwd    TD  TDmax  TDmin    TG  Time   WD  WDmax   WS  \
    time_obs                                                                      
    2020-09-01 00:00:00    7.9  27.6   27.6   27.6  27.4  2351  212    223  2.0   
    2020-09-01 00:10:00    8.7  27.7   27.8   27.6  27.9     1  208    195  1.9   
    2020-09-01 00:20:00    7.8  27.9   27.9   27.8  27.7    18  200    191  2.1   
    2020-09-01 00:30:00   10.4  27.9   27.9   27.8  27.6    21  199    192  1.9   
    2020-09-01 00:40:00    9.0  27.8   27.8   27.7  27.4    37  202    187  1.8   
    
                         WS1mm  WSmax  Ws10mm  
    time_obs                                   
    2020-09-01 00:00:00    2.2    2.7     2.0  
    2020-09-01 00:10:00    2.8    3.3     2.0  
    2020-09-01 00:20:00    2.6    3.0     2.1  
    2020-09-01 00:30:00    2.6    3.0     2.1  
    2020-09-01 00:40:00    2.5    3.0     1.9  


Loading data to the database.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to parse the data, converty it to parquet and load it to
the database.

This is useful when there is a large dataset. When the data is retrieved
from the database, it is returned as dask.

.. code:: ipython3

    loadedData = lowfreqDataLayer.getStationDataFromFile(pathToData="ims-results-json.json", 
                                                         fileFormat=lowfreqDataLayer.JSONIMS,
                                                         storeDB =True)

We got the file we loaded and it was added to the DB. If the station
already exists there, then the data is appended.

.. code:: ipython3

    print(loadedData.getData().head())


.. parsed-literal::

                          stn_name  stn_num      BP  DiffR  Grad   NIP  RH  Rain  \
    time_obs                                                                       
    2020-09-01 00:00:00  BET DAGAN       54  1001.4      0 -9999 -9999   0     0   
    2020-09-01 00:10:00  BET DAGAN       54  1001.3      0 -9999 -9999   0     0   
    2020-09-01 00:20:00  BET DAGAN       54  1001.2      0 -9999 -9999   0     0   
    2020-09-01 00:30:00  BET DAGAN       54  1001.2      0 -9999 -9999   0     0   
    2020-09-01 00:40:00  BET DAGAN       54  1001.2      0 -9999 -9999   0     0   
    
                         STDwd    TD  TDmax  TDmin    TG  Time   WD  WDmax   WS  \
    time_obs                                                                      
    2020-09-01 00:00:00    7.9  27.6   27.6   27.6  27.4  2351  212    223  2.0   
    2020-09-01 00:10:00    8.7  27.7   27.8   27.6  27.9     1  208    195  1.9   
    2020-09-01 00:20:00    7.8  27.9   27.9   27.8  27.7    18  200    191  2.1   
    2020-09-01 00:30:00   10.4  27.9   27.9   27.8  27.6    21  199    192  1.9   
    2020-09-01 00:40:00    9.0  27.8   27.8   27.7  27.4    37  202    187  1.8   
    
                         WS1mm  WSmax  Ws10mm  
    time_obs                                   
    2020-09-01 00:00:00    2.2    2.7     2.0  
    2020-09-01 00:10:00    2.8    3.3     2.0  
    2020-09-01 00:20:00    2.6    3.0     2.1  
    2020-09-01 00:30:00    2.6    3.0     2.1  
    2020-09-01 00:40:00    2.5    3.0     1.9  


Getting data from DB
--------------------

Getting the data from the DB is very similar to parsing and loading a
file.

Listing all the stations in the DB
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, we will list the stations in the DB:

.. code:: ipython3

    listStations = lowfreqDataLayer.listStations()
    print(listStations)


.. parsed-literal::

    ['BET DAGAN']


Getting data from the database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Get the data of the station from the database.

.. code:: ipython3

    datadb = lowfreqDataLayer.getStationDataFromDB(StationName=listStations[0])
    print(datadb.getData().head())


.. parsed-literal::

                          stn_name  stn_num      BP  DiffR  Grad   NIP  RH  Rain  \
    time_obs                                                                       
    2020-09-01 00:00:00  BET DAGAN       54  1001.4      0 -9999 -9999   0     0   
    2020-09-01 00:10:00  BET DAGAN       54  1001.3      0 -9999 -9999   0     0   
    2020-09-01 00:20:00  BET DAGAN       54  1001.2      0 -9999 -9999   0     0   
    2020-09-01 00:30:00  BET DAGAN       54  1001.2      0 -9999 -9999   0     0   
    2020-09-01 00:40:00  BET DAGAN       54  1001.2      0 -9999 -9999   0     0   
    
                         STDwd    TD  TDmax  TDmin    TG  Time   WD  WDmax   WS  \
    time_obs                                                                      
    2020-09-01 00:00:00    7.9  27.6   27.6   27.6  27.4  2351  212    223  2.0   
    2020-09-01 00:10:00    8.7  27.7   27.8   27.6  27.9     1  208    195  1.9   
    2020-09-01 00:20:00    7.8  27.9   27.9   27.8  27.7    18  200    191  2.1   
    2020-09-01 00:30:00   10.4  27.9   27.9   27.8  27.6    21  199    192  1.9   
    2020-09-01 00:40:00    9.0  27.8   27.8   27.7  27.4    37  202    187  1.8   
    
                         WS1mm  WSmax  Ws10mm  
    time_obs                                   
    2020-09-01 00:00:00    2.2    2.7     2.0  
    2020-09-01 00:10:00    2.8    3.3     2.0  
    2020-09-01 00:20:00    2.6    3.0     2.1  
    2020-09-01 00:30:00    2.6    3.0     2.1  
    2020-09-01 00:40:00    2.5    3.0     1.9  



