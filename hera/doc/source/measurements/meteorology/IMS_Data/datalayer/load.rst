
****************************
Load data
****************************

The data loading fundtion:
==========================

First, we need to import meteo from 'Hera':

.. code-block:: bash

    from hera import meteo

LoadData:
---------
The LoadData function organizes the data upload into the database with the appropriate format and structure.


Since this is particular data type, the following parameters are set by default for files coming from the IMS database
(of course they can be overwrite by the user):

- type='meteorological'
- DataSource='IMS'
- station_column='stn_name'
- time_coloumn='time_obs'

When calling this function, the following parameters must be set:

- **newdata_path:** The path to the folder with the new raw data
- **outputpath:** The path to the folder where the new arranged data will be saved
- **Projectname:** A project name string associated to th new loaded data. When the data is general or not explicitly linked to a particular project,we recommend using the general project 'IMS_Data'


When the data is associated with a metadata file, it can be loaded by passing path to the file with **'metadatafile'** variable.

Here is an example of how to use the LoadData function:

.. code-block:: bash

    meteo.IMS_dataloader.LoadData(newdata_path = '/raid/users/davidg/Project/2020/IMS_RAW/json/multiStations',
                                  outputpath = '/raid/users/davidg/Project/2020/IMS_RAW/parquet/',
                                  Projectname = 'IMS_Data',
                                  metadatafile = '/raid/users/davidg/Project/2020/IMS_RAW/meta.txt'
                                  )

The function performs the following procedures:

1. Converts the new data from json file to dask format
2. Adds the new data to the old data (if exists) by station name (each station saved to separate file)
3. Saves the unified data to a parquet files in the request 'outputpath' location
4. Saves new document to the metadata with the metadata added to the description (for new stations)





