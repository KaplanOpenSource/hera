
****************************
Read data
****************************



The data reading functions:
===========================

First, we need to import meteo from 'Hera'

.. code-block:: bash

    from hera import meteo

getDocFromDB:
-------------
The getDocFromDB function organizes the data reading from the database by the appropriate format and structure.
The function works so that it returns a list of records (documents) that match the requested query

Since this is particular data type, the following parameters are set by default for files coming from the IMS database
(of course they can be overwrite by the user):

- type='meteorological'
- DataSource='IMS'

When calling this function, the following parameters must be set:

- **Projectname:** A project name string associated to the data.

The query should be as accurate as possible so that the returned records are only those you want to work with:
It is recommended to use additional parameters appear in the description during the query request. For example, a StationName parameter.

Here is an example of how to use the getDocFromDB function:

.. code-block:: bash

    docList = meteo.IMS_datalayer.getDocFromDB('IMS_Data',StationName='AVNE ETAN')


getDocFromFile:
---------------
The function allows the user to load a document from a file on the disk and not from the metadata.
The function returns a list of 'metadata-like' documents and allows the user continue working in the same way
regardless of the data source.

When calling this function, the following parameters must be set:

- **path:** A path to the folder with the raw data

Here is an example of how to use the getDocFromFile function:

.. code-block:: bash

    docList=meteo.IMS_datalayer.getDocFromFile('/raid/users/davidg/Project/2020/IMS_RAW/json/BET_DAGAN/')

When calling the function, you can add description properties after the path. the build-in properties
('projectName', 'Resource', ect...) will be added to the document as in a metadata document,and all other properties
will be listed to the 'desc'.



