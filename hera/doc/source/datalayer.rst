Datalayer
*********

The data the that is produced in by the measurements or the simulation might be large
and therefore we usually store the data on the disk
and keep only the metadata in the database. Each record in the database (called a 'document') represents
a single piece of data (pandas, dask, xarray or any other data) associated with a **project**.

The meta data includes the fields:

 - The project name that contains the data.
 - The data (or a reference to it).
 - The format of the data on the disk.
 - fields that describe the data.


.. image:: Hera-DB.png


Each user has its own database that stores documents, but it is also possible to access
a different databases.

The datalayer also includes a command line interface to load, remove and move data around.


Setup
=====

The library uses a connection to the databse that is defined in the path
~/.pyhera/config.json. The structure of the config.json is

..  code-block:: javascript

    {
    <connection name 1> : {
        "dbIP": "DB IP",
        "dbName": "...",
        "password": "..." ,
        "username": "..."
    },
    <connection name 1> : {
        "dbIP": "DB IP",
        "dbName": "...",
        "password": "..." ,
        "username": "..."
    }
}

The default connection is the **default with the linux user name**.
It is possible to add other database connections. See below on how
to access other databases.

10-minute tutorial
===================

This tutorial demonstrate how to store and retrieve data from the database
with the default connection (the username of the linux system).

All pieces of data are related to a project.

Listing all projects
---------------------

List all the project by calling the geProjectList function.

.. code-block:: python

    projectList = hera.datalayer.getProjectList()

.. include:: datalayer/addgetdata.rst

Advanced Querying data
======================

Querying data relies on the mongoDB query language to
query the desc. Briefly, in the mongoDB language translates the JSON path to a
the list of keys seperated by '__'.

for example if the desc field of the document is

.. code-block:: javascript

    "a" : {
        "b" : {
            "c" : 1,
            "d" : [1,2,3],
            "e" : "A"
        },
        "b1" : 4
    }

The the path of the "c" key-path is 'a__b__c'.
So querying the documents where 'c' field is 1 is done as follows:

.. code-block:: python

    docList = datalayer.Measurements.getDocuments(projectName='projectName',a__b__c=1)

It is possible to add operators to query all the documents that fulfil a certain criteria.
For example add '__lt' to find all the documents that are less than a value.

.. code-block:: python

    docList = datalayer.Measurements.getDocuments(projectName='projectName',a__b__c__lt=1)
    docList = datalayer.Measurements.getDocuments(projectName='projectName',a__b__d__in=1)

The last line retrieves all the documents that the field 'd' includes the item 1 in it.





API
=====


.. toctree::
    :maxdepth: 2

    datalayer/api

Usage
=====

.. toctree::
    :maxdepth: 2

    datalayer/cli
    datalayer/usage
    datalayer/datahandlers
    datalayer/architecture