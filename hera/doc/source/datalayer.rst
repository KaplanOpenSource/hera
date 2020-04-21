datalayer
=========

The data the that is produced in by the measurements or the simulation might be large
and therefore we usually store the data on the disk
and keep only the metadata in the database. Each record in the database (called a 'document') represents
a single piece of data (pandas, dask, xarray or any other data) associated with a project.
The document also include information on the format the the data is stored and other fields to describe
the data.

.. image:: Hera-DB.png

Each user has its own database that stores documents, but it is also possible to access
a different database.

The datalayer also includes a command line interface to load, remove and move data around.



10-minute tutorial
------------------

This tutorial demonstrate how to store and retrieve data from the database.

Adding Data
***********

first, import the requested module:

.. code-block:: python

    from hera import datalayer

Each data is classified into one of the following categories.

- Measurements - Any acquisition of data from the 'real world'. Satellites, meteorological measurments and dispersion measurements and ect.
- Simulations  - Any output of a model. (OpenFOAM, WRF, LSM and ect).
- Analysis     - Any data that is created during work and analyis and needed to be cached to accelerate the computations.

Since each category can hold different types of data, each data document
holds its type in  a 'type' property
when saving a new document in the database, you must provide 'projectName', 'desc', 'dataFormat' and 'resource'

for example, a Measurements data from type Meteorology and parquet data format should be added like this:

.. code-block:: python

    datalayer.Measurements.addDocument(projectName='myProject',
                                       resource='path-to-parquet-files',
                                       dataFormat='parquet',
                                       type='Meteorology'
                                       desc={'property1': 'value1',
                                             'property2': 'value2'
                                             }
                                       )

Retrieving data
***************

When reading from the database, you have to provide the 'projectName' and create a query using the function 'getDocuments'.
First, you get all of the documents in the project that property1 equals 'value1'

.. code-block:: python

    docList = datalayer.Measurements.getDocuments(projectName='projectName',property1='value1')

Now, you can use the function getData() to retrieve data presented by the document.
For example, we took the first document from our docList and got its data.

.. code-block:: python

    doc = docList[0]
    data = doc.getData()


The getData method loads the data from the database accoding to the dataFormat.

Usage
-----

.. toctree::
    :maxdepth: 2

    datalayer/cli
    datalayer/usage
    datalayer/datahandlers
    datalayer/api
