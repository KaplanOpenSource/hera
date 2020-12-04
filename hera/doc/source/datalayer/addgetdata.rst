Adding data
===========

First lets define some mock-up data. Lets create 3 data from 3
distributions with different mean and std.

.. code:: ipython3

    import json
    import os
    import pandas
    import numpy
    from scipy.stats import norm
    import  matplotlib.pyplot as plt 
    from hera import datalayer
    
    x = numpy.linspace(norm.ppf(0.01), norm.ppf(0.99), 100)
    
    dataset1 = pandas.DataFrame(dict(x=x,y=norm.pdf(x,loc=0,scale=1)))
    dataset2 = pandas.DataFrame(dict(x=x,y=norm.pdf(x,loc=0,scale=0.5)))
    dataset3 = pandas.DataFrame(dict(x=x,y=norm.pdf(x,loc=0.5,scale=0.5)))
    
    print(dataset1.head())


.. parsed-literal::

    WARNING:root:FreeCAD not Found, cannot convert to STL


.. parsed-literal::

    You must install python-wrf to use this package 
              x         y
    0 -2.326348  0.026652
    1 -2.279351  0.029698
    2 -2.232354  0.033020
    3 -2.185357  0.036632
    4 -2.138360  0.040550


Now that we have data, we can save it. We would like to keep the
connection between the data and the parameters that generated it. so
that \* dataset1 will be described by loc=0 and scale = 1 \* dataset2
will be described by loc=0 and scale = 0.5 \* dataset3 will be described
by loc=0.5 and scale = 0.5

Therefore, we will save the loc and scale as metadata.

First, we create the file name that we we want to save to.

.. code:: ipython3

    workingdir = os.path.abspath(os.path.dirname(os.getcwd()))
    print(f"The current file directory is {workingdir}")


.. parsed-literal::

    The current file directory is /home/yehuda/Development/Hera/hera/doc/source


It is very important to save the **absolute path** (using
os.path.abspath) and not the relative path because that way the loading
process is independent of the location it is perfomed in.

.. code:: ipython3

    dataset1File = os.path.join(workingdir,"dataset1.parquet")
    dataset2File = os.path.join(workingdir,"dataset2.parquet")
    dataset3File = os.path.join(workingdir,"dataset3.parquet")

Now we can save the dataset. We choose parquet for convinience, but it
can be any other format.

.. code:: ipython3

    dataset1.to_parquet(dataset1File,engine='fastparquet',compression='GZIP')
    dataset2.to_parquet(dataset2File,engine='fastparquet',compression='GZIP')
    dataset3.to_parquet(dataset3File,engine='fastparquet',compression='GZIP')

When we save the data to the database we need to specify the project
name that the record is related to. we use

.. code:: ipython3

    projectName = "ExampleProject"

in our example.

Next, we add the documents to the database.

.. code:: ipython3

    datalayer.Measurements.addDocument(projectName=projectName,
                                       type="Distribution",
                                       dataFormat=datalayer.datatypes.PARQUET,
                                       resource=dataset1File,
                                       desc=dict(loc=0,scale=1))
    
    datalayer.Measurements.addDocument(projectName=projectName,
                                       type="Distribution",
                                       dataFormat=datalayer.datatypes.PARQUET,
                                       resource=dataset2File,
                                       desc=dict(loc=0,scale=0.5))
    
    datalayer.Measurements.addDocument(projectName=projectName,
                                       type="Distribution",
                                       dataFormat=datalayer.datatypes.PARQUET,
                                       resource=dataset3File,
                                       desc=dict(loc=0.5,scale=0.5))




.. parsed-literal::

    <Measurements: Measurements object>



The type of the document was chosen arbitrarily and can be any string.
This string helps in future queries of the data. It can also be an empty
string.

The desc property includes the metadata in a JSON format. It can be any
valid JSON.

Each data is classified into one of the following categories.

-  Measurements - Any acquisition of data from the ‘real world’.
   Satellites, meteorological measurments and dispersion measurements
   and etc.
-  Simulations - Any output of a model. (OpenFOAM, WRF, LSM and etc).
-  Cache - Any data that is created during analyis and needed to be
   cached to accelerate the computations.

Getting the data
================

Getting one record back
-----------------------

Now we will query the database for all the records in which loc=0 and
scale=1.

.. code:: ipython3

    List1 = datalayer.Measurements.getDocuments(projectName=projectName,loc=0,scale=1)
    
    print(f"The number of documents obtained from the query {len(List1)} ")
    item0 = List1[0]



.. parsed-literal::

    The number of documents obtained from the query 1 


Note that for consistency the query always returns a list.

The description of the record that matched the query is

.. code:: ipython3

    print("The description of dataset 1")
    print(json.dumps(item0.desc, indent=4, sort_keys=True))


.. parsed-literal::

    The description of dataset 1
    {
        "loc": 0,
        "scale": 1
    }


Now, we will extract the data.

.. code:: ipython3

    dataset1FromDB = item0.getData().compute()
    
    print(dataset1FromDB)


.. parsed-literal::

               x         y
    0  -2.326348  0.026652
    1  -2.279351  0.029698
    2  -2.232354  0.033020
    3  -2.185357  0.036632
    4  -2.138360  0.040550
    ..       ...       ...
    95  2.138360  0.040550
    96  2.185357  0.036632
    97  2.232354  0.033020
    98  2.279351  0.029698
    99  2.326348  0.026652
    
    [100 rows x 2 columns]


Getting multiple records back
-----------------------------

If the query is specified in a more general way. Lets get all the
records in which loc=0

.. code:: ipython3

    List2 = datalayer.Measurements.getDocuments(projectName=projectName,loc=0)
    
    print(f"The number of documents obtained from the query {len(List2)} ")


.. parsed-literal::

    The number of documents obtained from the query 2 


Updating the data.
==================

The hera system holds the name of the file on the disk and loads the
data from it. Therefore, if the datafile on the disk is overwitten, then
the data of the record is changed

Lets multiply dataset1 by factor 2. The file name is saved in the
resource attribute.

.. code:: ipython3

    dataset1['y'] *=2
    dataset1FileName = item0.resource 
    dataset1.to_parquet(dataset1FileName,engine='fastparquet',compression='GZIP',append=False)

.. code:: ipython3

    item0 = datalayer.Measurements.getDocuments(projectName=projectName,loc=0,scale=1)[0]
    dataset1FromDB = item0.getData().compute()
    print(dataset1FromDB)


.. parsed-literal::

               x         y
    0  -2.326348  0.053304
    1  -2.279351  0.059397
    2  -2.232354  0.066040
    3  -2.185357  0.073264
    4  -2.138360  0.081099
    ..       ...       ...
    95  2.138360  0.081099
    96  2.185357  0.073264
    97  2.232354  0.066040
    98  2.279351  0.059397
    99  2.326348  0.053304
    
    [100 rows x 2 columns]


Updating the metadata.
======================

Lets assume we want to add another property to the first record. To so
we wiill no update item0

.. code:: ipython3

    item0.desc['new_attribute'] = "some data"
    item0.save()




.. parsed-literal::

    <Measurements: Measurements object>



.. code:: ipython3

    item0_fromdb = datalayer.Measurements.getDocuments(projectName=projectName,loc=0,scale=1)[0]
    print(json.dumps(item0_fromdb.desc, indent=4, sort_keys=True))


.. parsed-literal::

    {
        "loc": 0,
        "new_attribute": "some data",
        "scale": 1
    }


Using Project
=============

Using the Project class simplifies the access to the different documents
of the project.

Define the project with

.. code:: ipython3

    from hera.datalayer import Project 
    
    p = Project(projectName=projectName)
    
    results = p.getMeasurementsDocuments(loc=0)
    [x.desc for x in results]




.. parsed-literal::

    [{'loc': 0, 'scale': 1, 'new_attribute': 'some data'},
     {'loc': 0, 'scale': 0.5}]



Deleting the metadata entry.
============================

We delete the metadata records similarly to the way we add them

The following will delete one record

.. code:: ipython3

    docdict = datalayer.Measurements.deleteDocuments(projectName=projectName,loc=0.5,scale=0.5)
    print("The deleted document")
    print(json.dumps(docdict[0], indent=4, sort_keys=True))


.. parsed-literal::

    The deleted document
    {
        "_cls": "Metadata.Measurements",
        "_id": {
            "$oid": "5f9475a6dbe374fa63669a8c"
        },
        "dataFormat": "parquet",
        "desc": {
            "loc": 0.5,
            "scale": 0.5
        },
        "projectName": "ExampleProject",
        "resource": "/home/yehuda/Development/Hera/hera/doc/source/dataset3.parquet",
        "type": "Distribution"
    }


Now we can erase the file from the disk. It is saved in the resource
property

.. code:: ipython3

    import shutil 
    
    if os.path.isfile(docdict[0]['resource']):
          os.remove(docdict[0]['resource'])
    else: 
        shutil.rmtree(docdict[0]['resource'])

Now, we can delete several documents

.. code:: ipython3

    docdictList = datalayer.Measurements.deleteDocuments(projectName=projectName,loc=0)
    
    for doc in docdictList:
        if os.path.isfile(doc['resource']):
            os.remove(doc['resource'])
        else: 
            shutil.rmtree(doc['resource'])


Using the project allows getting documents

