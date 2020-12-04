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

    The current file directory is /ibdata/yehudaa/Development/pyhera/hera/doc/source


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
    5  -2.091363  0.044787
    6  -2.044366  0.049358
    7  -1.997369  0.054276
    8  -1.950372  0.059551
    9  -1.903376  0.065196
    10 -1.856379  0.071218
    11 -1.809382  0.077625
    12 -1.762385  0.084421
    13 -1.715388  0.091610
    14 -1.668391  0.099192
    15 -1.621394  0.107164
    16 -1.574397  0.115521
    17 -1.527400  0.124256
    18 -1.480403  0.133356
    19 -1.433406  0.142806
    20 -1.386409  0.152590
    21 -1.339412  0.162683
    22 -1.292415  0.173062
    23 -1.245419  0.183696
    24 -1.198422  0.194554
    25 -1.151425  0.205599
    26 -1.104428  0.216792
    27 -1.057431  0.228089
    28 -1.010434  0.239446
    29 -0.963437  0.250814
    ..       ...       ...
    70  0.963437  0.250814
    71  1.010434  0.239446
    72  1.057431  0.228089
    73  1.104428  0.216792
    74  1.151425  0.205599
    75  1.198422  0.194554
    76  1.245419  0.183696
    77  1.292415  0.173062
    78  1.339412  0.162683
    79  1.386409  0.152590
    80  1.433406  0.142806
    81  1.480403  0.133356
    82  1.527400  0.124256
    83  1.574397  0.115521
    84  1.621394  0.107164
    85  1.668391  0.099192
    86  1.715388  0.091610
    87  1.762385  0.084421
    88  1.809382  0.077625
    89  1.856379  0.071218
    90  1.903376  0.065196
    91  1.950372  0.059551
    92  1.997369  0.054276
    93  2.044366  0.049358
    94  2.091363  0.044787
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
    5  -2.091363  0.089574
    6  -2.044366  0.098716
    7  -1.997369  0.108551
    8  -1.950372  0.119103
    9  -1.903376  0.130392
    10 -1.856379  0.142436
    11 -1.809382  0.155249
    12 -1.762385  0.168842
    13 -1.715388  0.183220
    14 -1.668391  0.198383
    15 -1.621394  0.214327
    16 -1.574397  0.231042
    17 -1.527400  0.248511
    18 -1.480403  0.266711
    19 -1.433406  0.285613
    20 -1.386409  0.305179
    21 -1.339412  0.325366
    22 -1.292415  0.346123
    23 -1.245419  0.367392
    24 -1.198422  0.389108
    25 -1.151425  0.411198
    26 -1.104428  0.433583
    27 -1.057431  0.456178
    28 -1.010434  0.478892
    29 -0.963437  0.501628
    ..       ...       ...
    70  0.963437  0.501628
    71  1.010434  0.478892
    72  1.057431  0.456178
    73  1.104428  0.433583
    74  1.151425  0.411198
    75  1.198422  0.389108
    76  1.245419  0.367392
    77  1.292415  0.346123
    78  1.339412  0.325366
    79  1.386409  0.305179
    80  1.433406  0.285613
    81  1.480403  0.266711
    82  1.527400  0.248511
    83  1.574397  0.231042
    84  1.621394  0.214327
    85  1.668391  0.198383
    86  1.715388  0.183220
    87  1.762385  0.168842
    88  1.809382  0.155249
    89  1.856379  0.142436
    90  1.903376  0.130392
    91  1.950372  0.119103
    92  1.997369  0.108551
    93  2.044366  0.098716
    94  2.091363  0.089574
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
            "$oid": "5f9527f84fde6232dbb04e67"
        },
        "dataFormat": "parquet",
        "desc": {
            "loc": 0.5,
            "scale": 0.5
        },
        "projectName": "ExampleProject",
        "resource": "/ibdata/yehudaa/Development/pyhera/hera/doc/source/dataset3.parquet",
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

