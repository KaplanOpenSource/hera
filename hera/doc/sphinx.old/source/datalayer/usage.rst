Adding data
-----------

This notebook shows how to add a document to database with the
datalayer.

First, we will create a synthetic dataframe that we want to save as a
database record:

.. code:: ipython3

    import pandas
    from hera import datalayer
    import json 

We now create a mock dataset:

.. code:: ipython3

    df = pandas.DataFrame(dict(date = ["2016-11-10", "2016-11-10", "2016-11-11", "2016-11-11","2016-11-11","2016-11-11","2016-11-11", "2016-11-11" ],
                               time = ["22:00:00", "23:00:00", "00:00:00", "01:00:00", "02:00:00", "03:00:00", "04:00:00", "05:00:00"],
                               value = [90, 91, 80, 87, 84,94, 91, 94]))
    df['date_time'] = pandas.to_datetime(df['date'] + ' ' + df['time'])
    df=df.set_index('date_time')
    print (df)


.. parsed-literal::

                               date      time  value
    date_time                                       
    2016-11-10 22:00:00  2016-11-10  22:00:00     90
    2016-11-10 23:00:00  2016-11-10  23:00:00     91
    2016-11-11 00:00:00  2016-11-11  00:00:00     80
    2016-11-11 01:00:00  2016-11-11  01:00:00     87
    2016-11-11 02:00:00  2016-11-11  02:00:00     84
    2016-11-11 03:00:00  2016-11-11  03:00:00     94
    2016-11-11 04:00:00  2016-11-11  04:00:00     91
    2016-11-11 05:00:00  2016-11-11  05:00:00     94


Several parameters must be given for any document, and must be defined
in order to add new data. These parameters are the ones given in the
next example. In addition, one may add any other parameters to the
document.

In order to add the data, we will determine all its properties:

The projectName must be a string.

.. code:: ipython3

    projectName = "addDataExample"

The document type is a string that is determined by the user. It is used
to simplify the querying of documents of a similar type.

.. code:: ipython3

    documentType = "ExampleData"

The description of the metadata is a dict. The dict can be any JSON
format, at any desired depth (but please be reasonable…).

.. code:: ipython3

    desc = dict(description_A="A", description_B="B")

The dataformat specifies what is the format of the data that the
document describes. In this example, we use JSON_PANDAS format which
means that the ‘resource’ attributes contains the JSON text of the
pandas.

.. code:: ipython3

    dataFormat = datalayer.datatypes.JSON_PANDAS 

In this example, we create the data using the ‘to_json()’ function.
However, the resource will typically include a parquet file.

.. code:: ipython3

    resource = df.to_json()

Now we add the document using hera.datalayer

.. code:: ipython3

    new_doc=datalayer.Measurements.addDocument(projectName=projectName, type=documentType, dataFormat=dataFormat, resource=resource, desc=desc)

We can see the document content by printing it as dict

.. code:: ipython3

    print(new_doc.asDict())


.. parsed-literal::

    {'_cls': 'Metadata.Measurements', 'projectName': 'addDataExample', 'desc': {'description_A': 'A', 'description_B': 'B'}, 'type': 'ExampleData', 'resource': '{"date":{"1478815200000":"2016-11-10","1478818800000":"2016-11-10","1478822400000":"2016-11-11","1478826000000":"2016-11-11","1478829600000":"2016-11-11","1478833200000":"2016-11-11","1478836800000":"2016-11-11","1478840400000":"2016-11-11"},"time":{"1478815200000":"22:00:00","1478818800000":"23:00:00","1478822400000":"00:00:00","1478826000000":"01:00:00","1478829600000":"02:00:00","1478833200000":"03:00:00","1478836800000":"04:00:00","1478840400000":"05:00:00"},"value":{"1478815200000":90,"1478818800000":91,"1478822400000":80,"1478826000000":87,"1478829600000":84,"1478833200000":94,"1478836800000":91,"1478840400000":94}}', 'dataFormat': 'JSON_pandas'}


Notice that the desc dictionary may not contain a key named “type”. The
allowed data formats are detailed in the hera.datalayer.datatypes
(partial list):

-  STRING : Any string.
-  TIME : any date/time object
-  HDF : a dask or pandas in hdf file format.
-  NETCDF_XARRAY : an xarray netcdf.
-  JSON_DICT : JSON as python dict
-  JSON_PANDAS : JSON as pandas.DataFrame
-  GEOPANDAS : a GIS-file format. returns as geopandas.GISDataFrame
-  PARQUET : dask or pandas in parquet format.
-  IMAGE : any Image data format. Preferably PNG.

They indicate how to read the data, and therefore must correspond to the
type of data located in the resource.

The added document can be loaded as presented in the “Getting data”
notebook.

Getting data
------------

Getting the data can achieved by using the getDocuments procedure of the
collection. Getting the data allows the user to `query the
database <#another_cell>`__ using the `mongo engine query
language <https://docs.mongoengine.org/guide/querying.html>`__.

For simplicity, lets retrieve the documents for which ‘description_A’
equals ‘A’.

.. code:: ipython3

    docList = datalayer.Measurements.getDocuments(projectName=projectName,description_A = "A")

The query returns a list.

The result obtained from the query is:

.. code:: ipython3

    print(len(docList))


.. parsed-literal::

    1


In order to see the details of the document, lets print it

.. code:: ipython3

    print(docList[0])


.. parsed-literal::

    Measurements object


.. code:: ipython3

    docList[0].resource




.. parsed-literal::

    '{"date":{"1478815200000":"2016-11-10","1478818800000":"2016-11-10","1478822400000":"2016-11-11","1478826000000":"2016-11-11","1478829600000":"2016-11-11","1478833200000":"2016-11-11","1478836800000":"2016-11-11","1478840400000":"2016-11-11"},"time":{"1478815200000":"22:00:00","1478818800000":"23:00:00","1478822400000":"00:00:00","1478826000000":"01:00:00","1478829600000":"02:00:00","1478833200000":"03:00:00","1478836800000":"04:00:00","1478840400000":"05:00:00"},"value":{"1478815200000":90,"1478818800000":91,"1478822400000":80,"1478826000000":87,"1478829600000":84,"1478833200000":94,"1478836800000":91,"1478840400000":94}}'



You can now read the data from the ‘doc’ and perform another query (for
example, on a date range):

.. code:: ipython3

    data=docList[0].getData()

The data we got is:

.. code:: ipython3

    print(data)


.. parsed-literal::

                              date      time  value
    2016-11-10 22:00:00 2016-11-10  22:00:00     90
    2016-11-10 23:00:00 2016-11-10  23:00:00     91
    2016-11-11 00:00:00 2016-11-11  00:00:00     80
    2016-11-11 01:00:00 2016-11-11  01:00:00     87
    2016-11-11 02:00:00 2016-11-11  02:00:00     84
    2016-11-11 03:00:00 2016-11-11  03:00:00     94
    2016-11-11 04:00:00 2016-11-11  04:00:00     91
    2016-11-11 05:00:00 2016-11-11  05:00:00     94


## Advanced Querying data

Querying the database uses the `mongo engine query
language <https://docs.mongoengine.org/guide/querying.html>`__. Briefly,
in the mongo enngine query language the JSON path translates to a the
list of keys seperated by ’__’.

for example if the desc field of the document is

.. code:: ipython3

    print(json.dumps({"a" : {"b" : {"c" : 1,"d" : [1,2,3],"e" : "A"},"b1" : 4}},indent=4))


.. parsed-literal::

    {
        "a": {
            "b": {
                "c": 1,
                "d": [
                    1,
                    2,
                    3
                ],
                "e": "A"
            },
            "b1": 4
        }
    }


The the path of the “c” key-path is ’a__b__c’. So querying the documents
where ‘c’ field is 1 is done as follows:

.. code:: ipython3

    tmp = datalayer.Measurements.getDocuments(projectName='projectName',a__b__c=1)

It is possible to add operators to query all the documents that fulfil a
certain criteria. For example add ’__lt’ to find all the documents that
are less than a value.

.. code:: ipython3

    tmp = datalayer.Measurements.getDocuments(projectName='projectName',a__b__c__lt=1)


To retrieve all the documents that the field ‘d’ includes the item 1 in
it.

.. code:: ipython3

    tmp = datalayer.Measurements.getDocuments(projectName='projectName',a__b__d__in=1)

Update data description
-----------------------

If the metadata changes, it is possible to update its value. For
example, the data before the update is:

.. code:: ipython3

    print('The resource is: %s' %docList[0].resource)
    print('The description is: %s' %docList[0].desc)


.. parsed-literal::

    The resource is: {"date":{"1478815200000":"2016-11-10","1478818800000":"2016-11-10","1478822400000":"2016-11-11","1478826000000":"2016-11-11","1478829600000":"2016-11-11","1478833200000":"2016-11-11","1478836800000":"2016-11-11","1478840400000":"2016-11-11"},"time":{"1478815200000":"22:00:00","1478818800000":"23:00:00","1478822400000":"00:00:00","1478826000000":"01:00:00","1478829600000":"02:00:00","1478833200000":"03:00:00","1478836800000":"04:00:00","1478840400000":"05:00:00"},"value":{"1478815200000":90,"1478818800000":91,"1478822400000":80,"1478826000000":87,"1478829600000":84,"1478833200000":94,"1478836800000":91,"1478840400000":94}}
    The description is: {'description_A': 'A', 'description_B': 'B'}


.. code:: ipython3

    docobj = docList[0]
    newdata1 = dict(docobj.desc)
    newdata1['description_C'] = "C1"
    resource1 = "resource1"
    
    
    newdata2 = dict(docobj.desc)
    newdata2['description_C'] = "C2"
    resource2 = "resource2"

**Method 1:** set the new attributes in the object and save.

.. code:: ipython3

    docobj.resource = resource1
    docobj.desc = newdata1
    docobj.save()




.. parsed-literal::

    <Measurements: Measurements object>



Now we check that the database was updated.

.. code:: ipython3

    after_update_docList = datalayer.Measurements.getDocuments(projectName=projectName,**desc)
    after_update_docobj = docList[0]
    print('The resource is: %s' %after_update_docobj.resource)
    print('The description is: %s' %after_update_docobj.desc)


.. parsed-literal::

    The resource is: resource1
    The description is: {'description_A': 'A', 'description_B': 'B', 'description_C': 'C1'}


**Method 2:** Using the update method

.. code:: ipython3

    docobj = docList[0]
    docobj.update(resource="resource2",desc=newdata2)




.. parsed-literal::

    1



Now we update the object and fetch the current values from the database:

.. code:: ipython3

    after_update_docList = datalayer.Measurements.getDocuments(projectName=projectName,description_A = "A")
    after_update_docobj = docList[0]
    print('The resource is: %s' %after_update_docobj.resource)
    print('The description is: %s' %after_update_docobj.desc)


.. parsed-literal::

    The resource is: resource2
    The description is: {'description_A': 'A', 'description_B': 'B', 'description_C': 'C2'}


However, the docobj still retains it old value.

.. code:: ipython3

    print('The after_update_docobjresource is: %s' %docobj.resource)


.. parsed-literal::

    The after_update_docobjresource is: resource1


In order to refresh the object in memory (i.e to reload from the
database), use the reload function:

.. code:: ipython3

    docobj.reload()




.. parsed-literal::

    <Measurements: Measurements object>



.. code:: ipython3

    print('The after_update_docobjresource is: %s' %docobj.resource)


.. parsed-literal::

    The after_update_docobjresource is: resource2


.. code:: ipython3

    copy1_docList = datalayer.Measurements.getDocuments(projectName=projectName,description_A = "A")
    copy2_docList = datalayer.Measurements.getDocuments(projectName=projectName,description_A = "A")
    
    copy1_doc0 = copy1_docList[0]
    copy2_doc0 = copy2_docList[0]

Now we update the first document in copy1_docList

.. code:: ipython3

    copy1_doc0.update(resource="resource3")




.. parsed-literal::

    1



The objects re remains the old value:

.. code:: ipython3

    print('Resource in copy 1: %s' % copy1_doc0.resource)
    print('Resource in copy 2: %s' % copy2_doc0.resource)


.. parsed-literal::

    Resource in copy 1: resource2
    Resource in copy 2: resource2


To update copy2 we need to reload it.

.. code:: ipython3

    copy2_doc0.reload()




.. parsed-literal::

    <Measurements: Measurements object>



Now the values in the instance copy2_doc0 are synchronized with the DB

.. code:: ipython3

    print('Resource in copy 1: %s' % copy1_doc0.resource)
    print('Resource in copy 2: %s' % copy2_doc0.resource)


.. parsed-literal::

    Resource in copy 1: resource2
    Resource in copy 2: resource3


Deleting a metadata document.
-----------------------------

There are 2 ways to delete documents from the DB.

The first deletes documents with the collection, and actually allows for
deletion of all the documents that satisfy a criteria (see `querying the
database <#another_cell>`__). The other method deletes one document
using the document object.

Note that the deletion deletes only the database document, and **not**
actual file on the disk (if a file is a resource).

In order to use this example, lets add 2 more documents to the database.
One with ‘description_A’ equals ‘A’ and the other with ‘description_A’
equals ‘B’.

.. code:: ipython3

    datalayer.Measurements.addDocument(projectName=projectName, type=documentType, dataFormat=dataFormat, resource=resource, desc=dict(description_A="A", description_B="C"))
    datalayer.Measurements.addDocument(projectName=projectName, type=documentType, dataFormat=dataFormat, resource=resource, desc=dict(description_A="B", description_B="D"))




.. parsed-literal::

    <Measurements: Measurements object>



Deleting documents using a query
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to delete all the documents that satisfy a criteria we use the
deleteDocuments method. For example, delete the documents for which
‘description_A’ equals ‘A’

.. code:: ipython3

    datalayer.Measurements.deleteDocuments(projectName=projectName,description_A='A')




.. parsed-literal::

    [{'_id': {'$oid': '63d76cf44afad317543e2840'},
      '_cls': 'Metadata.Measurements',
      'projectName': 'addDataExample',
      'desc': {'description_A': 'A', 'description_B': 'B', 'description_C': 'C2'},
      'type': 'ExampleData',
      'resource': 'resource3',
      'dataFormat': 'JSON_pandas'},
     {'_id': {'$oid': '63d76cf44afad317543e2841'},
      '_cls': 'Metadata.Measurements',
      'projectName': 'addDataExample',
      'desc': {'description_A': 'A', 'description_B': 'C'},
      'type': 'ExampleData',
      'resource': '{"date":{"1478815200000":"2016-11-10","1478818800000":"2016-11-10","1478822400000":"2016-11-11","1478826000000":"2016-11-11","1478829600000":"2016-11-11","1478833200000":"2016-11-11","1478836800000":"2016-11-11","1478840400000":"2016-11-11"},"time":{"1478815200000":"22:00:00","1478818800000":"23:00:00","1478822400000":"00:00:00","1478826000000":"01:00:00","1478829600000":"02:00:00","1478833200000":"03:00:00","1478836800000":"04:00:00","1478840400000":"05:00:00"},"value":{"1478815200000":90,"1478818800000":91,"1478822400000":80,"1478826000000":87,"1478829600000":84,"1478833200000":94,"1478836800000":91,"1478840400000":94}}',
      'dataFormat': 'JSON_pandas'}]



Deleting a single document:
~~~~~~~~~~~~~~~~~~~~~~~~~~~

To delete a single document, we get the document using the retrieve and
then delete it.

.. code:: ipython3

    document_to_delete = datalayer.Measurements.getDocuments(projectName=projectName,description_A = "B")

Now, we delete the document

.. code:: ipython3

    document_to_delete.delete()




.. parsed-literal::

    1



.. code:: ipython3

    docList = datalayer.Measurements.getDocuments(projectName=projectName,description_A = "B")

Now the database is empty

.. code:: ipython3

    print(docList)


.. parsed-literal::

    []

