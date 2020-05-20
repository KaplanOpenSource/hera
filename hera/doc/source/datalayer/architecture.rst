Datalayer Architecture
======================

The library uses the mongoengine library to access the mongoDB.

The library uses a object relational mapping (ORM) technology. That is,
One object represents the collection (a table) and each record
is represented as an object.

All the data is stored in the Metadata collection.

Each record has the following structure:

.. list-table:: Record

    * - Name
      - Description

    *  - projectName
       - The project that contains this record



Main classes
------------

The library has 2 main classes. The collection objects (defined in collection.py)
and the MetadataFrame object defined in document.metadataDocument and represents a single record.

Collection objects
^^^^^^^^^^^^^^^^^^

The collection objects manage the access to the documents.
Using these objects it is possible to query and add the documets.

When the collection objects are initialized they are related to a connection.
The default connection is the name of the user.

- Measurements_Collection
- Simulations_Collection
- Cache_Collection

In order to relate to anther user, use

.. code-block:: python

    measurement_newDB = Measurements_Collection(user="public")


MongoDB query language
----------------------





