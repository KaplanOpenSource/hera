Datalayer Architecture
======================

The datalayer module abstracts the access to the resource of the user.


Each record represents a certain piece of data and thereofre it holds the metadata that describes it
and a reference to the data itself.

The metadata is determined by the user that added that data and stored in the __desc__ property

The data is stores in the __resource__ attribute and its format is stored in the 'dataFormat'
property. A data can be either a link to a file on the disk, or the data itself.

The dataFormat determine how the data is retrived when using the getData function of the record.

When a user adds a data to the database, he describes the metadata and the data itself.
A data that is stored in the database can be either a file or the data itself.



MongoDB query language
***********************





Reference to code
-----------------

.. toctree::
    :maxdepth: 2

    collection
    datahandlers


