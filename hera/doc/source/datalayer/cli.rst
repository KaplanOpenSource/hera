Command line interface (CLI)
****************************

The CLI command allows the user to perform actions on the data.

.. code-block:: bash

    hera-datalayer [operation] [arguments]

Where operations are:

- list: list all the data.
- delete: Delete documents from the database with/without deleting data from the disk.
- load: Load documents to the database from JSON
- copyTo: copy data from or to the other databases.

Listing metadata documents
==========================

with the following command lists all the documents according to

.. code-block:: bash

    hera-datalayer list [projectName] <query 1> <query 2> ... <query N>

If project name is not present, then assume all the arguments are query.

Each query pharse is a mongoDB query string.
 That is [field name] = value, or more complicated queries.

Example
-------

For example in order to list all the experimental metadata of station 'Check_Post' in project 'Haifa':

hera-datalayer list projectName="'Haifa'" type="'meteorological'" station="'Check_Post'"

Counting metadata documents
===========================

In order to count the number metadata documents according to a wanted query you should use the CLI(command line interface) with the following command:

.. code-block:: bash

    hera-datalayer count [projectName] <query 1> <query 2> ... <query N>

If project name is not present, then assume all the arguments are query.

Each query pharse is a mongoDB query string.
 That is [field name] = value, or more complicated queries.


Deleting metadata from the database
===================================

In order to delete metadata from the database you should use the CLI(command line interface).

There are 2 steps :
1. Query the metadata you want to delete.
2. Delete the result of the query.

1. Query the metadata
---------------------

In order to query the metadata from the database you should use the following command:

.. code-block:: bash

    hera-datalayer delete <query>

The query should filter the metadata documents you want to delete.
After you execute the command, it will save a json file called 'docsToDelete.json' with those metadata documents to your current directory.

Note that this json file is in the same format as the json for loading metadata through the CLI.

2. Delete the results of the query
----------------------------------

In order to delete the results of the query you should use the following CLI command:

.. code-block:: bash

    hera-datalayer delete <JSON>

The <JSON> should be the json file of the results of the query which called 'docsToDelete.json' unless you changed its name.

Example
^^^^^^^

For example in order to delete all the experimental metadata of station 'Check_Post' in project 'Haifa' :

1. Query the metadata

.. code-block:: bash

    hera-datalayer delete projectName="'Haifa'" type="'meteorological'" station="'Check_Post'"

2. Delete the results of the query

.. code-block:: bash

    hera-datalayer delete docsToDelete.json



Loading metadata to the database
================================

Loading data into the database from the CLI is done by the following command:

.. code-block:: bash

    hera-datalayer load docsToLoad.json

docsToload.json should contain the metadata to load, and should be in the form of:

.. code-block:: javascript

    {
    "Measurements": [{"projectName='myProject',
                      "resource='path-to-parquet-files',
                      "dataFormat='parquet',
                      "type='meteorology'
                      "desc={"property1": "value1",
                             "property2": "value2"
                             }
                      },
                     measurementsDoc2,
                     .
                     .
                     .
                     ],
    "Simulations": [simulationsDoc1,
                    simulationsDoc2,
                    .
                    .
                    .
                    ],
    "Analysis": [analysisDoc1,
                 analysisDoc2,
                 .
                 .
                 .
                 ]
    }


Note that a metadata document should contain the following keys:

#. "projectName"
#. "resource"
#. "dataFormat"
#. "type"
#. "desc"


Copy metadata documents To/From other database
==============================================

In order to copy metadata documents to/from others databases you should use the CLI(command line interface) with the following commands:

1. Copy to
----------

hera-datalayer copyTo <others database login info> <query 1> <query 2> ...  <query N>

Each query pharse is a mongoDB query.

2. Copy from
------------

.. code-block:: bash

    hera-datalayer <others database login info> <query 1> <query 2> ...  <query N>

Where:

- <others databse login info> should be username:password@IP/dbName

        username, password and dbName belongs to the person you want to copy from/to.
        IP is the ip where the database.

- <query i> are query terms.