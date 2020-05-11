****************************
IMS Analysis layer
****************************

The imsdata analysis class organizes the data post process.
The class provides some tools for meteorological data processing as well as tools for organizing the data in a way
that is required for further work.



Add date and time columns
-------------------------

addDatesColumns function adds date and time columns to the dataframe from the date index column. These fields can,
for example, when you want to extract data by different time slices (by season, year, hour, etc.)

.. code-block:: bash

    data=docList[0].getData()
    data=meteo.IMS_analysis.addDatesColumns(data)

By default, the function takes the 'index' column as the source for calendar information.
The user can give his own date and month columns using the 'datecolumn' and 'monthcolumn' variables.
In addition, the default function handles dask-type data. You can change the type to 'pandas' with the 'dataType' variable.

see: :ref:`Analysis-reference-label`
