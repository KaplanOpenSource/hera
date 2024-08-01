
Python
*******

Formatting string
=================

    In python >= 3.6 it is possible to format a string with the 'f""' directive:

.. code-block:: python

    a = 5
    s = f"The value of a is {a}"


will result in ``s`` holding ``"The value of a is 5"``

Python 3.8 added the ``=`` suffix for debug printings:

.. code-block:: python

    a = 5
    s = f"Now {a=}"

will result in an error with Python < 3.8, but newer
Pythons will give ``"Now a=5"``

Debug logging
=============

To enable debug in specific loggers, see
:py:func:`hera.utils.logging.helpers.initialize_logging`.

.. automodule:: hera.utils.logging.helpers


Creating objects during code execution
=======================================

Sometimes (very rarely) it is important to derive objects in realtime.
In Hera we use it to create object for the ORM (object relation mapping) of MongoDB.

For example to create a new class that derives form myFather:

.. code-block:: python

    newClass = type('mynewclass', (Son,Father,), {})

    newClassInstance = newClass()

is Equivalent to

.. code-block:: python

    class mynewclass(Son,Father):
        ...


How to convert HDF file to Parquet format:
============================================

 like this:

.. code-block:: python

    import pandas
    import os
    import dask.dataframe
    datafold='/data3/Campaigns_Data/hdfData/2006-12_IMS/'
    newdatafold='/data3/Campaigns_Data/parqueData/IMSdata/'
    datafilename='IMS.h5'
    if not os.path.exists(newdatafold):
        os.makedirs(newdatafold)

    np=1

    with pandas.HDFStore (datafilename) as file:
        keys = [x for x in file]

    del file

    for key in keys:
        newfold = os.path.join(newdatafold, key[1:])
        if not os.path.exists(newfold):
            print(key)
            os.makedirs(newfold)

            tmppandas=pandas.read_hdf(datafilename,key=key)
            tmpdata = dask.dataframe.from_pandas(tmppandas, npartitions=np)
            FileNameToSave=key[1:]+'.parquet'
            tmpdata.to_parquet(os.path.join(newfold,FileNameToSave), engine='pyarrow')
        else:
            print(key + ' fold exist')

How to correct HDF file data with bad/broken data cells:
========================================================
If your data have a mix of float and string data (for any reason) the '.parquet' save will fail.
To solve it, you should replace all bad values with 'None'.

First, find all bad values in the mixed column data (it will appear as an error when the '.parquet' saving process will fail).
In the example, saving failed in the following stations and columns:

.. code-block:: python

    badStationList=['Metzoke_Dragot','Hazeva','Tavor_Kadoorie']
    badcolumns=['WS1mm','Grad','WSmax']

you can use this small function to list all bad values:

.. code-block:: python

    def convetToFloat(x):
         try:
             float(x)
         except:
             print("---%s---" % x)

and then:

.. code-block:: python

    datafold='/data3/Campaigns_Data/hdfData/2006-12_IMS/'
    datafilename='IMS.h5'

    for badStation,col in product(badStationList,badcolumns):
    tmppandas = pandas.read_hdf(os.path.join(datafold,datafilename),key=badStation)
    if col in tmppandas.columns:
        tmppandas[col].apply(convetToFloat)

the bad values will appear in the head of the print. collect all bad values in list:

.. code-block:: python

    badvalues=['ST_A','ST_B','ST_C','ST_D','<Samp']

now, we will find all bad values and replace them with 'None'

.. code-block:: python

    newdatafold='/data3/Campaigns_Data/parquetData/IMSdata/'
    np=1

    for stn in badStationList:
        tmppandas = pandas.read_hdf(os.path.join(datafold,datafilename), key=stn)
        newfold = os.path.join(newdatafold, stn)
        if not os.path.exists(newfold):
            print(key)
            os.makedirs(newfold)
        for val,col in product(badvalues,badcolumns):
            print (val,col)
            if col in tmppandas.columns:
                tmppandas[col]=tmppandas[col].replace(val,None)

change column type to float:

.. code-block:: python

        for col in badcolumns:
            tmppandas[col]=tmppandas[col].astype(float)
        tmpdata = dask.dataframe.from_pandas(tmppandas, npartitions=np)

and save the data:

.. code-block:: python

        FileNameToSave = stn + '.parquet'
        tmpdata.to_parquet(os.path.join(newfold, FileNameToSave), engine='pyarrow')


