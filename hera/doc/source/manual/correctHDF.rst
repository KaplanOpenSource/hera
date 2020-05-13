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

