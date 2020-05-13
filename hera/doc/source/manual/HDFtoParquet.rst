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