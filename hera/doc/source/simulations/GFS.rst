GFS
===

GFS data
********
GFS data can be taken online for the last two weeks. it based on a global forecast model with 0.25 degrees resolution.
The data are located at:
https://www.nco.ncep.noaa.gov/pmb/products/gfs/ and should be downloaded manually (we have to write a script that will download it)

To get the list of the fields

.. code-block:: python

    mygfs = gfs()
    filename2 = '/ibdata2/nirb/testme.grib'
    gfs.get_gfs_list(filename2)
    gfs.get_gfs_data(filename2, lat=31.76, lon=35.21, band_num=442)

* We still have to validate the coordinate to index code at get_gfs_data
* The GFS has wind / temperature data at different levels and get can get the profile out of it


