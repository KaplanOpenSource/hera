WRF
===

The WRF module is used to manage wrf simulation results.
It helps converting the results into a convinient pandas dataframe.


10-minute tutorial
------------------

The module's main function is called getPandas.
It converts a requested part of the wrf results into a pandas dataframe.
The requested part is a set of points closest to a requested longitude or lantitude.
The longtitude or lantitude can be given in either ITM or wgs84 units.

For example, suppose we are intersted in the longitude 732000.
We can get the data from the wrf simulation this way:


.. code-block:: python

    from hera import WRF
    wrfDatalayer = WRF.wrfDatalayer()
    path = "somepath"
    longtitude = 200000

    data = wrfDatalayer.getPandas(datapath=path, lon=longitude)

Additional examples for more abilities of the function may be found in the Get Data page.

Usage
-----

.. toctree::
    :maxdepth: 1

    WRF/GetData