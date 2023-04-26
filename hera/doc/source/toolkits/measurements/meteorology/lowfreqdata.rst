Low Frequency data tool
=======================

The low frequency tool is designed to load, analyze and present meteorological data that was
averaged over several minutes.




Accessing data
^^^^^^^^^^^^^^^^

.. include:: lowfreqdata/lowfreqdata.rst


Analysis
^^^^^^^^^^^^^^^^

.. include:: lowfreqdata/Analysis.rst


Plots
^^^^^^^^^^^^^^^^

See the gallery for details of the available plots

Supported formats
^^^^^^^^^^^^^^^^

Currently supported formats are:

* IMS-JSON format.

If the data is stored to the database, then it converted to parquet format.


Low frequency data API
^^^^^^^^^^^^^^^^^^^^^^^


.. _datalayer-reference-label:

datalayer class
###############

.. autoclass:: hera.measurements.meteorology.lowfreqdata.datalayer.DataLayer
    :members:
    :undoc-members:


.. autoclass:: hera.measurements.meteorology.lowfreqdata.parsers.Parser_JSONIMS
    :members:
    :undoc-members:


.. _Analysis-reference-label:

Analysis class
###################


.. autoclass:: hera.measurements.meteorological.imsdata.Analysis
    :members:
    :undoc-members:

.. _Plot-reference-label:

Plot class
##########

.. autoclass:: hera.measurements.meteorological.presentationlayer.dailyClimatology.Plots
    :members:
    :undoc-members:

.. autoclass:: hera.measurements.meteorological.presentationlayer.dailyClimatology.SeasonalPlots
    :members:
    :undoc-members:

.. autoclass:: hera.measurements.meteorological.presentationlayer.dailyClimatology.DailyPlots
    :members:
    :undoc-members:



