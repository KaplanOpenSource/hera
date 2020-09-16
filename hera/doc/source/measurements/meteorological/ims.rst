IMS Data
++++++++

The imsdata module designed to provide tools required to work with data obtained
from Israel Meteorological Service (IMS).
It used for loading and reading IMS data from and to the the database,
provide some analysis tools and helps to ease meteorological data plots

Some of the tools in the module also allow you to work with data coming from other sources than IMS,
where its structure is similar (meteorological log with timestamp)
For uniformity, the different uses of the 'imsdata' module will be with the prefix **'IMS_'**

Datalayer
---------

The imsdata datalayer class organizes the data loading into the database
or reading data from the database. all databse properties is defined for a meteorological data format from ims source.

Load data
=========

First, we need to import meteo from 'Hera':

.. code-block:: bash

    from hera import meteo


The LoadData function organizes the data upload into the database with the appropriate format and structure.

Since this is particular data type, the following parameters are set by default for files coming from the IMS database
(of course they can be overwrite by the user):

- type='meteorological'
- DataSource='IMS'
- station_column='stn_name'
- time_coloumn='time_obs'

When calling this function, the following parameters must be set:

- **newdata_path:** The path to the folder with the new raw data
- **outputpath:** The path to the folder where the new arranged data will be saved
- **Projectname:** A project name string associated to th new loaded data. When the data is general or not explicitly linked to a particular project,we recommend using the general project 'IMS_Data'


When the data is associated with a metadata file, it can be loaded by passing path to the file with **'metadatafile'** variable.

Here is an example of how to use the LoadData function:

.. code-block:: bash

    meteo.IMS_dataloader.LoadData(newdata_path = '/raid/users/davidg/Project/2020/IMS_RAW/json/multiStations',
                                  outputpath = '/raid/users/davidg/Project/2020/IMS_RAW/parquet/',
                                  Projectname = 'IMS_Data',
                                  metadatafile = '/raid/users/davidg/Project/2020/IMS_RAW/meta.txt'
                                  )

The function performs the following procedures:

1. Converts the new data from json file to dask format
2. Adds the new data to the old data (if exists) by station name (each station saved to separate file)
3. Saves the unified data to a parquet files in the request 'outputpath' location
4. Saves new document to the metadata with the metadata added to the description (for new stations)


Analysis
--------

The imsdata analysis class organizes the data post process.
The class provides some tools for meteorological data processing as well as tools for organizing the data in a way
that is required for further work.


Add date and time columns
=========================

addDatesColumns function adds date and time columns to the dataframe from the date index column. These fields can,
for example, when you want to extract data by different time slices (by season, year, hour, etc.)

.. code-block:: python

    data=docList[0].getData()
    data=meteo.IMS_analysis.addDatesColumns(data)

By default, the function takes the 'index' column as the source for calendar information.
The user can give his own date and month columns using the 'datecolumn' and 'monthcolumn' variables.
In addition, the default function handles dask-type data. You can change the type to 'pandas' with the 'dataType' variable.


IMS plots
--------

The imsdata plots classes arrange the use of some conventional meteorological plots.
The main plot class organizes the common settings for all plots (colors, labels, axes, etc.).
Then, sets of plots are grouped into a subplot classes by common characteristics:


****************************
Daily plots
****************************

A single plot of a meteorological variable during the day (24 hours). The plot can be single day data, longer period data, or a period statistics. The data will be drawn vs. the time during the day.
**When calling one of the plot functions, many properties related to the type of drawing or to th axes can be controlled and changed (line color or mark, with or without scatter, contour levels, labels, etc.). It is highly recommended to look at the API for complete and detailed documentation**

some of the following functions by default combine some of the plot functions

First, we need to import meteo from 'Hera':

.. code-block:: bash

    from hera import meteo


Line Plot
=====================
Make a simple line plot for selected plot field and date.

When calling this function, the following parameters must be set:

- **Data:** The data set to plot
- **plotField:** A string of the field to plot. must be one of the column names from data
- **date:** The date (in format of 'YYYY-MM-DD') to plot.

Here is an example of how to use the dateLinePlot function:

.. code-block:: bash

    Data = MyDataSet
    plotField='plotFieldColumnName'
    dateToPlot='2019-01-01'

    meteo.IMS_dailyplots.dateLinePlot(Data, plotField = plotField, date = dateToPlot)


probability Contourf Plot
=========================
Make a probability contour plot for selected plot field

When calling this function, the following parameters must be set:

- **Data:** The data set to plot
- **plotField:** A string of the field to plot. must be one of the column names from data

Here is an example of how to use the plotProbContourf function:

.. code-block:: bash

    Data = MyDataSet
    plotField='plotFieldColumnName'

    meteo.IMS_dailyplots.plotProbContourf(Data,plotField)


Scatter Plot
=====================
Make a scatter plot for selected plot field.

When calling this function, the following parameters must be set:

- **Data:** The data set to plot
- **plotField:** A string of the field to plot. must be one of the column names from data

Here is an example of how to use the plotScatter function:

.. code-block:: bash

    Data = MyDataSet
    plotField='plotFieldColumnName'

    meteo.IMS_dailyplots.plotScatter(Data, plotField)



****************************
Seasonal plots
****************************

A 2x2 plot of a meteorological variable during the day (24 hours), grouped by season.
Each data set will be grouped by season.
Each season's data will be drawn in a sub-plot within the main plot, and will be drawn vs. the time during the day

**When calling one of the plot functions, many properties related to the type of drawing or to th axes can be controlled and changed (line color or mark, with or without scatter, contour levels, labels, etc.). It is highly recommended to look at the API for complete and detailed documentation**


Probability conturf Plot
========================
Make a 2x2 probability contour plot for selected plot field by season

When calling this function, the following parameters must be set:

- **Data:** The data set to plot
- **plotField:** A string of the field to plot. must be one of the column names from data

Here is an example of how to use the plotProbContourf function:

.. code-block:: bash

    Data = MyDataSet
    plotField='plotFieldColumnName'

    meteo.IMS_seasonplots.plotProbContourf_bySeason(Data, plotField)





*************
imsdata API
*************


.. _datalayer-reference-label:

datalayer class
###############

.. autoclass:: hera.measurements.meteorological.imsdata.DataLayer
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



