****************************
Daily plots
****************************

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


