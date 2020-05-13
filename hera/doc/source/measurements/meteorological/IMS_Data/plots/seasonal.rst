****************************
Seasonal plots
****************************
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





