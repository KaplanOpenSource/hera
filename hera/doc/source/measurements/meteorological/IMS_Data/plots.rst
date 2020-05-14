***********
IMS plots
***********

The imsdata plots classes arrange the use of some conventional meteorological plots.
The main plot class organizes the common settings for all plots (colors, labels, axes, etc.).
Then, sets of plots are grouped into a subplot classes by common characteristics:

- Daily plots: A single plot of a meteorological variable during the day (24 hours). The plot can be single day data, longer period data, or a period statistics. The data will be drawn vs. the time during the day.

- Seasonal plots: A 2x2 plot of a meteorological variable during the day (24 hours), grouped by season. Each data set will be grouped by season. Each season's data will be drawn in a sub-plot within the main plot, and will be drawn vs. the time during the day

For quick use, we recommend using Hera's gallery

.. toctree::
    :maxdepth: 1

    plots/daily
    plots/seasonal