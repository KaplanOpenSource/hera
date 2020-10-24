.. pyhera documentation master file, created by
   sphinx-quickstart on Tue Nov 26 15:51:03 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyhera's documentation!
##################################

Overview
========

The 'Hera' project aims to make the departmental data more accessible and to standardize its storage,
and provide tools for analysis and presentation.

The library manages data by providing interface to store pieces of data
and their metadata (i.e description of the data). Then, the library
provides tools that help to load, analysis and presentations for specific data types.


.. image:: Hera-Layers.png



Datalayer
   Tools to stome the metadata (data that describes the data itself) of each piece of data as a document in the database.
   A general layer that helps to manage the metadata records.
               It can be used directly through other layers.

Cache Caching data after long computation.

Measurements
   Provide tools to manage data that originates from measurements.

   Contains two groups of tools

   - GIS
         Manages GIS data and has the following tools:
         - buildings
            Manages GIS layer of buildings. Allows slicing and performing calculations of
            the porosity and conversion to STL.

         - image
            Manages Images of locations. Saves the image with its dimensions.

         - shape
            Manages spatial polygons and data.

         - topography
            Manages the topography. Slicing areas from the global database (such as BNTL),
            computing heights, and converting to STL.

   - Meteorological
      Manages meteorological data, provide tools to analyze:
      - High frequency data
         This refers to the raw data from the devices.
         Currently we have implemented:

         - Sonic data
               Analyze and calculate turbulence properties of raw sonic data (20-30Hz).
               Loads data from several formats.

         - Temerature data
               Analyze and calculate raw data from temperature measurements (1Hz).


      - Low frequency data
            Usually data that was processes and averaged to several minutes.
            Loads data from several formats.

Simulations
   Manages and analyzes data that originates from simulations.

   Contains several tools:

       - OpenFOAM
               Manage openFOAM simulations.

       - LSM
            Manage and run Lagrangian stochastic models.



Risk assessment
   Provides tools to perform risk assessment. That is project the concentration and estimate the consquences
   on the health of the population.

Installing & setup
==================

We currently recommend using the package as a development package.
You must have MongoDB installed.

# 1. Download the package from the git repository.

::

     git clone [path to repository] path-to-save/pyhera

# 2. Add the path to the to the PYTHONPATH.

::

   export PYTHONPATH=$PYTHONPATH:path-to-save/pyhera/hera

   we recommend to add it to the .bashrc.

#3. Create a configuration file in $HOME/.pyhera

.. literalinclude:: examples/configfile/config.json

example configuration file  can be downloaded :download:`here <examples/configfile/config.json>`.

The default database name is the linux username and it **must** exist.


Usage and API
=============


.. toctree::
   :maxdepth: 3
   :caption: Contents:

   datalayer
   measurements
   simulations
   riskassessment
   instructions_for_dev
   how_to
   auto_examples/index


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
