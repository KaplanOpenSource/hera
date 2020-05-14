.. pyhera documentation master file, created by
   sphinx-quickstart on Tue Nov 26 15:51:03 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyhera's documentation!
==================================

The 'Hera' project aims to make the departmental data more accessible and to standardize its storage,
and provide tools for analysis and presentation.

The library stores the metadata (data that describes the data itself) of
each piece of data as a document in the database. The document that holds the
metadata includes:
 - The project name that contains the data.
 - The data (or a reference to it).
 - The format of the data on the disk.
 - fields that describe the data.

In order to organize the data, the library is comprised of 4 parts.


- Datalayer  : A general layer that helps to manage the metadata records.
               It can be used directly through other layers.

- Measurements : This layer manages data that originates from measurements.
                 Contains several specializations:

                 - GIS: Manages GIS data. Allows to slice and store regions
                        as well as raster data.

                 - Meteorological: Manages meteorological data, provide tools to analyze:
                     - Raw sonic data
                     - Data or mean data from IMS.

- Simulations : Data that originates from simulations.
                Contains several specializations:

                - OpenFOAM: manage openFOAM simulations.
                - LSM: Manage Lagrangian stochastic models.


.. toctree::
   :maxdepth: 3
   :caption: Contents:

   datalayer
   measurements
   simulations
   instructions_for_dev
   how_to
   auto_examples/index


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
