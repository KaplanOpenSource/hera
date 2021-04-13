.. _toolkitPage:

Toolkits
********

Toolkit is a code to load, analyze and present a certain type of data.

Generally, each toolkit has 3 components.

    - **Datalayer**: Manages loading and parsing data.
                     Also, manages the addition of the data to the database.

    - **Analysis**: Common analysis and computational functions.

    - **Presentation Layer**: Common graphs and other presentation layer.



Managing toolkits
================

Access to the toolkits is done through the tookitHome class:

.. code-block:: python

    from hera import toolkitHome

Getting a toolkit
-----------------

.. code-block:: python

    projectName = "The-Project-Name"
    newtoolkit  = toolkitHome.get(toolkitName=toolkit-name,
                                   projectName=projectName,
                                  ... extra toolkit specific parameters ...)


Listing toolkits
----------------
    In order to list the available


+--------------------------------------------+--------------------------------------------+
|           .. centered:: GIS                                                             |
+--------------------------------------------+--------------------------------------------+
|Toolkit                                     |    Usage                                   |
+============================================+============================================+
|:ref:`Buildings <BuildingsPage>`            |    Manages building shapes (vector data)   |
+--------------------------------------------+--------------------------------------------+
|:ref:`Shapes   <ShapesPage>`                |    Manages location shapes (vector data)   |
+--------------------------------------------+--------------------------------------------+
|:ref:`Raster <RasterPage>`                  |    Manages Raster of geogrpaphic places    |
+--------------------------------------------+--------------------------------------------+
|:ref:`Topography <TopographyPage>`          |    Manages topography      (vector data)   |
+--------------------------------------------+--------------------------------------------+
|:ref:`Demography <DemographyPage>`          |    Manages demography      (vector data)   |
+--------------------------------------------+--------------------------------------------+
|           .. centered:: Simulations                                                     |
+--------------------------------------------+--------------------------------------------+
|:ref:`LSM <LSMPage>`                        | Stochastic lagrangian simulation (fortran) |
+--------------------------------------------+--------------------------------------------+
|           .. centered:: Risk assessment                                                 |
+--------------------------------------------+--------------------------------------------+
|:ref:`Risk assessment <RiskAssessmentPage>` |    Estimating the effects of dispersion    |
+--------------------------------------------+--------------------------------------------+


Toolkit architecture
=====================

Toolkits has datasource

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   toolkits/GIS_Buildings
   toolkits/GIS_Raster
   toolkits/GIS_Demography
   toolkits/GIS_Shape
   toolkits/GIS_Topography
   toolkits/LSM
   toolkits/riskAssessment


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`



