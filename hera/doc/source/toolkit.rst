.. _toolkitPage:

Toolkits
********

Overview
========

A Toolkit is collection of procedures that are used work with a dtat of a certain type
(like Urban GIS, topography, simulations and ect.).

The structure and function of  toolkits
========================================

Because all the data in hera is attached to project, it is necessary
to supply a projectName when initializing one. Then, the functions
of the toolkit will be performed on the data that exists in that project.

A toolkit is a component that can handle a certain type of data. For example a
BUILDING toolkit will allow the user to use the GIS function that were programmed
in hera, like compting :math:`\lambda_f` or :math:`\lambda_p` of a city.

Generally, each toolkit has 3 components.

    - **Datalayer**: Manages loading and parsing data.
                     Also, manages the addition of the data to the database.

    - **Analysis**: Common analysis and computational functions.

    - **Presentation Layer**: Common graphs and other presentation layer.






How to use toolkits
================

Access to the toolkits is done through the tookitHome class:

.. code-block:: python

    from hera import toolkitHome

Getting a toolkit
-----------------

In order to get a toolkit use the following

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
|:ref:`Hermes workflow <_HermesWorkflow>`    | The hermes workflow toolkit                |
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
   toolkits/simulations/HermesWorkflow
   toolkits/riskAssessment


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`



