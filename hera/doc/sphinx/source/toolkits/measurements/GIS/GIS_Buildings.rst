.. _BuildingsPage:


Buildings
=========

A toolkit to manage the buildings and geometric properties that influence the flow.

The toolkit also supports the definition of regions that are 'cut' from a file that contain an extensive region.

Currently, the hera-data contains the BNTL database that encompase all Israel. The toolkit provides tools to cut a sub-region and
store it in the database for future use.


First usage
-----------

If the project requires the use of the buildings then it is necessary
to load the BNTL datasource to the project.

To do so use the line

.. code-block::

    hera-data-load <project Name> Buildings BNTL

In the hera-data project directory.

Usage
-----


.. toctree::
    :maxdepth: 3
    :caption: Contents:

    Buildings/arrangingData
    Buildings/ConvexPolygons


API
---

.. autoclass:: hera.measurements.GIS.vector.buildings.toolkit.BuildingsToolkit
    :members:
    :undoc-members:
    :inherited-members:
