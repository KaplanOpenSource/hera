OPENFOAM utilities
==================

extracting cell centers data
----------------------------

If one needs to extract the cell center locations from the openfoam run, the 'centersToPandas' function can be used.

before using it, run this command using the command line at the appropriate folder:

.. code-block:: python

    postProcess -func writeCellCentres

The cell centers data is now extracted to C, Cx, Cy and Cz files.
now we can convert it to pandas dataframe using centersToPandas function:

.. code-block:: python

    from hera import openfoam
    celldata=openfoam.centersToPandas()


see: :ref:`utils-reference-label`
