OpenFOAM
========


### Note! the new  Hermes workflow (previously pipeline) and hera pipeline usage is in
'Run Hermes and Hera to get new OF simulation' bellow. The relevant patrs in the old  usage should be merge into it ###

OpenFOAM simulations results can be managed using hera.
The results can be converted to convinient pandas dataframe and be analysed
using built-in plotting functions.

Plotting can be done using the plotting tool.
Examples can be find in hera-gallery. The tool is called plotting.


.. code-block:: python

    from hera import openfoam
    plotting = openfoam.Plotting()

Usage
-----

.. toctree::
    :maxdepth: 2

    openFoam/StepByStep/HermesHeraToOF
    openFoam/ManagingData
    openFoam/analysis
    openFoam/utils

    openFoam/preProcess
    openFoam/meshUtilities
    openFoam/tests
    openFoam/api
