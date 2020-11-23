LSM
===

The pre-process tool is used in order to load openFoam LSM simulation results to
the database, and extracting the concentration from the results.

Loading Results
---------------

The locations and velocities of the particles in an OpenFOAM LSM simulations are extracted
to a dask dataframe using the following function:

.. code-block:: python

    from hera.simulations.openfoam.lsm.preProcess import preProcess
    pre = preProcess(projectName="LSMOpenFOAM")
    directory = "/home/ofir/Projects/2020/LSM/horizontal500000" # The directory of the case.
    file = "/home/ofir/Projects/2020/LSM/horizontal500000/results.parquet"
    data = pre.extractRunResult(casePath=directory,cloudName="kinematicCloud",
                                file=file,save=True,addToDB=True)

The casePath and cloudName parameters are mandatory. They are the path of the case
directory and the name of the cloud directories respectively.
The save and addToDB parmaeters may be set to True in order to save the
dataframe to the disc and add it to the database.
The new file name and path may be given as the file paramters. The default path is the
current working directory and the default name is "Positions.parquet".
Any additional descriptors for the new document in the database may be given as **kwargs.

Getting the concentration
-------------------------

The concentration is achieved using the following function:

.. code-block:: python

    from unum import *
    Concentration = pre.getConcentration(data, Q=1*kg, measureTime=10*s,
                                         measureX=10*m,measureY=10*m, measureZ=10*m,
                                         Qunits=mg,lengthUnits=m,timeUnits=s)

The data parameter is a dask or pandas dataframe that was extracted using the
function descussed in the previous section.
All other parameters are optional; the values above are the default value.
Q is the total mass of the particles.
The "measure" parameters are the dimensions of the cells and time steps that are used
to calculate the concentration.
The "units" parameters are the units of the concentration (Qunits/(lengthUnits**3)
and time in the final dataframe.