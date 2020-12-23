LSM
===

The pre-process tool is used in order to load openFoam LSM simulation results to
the database, and extracting the concentration from the results.
In addition, it may be used to create sources of the particles.

.. code-block:: python

    from hera.simulations.openfoam.lsm.preProcess import preProcess
    pre = preProcess(projectName="LSMOpenFOAM",casePath = "The case path", cloudName = "The cloud name")

When initializing the tool, the project name, the path of the LSM case and the name of the cloud, defined in the case,
must be given.
They may be changed later. for example, one may get the current cloud name and update it this way:

.. code-block:: python

    oldCloudName = pre.cloudName #gets the old name
    pre.cloudName = "newCloudName" #updates to the new name
    
Getting the concentration
-------------------------

The concentration is achieved using the following function:

.. code-block:: python

    from unum.units import *
    Concentration = pre.getConcentration(endTime=100,file="somedirectory/somefilename.parquet"
                                         startTime=1,Q=1*kg, dt=10*s,dx=10*m,dy=10*m, dz=10*m,
                                         Qunits=mg,lengthUnits=m,timeUnits=s,nParticles=None,save=False,addToDB=True)

endTime is the final time step to read.
All other parameters are optional; the values above are the default value, except for file, whose default is None.

Q is the total mass of the particles.
In case of an instantaneous release, one should deliver a value in units of mass, as shown above.
In case of a continueous release with constant rate, one may either deliver a value in units of mass or mass over time.
If the units are of mass, this is the mass which is released every time step.
If the units are mass over time, Q defines the rate of the release.
In case of a continueous release which is not constant with time, one should deliver a list of values in units of mass,
which is equal in length to the number of time steps.

The dx, dy, dz, dt parameters are the dimensions of the cells and time steps that are used
to calculate the concentration.
The "units" parameters are the units of the concentration (Qunits/(lengthUnits**3)
and time in the final dataframe.

nParticles is the number of particles released in a single time step. If nParticles is None, the number of particles is found using the kinematicCloudPositions file;
if the function can't find the file it may ask the user to deliver the number to the function.

If save is set to True, the concentration is saved under the file parameter name.
If this parameter is None, it is saved in the current working directory under the name "Concentration.parquet".
If addToDB is True, the data is added to the database. Any aditional descriptors may be given as **kwargs.

Loading positions and velocities
--------------------------------

The locations and velocities of the particles in an OpenFOAM LSM simulations are extracted
to a dask dataframe using the following function:

.. code-block:: python

    file = "/home/ofir/Projects/2020/LSM/horizontal500000/results.parquet"
    data = pre.extractRunResult(times=[i for i in range 100], withVelocities=False, file=file, save=False, addToDB=True)

All parameters are optional.
times is a list of time steps to read, in this case, the first hundred.
If it is not given, all time steps in the casePath are extracted.
By default, only the positions of the particles are extracted.
Their velocities may be extracted too by setting the withVelocities parameter to True.

The save and addToDB parmaeters may be set to True in order to save the
dataframe to the disc and add it to the database.
The new file name and path may be given as the file paramters. The default path is the
current working directory and the default name is "Positions.parquet".
Any additional descriptors for the new document in the database may be given as **kwargs.

Creating source files
---------------------

As of today, only instantaneous, point source files can be created using the tool.
They are saved under the name "kinematicCloudPositions" at the constant subdirectory of the case path.
They can be created this way:

.. code-block:: python

    x = 10 # the x coordinate
    y = 10 # the y coordinate
    z = 10 # the z coordinate
    nParticles = 100000 # the number of particles
    pre.makePointSource(x,y,z,nParticles)