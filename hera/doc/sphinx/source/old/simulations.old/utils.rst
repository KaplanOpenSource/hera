Utils
=====

Coordinate Handler
------------------

Hera has a class that helps interpulate values of a simulation in new coordinates.
It can be usefull, for example, in order to transform openFOAM simulations results from
the original grid to a regularized grid or "sigma" topography-following coordinates.
These are the two grid types that are supported at the moment.
The class supports two-dimension slices only.

The default grid type is regularized grid, with x and z coordinates,
and variables called "U_x" and "U_z".
It can be achieved this way:

.. code-block:: python

    from hera import coordinateHandler
    data = datalayer.Simulations.getDocuments(projectName="Some2DsimulationData")[0].getData()
    dataRegular = coordinateHandler.regularizeTimeSteps(data=data)

The default number of grid points in each direction is 600, but it may be changed.
The minimum and maximum values are the minimum and maximum of the original dataframe,
but other limits may be requested.

In addition, if no time steps are requested, the function assumes that all the data
belongs to the same time steps. Otherwise, timesteps should be handed.

For example, if we wish to use other coordinates, limits, number of points, timesteps
and variables, it can be done this way:

.. code-block:: python

    coord1 = "x"
    coord2 = "y"
    coord1Lim = (0,10) # Min and max values for "x".
    coord2Lim = (0,10) # Min and max values for "y".
    n = (100,100) # number of grid points in x and y directions.
    timelist = [0,1]
    fieldList = ["firstVariable","secondVariable"]
    dataRegular = coordinateHandler.regularizeTimeSteps(data=data, fieldList=fieldList,
                                                        timelist=timelist,n=n,coord1=coord1,
                                                        coord2=coord2, coord1Lim=coord1Lim,
                                                        coord2Lim=coord2Lim)


If we would like to use sigma coordinates, we have two options.
If the dataframe has a column with the heights of the topography, it can be done like this:

.. code-block:: python

    dataRegular = coordinateHandler.regularizeTimeSteps(data=data, coordinateType="Sigma",
                                                        terrain="topographyColumnName")

The function uses the original x coordinate values as x coordinate grid points, and
creates by default 10 sigma values between 0 and 1.

If the topography is not defined, the user should first build a dataframe with the topography's
x and z coordinates values. Then, the x values of this dataframe would be used for the new,
sigma-coordinates grid.

For example, for transforming the data to sigma coordinates with a hundred sigma values,
in the case of a sinusoid wave:

.. code-block:: python

    terrain = pandas.DataFrame({"x":[x for x in range(200)],"z":[numpy.sin(x*numpy.pi/100) for x in range(200)]})
    dataRegular = coordinateHandler.regularizeTimeSteps(data=data, coordinateType="Sigma",
                                                        terrain=terrain,n=100)
