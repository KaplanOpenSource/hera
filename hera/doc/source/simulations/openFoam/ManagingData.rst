Managing Data
=============

This page shows the application of the dataManipulations module.

Arranging data
--------------

There are built-in functions that arrange slices that are perpendicular to the z axis or
3D clips. Their main advantage is that they add the terrain height and height from terrain
to the data. They also add the velocity's magnitude and in the case of slices, the distance downwind.

The slice function can be used like this:

.. code-block:: python

    from hera import datalayer
    from hera import openfoam

    op = openfoam.dataManipulations(projectName="Example")
    slice1 = datalayer.Measurements.getDocuments(projectName="Examples",
             type="OFsimulation", filter="Slice1")[0].getData().compute()
    data = op.arrangeSlice(data=slice1, ydir=False)

ydir=False means that the y component of the velocity is negative.
The data was loaded from the database after it was added as demonstrated in the
10 minutes tutorial.

The function that arrange clips works similarly.

.. code-block:: python

    clip1 = datalayer.Measurements.getDocuments(projectName="Examples",
             type="OFsimulation", filter="Clip1")[0].getData().compute()
    data = op.arrangeClip(data=clip1)

Both functions can be used to save the data to the database.
The procedure is the same for both functions.

.. code-block:: python

    data = op.arrangeClip(data=clip1, save=True,
                          path="/home/ofir/Projects/clip.hdf", key="Clip", addToDB=True)

The path and key are optional. The default path is the current path,
in which case the name of the file is ArrangedClip.hdf or ArrangedSlice.hdf.
The default key is "Clip" or "Slice".
In addition, any additional kwargs may be added, and they will be used
as descriptors for the new document.

Interpolate values in height
----------------------------

The module can be used to intepolate wanted values in a specific height over terrain.
For example, in order to find the velocity in height 20 m above the terrain in a slice,
we can use the next function:

.. code-block:: python

    data = op.makeSliceHeightData(data=Slice1, variable="Velocity", height=20, limit=100)

The limit is the minimal distance between points.
The function returns a dataframe with the two columns, one for the distance downwind
and one for the variable.

A similar function exists for clips:

.. code-block:: python

    data = op.makeClipHeightData(data=Clip1, variable="Velocity", height=20, limit=100)

This function actually uses the one for slices for different fixed x coordinates.
In order to save time, the function doesn't use all possible x coordinates, but skips
with fixed gaps. The default gap is 5 steps, but it can be changed.
Because this procedure takes a long time, the outcome can be saved and add to the database.
For example, for adding a calculation with gap of 2 steps, one can use:

.. code-block:: python

    data = op.makeClipHeightData(data=Clip1, variable="Velocity", height=20, limit=100,
                                skips=2,  path="/home/ofir/Projects/clip20.hdf", key="Clip", addToDB=True)

The function adds the height as a descriptor in the document by default.
The outcome of this function is a dataframe with columns for x coordinates,
y coordinates and requested variable values.

Finding data for height plot
----------------------------

One may wish to plot a variable in a fixed point in different heights.
However, for many points, there is data for limited heights.
The next function can help find detailed points.
It returns a list of points for which there are at least n different heights,
n is a parameter set by the user. The default is 10.

.. code-block:: python

    data = op.findDetaiedLocations(data, n=20)
