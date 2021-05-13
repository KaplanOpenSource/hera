Tests
=====

The tests tool may support performing tests of the validity of the simulation.
Note that it should be executed in python 2 environment.

.. code-block:: python

    from hera import openfoam
    path = "/home/ofir/Projects/2020/Carmel5" #The path of the case.
    tests = openfoam.tests(casePath=path)

Get a Z-Axis Slice
..................
It is usefull to examine the velocity at height of about 90% of the total height.
One may use the next function:

.. code-block:: python

    data = tests.getHeightSlice()

The function returns a pandas dataframe of a slice perpendicular to the z axis, with the
velocity's components and total magnitude.
The default height is 90% of the total height, but may be altered.
Additional fields may be requested too; however, the velocity should remain. for example:

.. code-block:: python

    data = tests.getHeightSlice(percentage=80, fields=["U","p"])

It is convinient to examine the data using a contour:

.. code-block:: python

    import matplotlib.pyplot as plt
    plt.tricontourf(data.x, data.y, data.Velocity)

In addition, the relative difference between the maximum and minimum velocities
may indicate how valid is the simulation. The difference in percentage can be calculated
like this:

.. code-block:: python

    difference = (data["Velocity"].max()-data["Velocity"].min())/data["Velocity"].min()*100

Get The Difference Of Velocity With Height
..........................................

It might be usefull to find whether there are coordinates for which the velocity decreases
with height. The physical behavior we expect is an increase of velocity with height,
and a decrease may indicate errors.
There ia a function that returns a dataframe with the difference of the velocity in height.

.. code-block:: python

    data = tests.changeOfVelocityInHeight()

Then, we can check whether there are coordinates for which the velocity decreases:

.. code-block:: python

    decreasingVelocityPoints = data.loc[data.diffs<0]

Perform all tests
.................

There is also a function that activate both functions discussed above.
It returns a dictionary with two keys. The first key is called UalongXYdifference.
It is the size of the difference between the minimum amd maximum values of U
at a chosen percentage height.
The second key is called UalongZconsistent. Its value is True if there isn't an area in
which the velocity decreases with height, and False otherwise.

.. code-block:: python

    tests.performAllTests(percentage=80,fields = "U")

