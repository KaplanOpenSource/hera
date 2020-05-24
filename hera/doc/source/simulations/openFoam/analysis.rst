Analysis
========

This page shows some basic analysis processes that hera supports.

Tests
-----
The tests module may support performing tests of the validity of the simulation.
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