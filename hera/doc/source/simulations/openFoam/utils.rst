Utils
=====

The utils are used to
extract data from OpenFOAM files or create new OpenFOAM files.

Height from ground
------------------
It may be usefull to create a field with the height from ground at each point.
This field is mandatory in order to run an OpenFOAM LSM simulation, for example.
A file which holds this field may be made using the next function.

It should be noted that before using the function,
the cell centers must be extracted to a file.
This can be done by running the next command at the case path:

.. code-block:: python

    cellCenters

Then, the function may be called.

.. code-block:: python

    from hera.simulations.openfoam import makeCellHeights
    casePath = "the case path"
    times = [0,10000]
    makeCellHeights(casePath,times, ground="ground",fileName="cellHeights",fillna=True,fillVal=0)

casePath is the full path of the case.
times is a list of time steps directories in which the file is saved.
ground is the name of the patch from which the vertical distance is calculated;
its default value is "ground".
fileName is the name of the new file;
its default value is "cellHeights".

The height of each cell is calculated using an interpulation of the z coordinate of the
ground patch at the cell's x and y coordinates.
However, there may be points for which these coordinates are outside
the spectrum of the ground patch.
For these points, an interpulation cannot be made.
By default, the nan values found for the z coordinate for these points
is filled with the value 0.
One may deactivate the fill by setting the parameter fillna to False,
and change the value used for the filling using the parameter fillVal.

The file which is written is a list of vectors for all cells, sorted by the
cells order in the mesh.
The first component of each vector is the cell's x coordinate,
the second is the cell's y coordinate,
and the third is its "height",
the vertical distance from the "ground" patch.