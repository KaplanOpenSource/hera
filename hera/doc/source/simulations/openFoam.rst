OpenFOAM
========

OpenFOAM simulations results can be managed using hera.
The results can be converted to convinient pandas dataframe and be analysed
using built-in plotting functions.

Ten minutes tutorial
--------------------


This tutorial first shows how to load openFoam simulation results to the database.

The process is devided to two parts: the first executes operations on the raw data and saves it to the disk,
and the second loads the data to the database.
The first part is conducted in a python-2 environment, and the second one in python-3.
Both operations use a command line interface.

Building a JSON pipeline file
.............................

In order to execute operations on the data and saving it to the disk, one has to built a json file that
specifies the operation.
The file needs to have the next structure:

.. code-block:: python

    "metadata" : {
               "datadir"    : "/home/ofir/Projects/2020/Carmel5/results",
               "timelist"   : [1000],
               "fields"     : ["U"],
               "MeshRegions": ["internalMesh"]
        },
        "pipelines" : {
                         "Slice1" : {
                             "type" : "Slice",
                             "write"   : "hdf",
                             "params" : [["SliceType.Normal",[1,-0.6,0]], ["SliceType.Origin", [598167.367, 6339602.5, 0]]],
                             "downstream" : {
                                             "Slice2" : {
                                                           "type" : "Slice",
                                                           "write"   : "hdf",
                                                           "params" : [["SliceType.Normal",[1,0.6,0]], ["SliceType.Origin", [598167.367, 6339602.5, 0]]]
                                                   }
                         }
                     }
    }


The json file consists a part called metadata, and a part called pipeline.
The metadata may consist the parameters shown above, that affect the operation:
they controll which fields, timesteps and mesh regions would be used. These parameters are not mandatory. It may also consist additional parameters that would be used as descriptors in the hera database.

Executing operations on the data
................................

The execution is done using the command line interface:

.. code-block:: python

    hera-OpenFOAM executePipeline [JSONpath] [name] [casePath] [caseType] [sourcename]

JSONpath is the path of the json file,
name is a name used for the new files, casePath is the directory of the openFOAM project.
caseType and sourcename are optional.
caseType is either 'Decomposed Case' for parallel cases or 'Reconstructed Case'
for single processor cases, the default is 'Decomposed Case'.
sourcename is a connection string to the paraview server.
The default is None, which work locally.

For example,

.. code-block:: python

    hera-OpenFOAM executePipeline "/home/ofir/Projects/openFoamUsage/askervein/test.json" "test" "/home/ofir/Projects/openFoamUsage/askervein" "Reconstructed Case"


Loading the data to the database
................................

The loading is done using a command line interface:

.. code-block:: python

    hera-OpenFOAM [path] [name] [projectName] --keepHDF

The path is the full directory of the directory
specified as the metadata "datadir" in the json file.
The name is the name used for the executePipeline.
The projectName is used as the project name in the database.
for example,

.. code-block:: python

    hera-OpenFOAM load "Development/Hera/hera/simulations/openfoam/postprocess/dir" "test" "Example"

This command saves the results of each filter in a parquet file.
A document that links to the parquet is added to the database.
The type indicated in the metadata is "OFsimulation".
In addition, a descriptor called "filter" holds the name of the filter,
for example, "Slice2", and a parameter called "pipeline" holds the whole pipeline,
for example, "Slice1_Slice2".

The operation is deleting the hdf files that the vtkpipe.execute function has built.
All the data that was held in the hdf files is now saved in a more organized manner in the parquet files.
If one wishes to keep the hdf files, it can be done by adding "--keepHDF" at the end of the command:

.. code-block:: python

    hera-loadOF load "Development/Hera/hera/simulations/openfoam/postprocess/dir" "test" "Example" --keepHDF

Arranging data
..............

The openfoam module has built-in functions which are very usefull to use.
They may add the velocity magnitude, terrain height, height from terrain
and distance downwind, or interploate values for a specific height.
For example, we will load from the databae and arrange the slice that was
added by the lines above by these lines:

.. code-block:: python

    from hera import datalayer
    from hera import openfoam

    op = openfoam.dataManipulations(projectName="Example")
    slice1 = datalayer.Measurements.getDocuments(projectName="Examples",
             type="OFsimulation", filter="Slice1")[0].getData().compute()
    data = op.arrangeSlice(data=slice1, ydir=False)

ydir=False means that the y component of the velocity is negative.
A full description of the dataManipulations module can be found ahead.

Plotting
........

Plotting can be done using the plotting module.
Examples can be find in hera-gallery. The module is called plotting.


.. code-block:: python

    plotting = openfoam.Plotting()

Usage
-----

.. toctree::
    :maxdepth: 1

    openFoam/ManagingData
    openFoam/analysis