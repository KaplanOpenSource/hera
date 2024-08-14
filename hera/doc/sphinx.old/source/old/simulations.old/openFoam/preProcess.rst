Pre-process
===========

The pre-process tool is used in order to load openFoam simulation results to
the database.

This process is divided to two parts: the first executes operations on the raw data and saves it to the disk,
the second loads the data to the database.
The first part is conducted in a python-2 environment, and the second one in python-3.
To switch to a Python 2 environment and export the appropriate folders you need to check how exactly it works on your
computer, but it should be like:

.. code-block:: python

    conda deactivate
    conda deactivate


    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/raid/opt/OpenFOAM/ThirdParty-7/build/linux64Gcc/ParaView-5.4.1/lib

    export PYTHONPATH=/home/davidg/.local/lib/python2.7/site-packages
                     :/home/davidg/Development/pyhera
                     :/raid/opt/OpenFOAM/ThirdParty-7/platforms/linux64Gcc/ParaView-5.4.1/lib/paraview-5.4/site-packages
                     :/raid/opt/OpenFOAM/ThirdParty-7/platforms/linux64Gcc/ParaView-5.4.1/lib/paraview-5.4/site-packages/vtk

    source $HOME/OpenFOAM/OpenFOAM-7/etc/bashrc WM_LABEL_SIZE=64 FOAMY_HEX_MESH=yes

Both operations use a command Line Interface (CLI).

Building a JSON pipeline file
.............................

In order to execute operations on the data and saving it to the disk, one has to built a json file that
specifies the operation.
The file needs to have the next structure:

.. code-block:: python

    {"metadata" : {
               "datadir"    : "/home/ofir/Projects/2020/Carmel5/results",
               "timelist"   : [1000],
               "fields"     : ["U"],
               "MeshRegions": ["internalMesh"],
               "CaseDirectory": "/home/ofir/Projects/2020/Carmel5"
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
they control which fields, time steps and mesh regions would be used. These parameters are not mandatory.
It may also consist additional parameters that would be used as descriptors in the hera database.
The only mandatory parameter is "CaseDirectory", which is the case directory of the simulation.

The pipeline part may also consist a 'write' parameter. If write exists then the filter will be saved
to a file in the specified format ("hdf"/"netcdf")

Executing operations on the OpenFOAM data
.........................................

The execution is done using the CLI:

.. code-block:: python

    hera-OpenFOAM executePipeline [JSONpath] [name] [casePath] [caseType] [servername] [tsBlockNum]

JSONpath is the path of the json pipeline file,
name is the name used for the target resulted folder and for the new files (if hdf format selected),
casePath is the path to the directory of the openFOAM project.
caseType, servername  and tsBlockNum are optional.
The CLI is positional order depended, so if one is defined- all the optional before must be defined as well
caseType is either 'Decomposed Case' for parallel cases or 'Reconstructed Case'
for single processor cases, the default is 'Decomposed Case'.
servername is a connection string to the paraview server.
The default is None, which work locally.
tsBlockNum is the number of time steps will be saved to single file. The default is 100

For example,

.. code-block:: python

    hera-OpenFOAM executePipeline "/home/ofir/Projects/openFoamUsage/askervein/test.json" "test" "/home/ofir/Projects/openFoamUsage/askervein" "Reconstructed Case"

The resulted files is now saved (in the format specificated in the 'write' property at the pipeline)
in the 'name' folder

Loading the data to the database
................................

The loading is done using a CLI:

.. code-block:: python

    hera-OpenFOAM load [JSONpath] [path] [name] [projectName] -keepHDF

JSONpath is the path of the json file
The path is the full directory of the directory
specified as the metadata "datadir" in the json file.
The name is the name used for the executePipeline.
The projectName is used as the project name in the database.
for example,

.. code-block:: python

    hera-OpenFOAM load "/home/ofir/Projects/openFoamUsage/askervein/test.json" "Development/Hera/hera/simulations/openfoam/postprocess/dir" "test" "Example"

This command saves the results of each filter that has a 'write' property in the pipeline to the database.
hdf files will be converted to parquet file before loading.
A document that links to the filter data is added to the database.
The type indicated in the metadata is "OFsimulation".
In addition, a descriptor called "filter" holds the name of the filter,
for example, "Slice2", and a parameter called "filterpipeline" holds the whole pipeline,
for example, "Slice1.Slice2".

The operation is deleting the hdf files that the vtkpipe.execute function has built.
All the data that was held in the hdf files is now saved in a more organized manner in the parquet files.
If one wishes to keep the hdf files, it can be done by adding "-keepHDF" at the end of the command:

.. code-block:: python

    hera-loadOF load "/home/ofir/Projects/openFoamUsage/askervein/test.json" "Development/Hera/hera/simulations/openfoam/postprocess/dir" "test" "Example" -keepHDF
