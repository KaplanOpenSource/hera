Load Simulation Results
=======================

The next notebook shows how to load openFoam simulation results to the database.

The process is devided to two parts: the first executes operations on the raw data and saves it to the disk,
and the second loads the data to the database.
The first part is conducted in a python-2 environment, and the second one in python-3.

Executing operations on the data
--------------------------------

In order to execute operations on the data and saving it to the disk, one has to built a json file that
specifies the operation.
The file needs to have the next structure:

.. code-block:: python

    "metadata" : {
               "datadir"    : "results",
               "timelist"   : [1000],
               "fields"     : ["U"],
               "MeshRegions": ["internalMesh"]
        },
        "pipeline" : {
                         "type" : "Slice",
                         "guiname" : "Slice1",
                         "write"   : "hdf",
                         "params" : [["SliceType.Normal",[1,-0.6,0]], ["SliceType.Origin", [598167.367, 6339602.5, 0]]],
                         "downstream" : {
                                         "pipeline" : {
                                                           "type" : "Slice",
                                                           "guiname" : "Slice2",
                                                           "write"   : "hdf",
                                                           "params" : [["SliceType.Normal",[1,0.6,0]], ["SliceType.Origin", [598167.367, 6339602.5, 0]]]
                                                       }
                         }
                     }
    }


The json file consists a part called metadata, and a part called pipeline. The metadata may consist the parameters shown above, that affect the operation: they controll which fields, timesteps and mesh regions would be used. These parameters are not mandatory. It may also consist additional parameters that would be used as descriptors in the hera database.

The operation uses the vtkPipeline module:

.. code-block:: python

    import...
    import json
    with open('test.json') as json_file:
         data = json.load(json_file)

    name = "test" # A name used for the files.
    casePath = "/home/ofir/Projects/openFoamUsage/askervein" # The path of the simulation.
    caseType = "Reconstructed Case"
    servername = None

    vtkpipe = VTKpipeline(name=name, pipelineJSON=data, casePath=casePath, caseType=caseType, servername=servername)

The operation is execuded using this line:

.. code-block:: python

    vtkpipe.execute("mainReader")


Loading the data to the database
--------------------------------

The loading is done using a command line interface:

.. code-block:: python

    hera-loadOF [path] [name] [projectName] --keepHDF

The path is the full directory of the directory
specified as the metadata "datadir" in the json file.
The name is the name used for the VTKpipeline method.
The projectName is used as the project name in the database.
for example,

.. code-block:: python

    hera-loadOF "Development/Hera/hera/simulations/openfoam/postprocess/dir" "test" "Example"

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

    hera-loadOF "Development/Hera/hera/simulations/openfoam/postprocess/dir" "test" "Example" --keepHDF