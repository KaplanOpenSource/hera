.. _openFOAM_simpleFoam:

simpleFoam
##########

Introduction
-------------

The simpleFoam toolkit aids in executing steady-state incompressible solvers.


Usage
-----

The stages required to run a simpleFoam with the hera toolkit:

1. Copy the basic hermes workflow from 'hermes/examples/openFOAM/simpleFOAM.
2. Customize the workflow (the instructions here are for indoor workflow, outdoors will be build later):

    1. Write the name of the object file in the Parameters.Execution.input_parameters.objectFile path

      - If you don't need snappyHexMesh, you can remove this, and remove all the refereces to the obj file,
        snappy and  surfaceFeatures exatraction in the workflow.

    2. Update the directory to the directory where the simulations are executed in the fileWriter node.

    3. Update the regions in the snappy node, and update the locationInMesh in snappyDict.

    3. Run the hera-openfoam to create the boundary conditions
       and the blockmesh vertices.

..  code-block::

    >> hera-openfoam objects createVerticesAndBoundary <object file name> --fields [list of field names]

    The list of field names is required to create the boundary conditions for the changeDictionary and node.


    4. Replace the vertices in the blockmesh node.

    5. Change the parameters in the snappyHexMesh. Note that the names of the boundaries were given to you
       in step 2, and that the name of the obj file in the simulation is always building.obj

3. Execute the workflow:
    We can either execute the workflow directly with hermes or with hera support.

    If the simulation is intended to be used for dispersion, there is an advantage to use it with hera,
    since the name of the flow will be used as an original flow field (OFF).

    Executing with the Hera workflow also allows you to avoid running multiple simulations.

    To execute with hermes:

..  code-block::

    >> hermes-workflow buildExecute <flowField name>.

    To execute workflow and add to the db (if not there):

..  code-block::

    >> hera-workflows execute <workflowname> [--projectName <projectName>] [--workflowGroup <workgroup name>] [--overwrite] [--assignName] [--allowDuplicate]

    Where:
    * projectName <projectName>: specifies the project name to use. If not specified, tries to load from the caseConfiguratino file.
    * workflowGroup <workgroup name>: specifies the workgroup name
    * overwrite: overwrite the current values over the name in the database if exists.
    * assignName : finds the next free ID and uses it as the name for the simulation. Useful mainly for automated creation of workflows.
    * allowDuplicate: allows duplicate identical workflows in the database

    Just adding to the DB, without executing:

..  code-block::

    >> hera-workflows add <workflowname> [--projectName <projectName>] [--workflowGroup <workgroup name>] [--overwrite] [--assignName] [--allowDuplicate]

    Where:
    * projectName <projectName>: specifies the project name to use. If not specified, tries to load from the caseConfiguratino file.
    * workflowGroup <workgroup name>: specifies the workgroup name
    * overwrite: overwrite the current values over the name in the database if exists.
    * assignName : finds the next free ID and uses it as the name for the simulation. Useful mainly for automated creation of workflows.
    * allowDuplicate: allows duplicate identical workflows in the database




