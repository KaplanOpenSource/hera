.. _HermesWorkflow:

Hermes-workflow toolkit
========================

Overview
--------

The hermes workflow is used to automate the creation of configuration files and running files
for running general application.

Each workflow is stored as a JSON format. In order to execute a workflow, the JSON is translated  to python program by the
Hermes package. The Hermes-workflow toolkit manages the toolkits in the project. That is, the toolkit allows
the user to add workflows to the project, check if a simulation exists, and compare two simulations.

Each workflow has a name (the name of the workflow json file) and belong
to a group of simulations (simulationGroup). Generally, the names of the simulations
are [group name]_[group id]. However, the user can determine names that do not follow
this convention.

To simplify the use of the toolkit, the toolkit also supplies command line interface (CLI)
that allows the user to perform all the operations.

Workflow JSON structure
------------------------

The  Hermes-workflow JSON has two parts:

..  code-block:: javascript

    {
        "workflow": {....} ,
        "heraMetaData": {...}
    }

1. workflow     - The definition of the workflow.
                  The structure is delineated in the Hermes package documentation.

2. heraMetaData - A description of the metadata that is used to describe.
                  The structure of the heraMetaData JSON is:

..  code-block:: javascript

    {
        "workflowType": "OF_FlowField" ,
        "projectName": "NTA2022",
        "simulationGroup": "simpleStation",
        "caseExecution": {...}
    }

Where:

*  workflowType is the type of the workflow. This parameter will
be used in for more specialized workflows.

Currently we support:

#. OF_FlowField  : A specialized workflow for solving flow fields with OpenFOAM.
#. OF_Dispersion : A specialized workflow for solving lagrangian dispersion

* projectName is the default project name to use

* simulationGroup : The name of the group to add the simulations to.

* caseExecution : Information on how to run the simulation (and instruction of
how to build the runtime).

..  code-block:: javascript

    {
        "caseExecution": {
            "parallelCase": true,
            "runFile": [
                {
                    "name": "blockMesh",
                    "couldRunInParallel": false,
                    "parameters": null
                    "foamJob": true
                },
                .
                .
                .
        }
    }

Where:

* parallelCase : a flag to determine whether to use parallel execution if the node supports it.

* runFile : A list of programs to run.
        each program is determined by:

    * name: The name of the program to run.
    * couldRunInParallel: If true and the case can run in parallel, then
                          use foamJob -parallel
    * parameters: A string of parameters to add to the execution. Currently it is not dynamic.
    * foamJob: if true use foamJob, else just use the string in the name.

Usage (CLI)
-----------

In this section we will describe how to use the toolkit with the CLI.

The CLI allows the user to :

#. Add workflows to the project.
#. List workflows in the project.
#. Build the python execution file.
#. Run the python execution file.
#. Build and Run the python execution file.
#. Compare simulations.

Add workflows to the project
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adding a workflow to the project using the CLI has  ... stages.

#.  Determine the simulation and group names.
    The default behaviour assumes the workflow file name has the format
    [group name]_[group id].

    Then, the default is use the workflow file name as the simulation name,
    and parse it to get the group name and id.

    However, when using the CLI the user can determine the group name
    and can set the simulation name to be of the default format with the
    next available ID in the group.

    Note: If the simulation name is not [group name]_[group id],
          then the group-id of the simulation will be None.

#. Add the simulation to the database.
   If the name exists, or if the workflow already exists in the DB (possibly
   with another name) then it will raise an error.

   If the name of the simualation exists,
   use --overwrite to update the value of the simulation with the given workflow

   If the simulation data already exists in the DB, use --force
   to add it again with the new name.

   The record in the DB has the following fields:

..  code-block:: javascript

    {
        groupName : <group name>,
        groupID : <group ID>,
        simulationName:  <simulationName>,
        workflow    : workflow JSON,
        parameters: <The parameters of all the nodes>
    }

    The resource of the document is the dicrecotry of the simulation, the type is STRING
    and the type is the type of the workflow.

#. Perform addition actions that the user requested (using the


Using the CLI is as follows:

.. code-block::

    >> hera-workflows add <workflow file>
                         [--projectName <projectName>]
                         [--groupName <groupName>]
                         [--overwrite]
                         [--force]
                         [--assignName]
                         [--action Add|AddBuild|AddBuildRun]

Adds the workflow with the name of the workflow file.

* If --groupName appears use the name supplied as the group name.
  Otherwise deduce the groupname from the workflow file name.

* If --overwrite exists than overwite the DB document with the contents
  of the file.

* If --force exists than allow the addition of workflow that exists in the DB under a different name.

* If --assignName exists then find the next available ID in the group and use it.

* Use the --action to add, add and build the python execution or add, build the execution python and
then execute it.

List workflows in the project
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

List all the simulations







