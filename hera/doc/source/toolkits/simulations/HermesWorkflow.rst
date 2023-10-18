.. _HermesWorkflow:

Hermes-workflow toolkit
========================

Overview
--------
The Hermes workflow is designed to automate the creation of configuration files and the execution of files for running general applications.

Each workflow is stored in JSON format. To execute a workflow, the JSON is translated into a Python program using the Hermes package.
The Hermes-workflow toolkit manages the toolkits within a hera project.
This toolkit enables users to add workflows to the project, check for the existence of a workflow based on its parameters, and compare different
workflows.

Each workflow has a name, typically based on the name of the JSON workflow file, although this is not mandatory.
Additionally, each workflow belongs to a simulation group (simulationGroup). The toolkit allows the users
to simply compare the parameters of the different workflows within the group.
Generally, the simulations are named [group name]_[group id]. However, users can choose names that do not follow this convention.
The groups are defined dynamically, that is, if there is a workflow with a group defined.

The toolkit can be used as a library from code, or directly from a command-line interface (CLI). This CLI enables users to perform all operations conveniently.

Basic usage
-----------

The toolkit allows the user to:

#. Add a workflow to a group in a project.
#. List the simulation groups in a project.
#. List the names of the workflows in group.
#. Compare the parameters of different simulations.
#. Export the workflow from the database to a file in the directory.
#. Build the python execution file.
#. Run the python execution file.
#. Build and Run the python execution file.

All these actions can be achieved with the CLI.

Preparing usage
---------------

Preparing the usage of the simulations has two steps:

#. Building the workflow template(s)
#. Setting up default configuration options.

Building the workflow template
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Building the workflow template can be achieved by :

#. Building the workflow from scratch.
#. Adjust an existing template manually
#. Adjust an existing template using case configuration.

Building a workflow from scratch and adjusting it manually is documented in the hermes workflow.

Adjusting the existing template using case configurateion is covered here.
We note that sometimes it is necessary to adjust the results manually (or using the GUI in the hermes workflow).

The specialization of the workflow perfomed by setting the values of different parameters in the template
to reflect the needs of the specific simulation. For example, in a template for wind calculation using OpenFOAM,
the specialization is setting up the topography and the blockmesh (and possible the urban region, if needed).

The current [planned] specializations are:

* Flow field - indoors
* Flow field - outdoors (topography and urban)
* Stochastic Lagrangian Dispersion

[The explanation of the specialization for the different cases will be described in the future]

Setting up default configuration options.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The defaults are stored in the **caseConfiguration.json** file.
Specifically, the **projectName** that is used when no project name is supplied to the CLI
is given in:


..  code-block:: javascript

    {
        "projectName" : <project name>
    }

Command Line Interface
----------------------

In this section we will describe how to use the toolkit with the CLI.

List the workflow groups in a project.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

List all the workflow groups in the project.
A workflow groups is defined when a simulation was added to that group.
A group is deleted when all the simulation that belong to the group were
deleted from the project.

.. code-block::

    >> hera-workflows list group <workflow file>
                         [--projectName <projectName>]

List workflows in a group
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Listing all the workflows in the group:

.. code-block::

    >> hera-workflows list group <workflow file>
                         [--projectName <projectName>]
                         [--nodes] or [--parameters]

* --nodes flag will print a list of the nodes for each workflow.
* --parametrs flag will also print the list of parameters.

Add/update workflows
^^^^^^^^^^^^^^^^^^^^^

A workflow that was added to the project, belongs to a workflow group automatically.
The user can supply the workgroup, or let hera determine it from the code.

When executing, the code automatically updates the python-workflow, removes the old dependecy files and executes the workflow.

.. code-block::

    >> hera-workflows add <workflow file>
                         [--projectName <projectName>]
                         [--groupName <groupName>]
                         [--overwrite]
                         [--force]
                         [--assignName]
                         [--execute]

Adds the workflow with the name of the workflow file.

* if --projectName is not supplied, the try to read it from the caseConfiguration.json file.

* If --groupName appears use the name supplied as the group name.

  Otherwise deduce the groupname from the workflow file name.
  That is, we assume that the name of the workflow is <groupname>_<id>.json

* If --overwrite exists than overwite the DB document with the contents
  of the file. This allows the update of the workflow

* If --force exists than allow the addition of workflow that exists in the DB under a different name.

* If --assignName exists then find the next available ID in the group and use it.

* Use the --execute to build and execute the workflow.

Execute and add workflows
^^^^^^^^^^^^^^^^^^^^

The execute commands is similar to the add command, but it executes the workflow.

Remember that you can also execute a workflow using the hermes-workflow interface, and bypass the
hera mechanism with the projects. This could be useful to test a workflow, or as an alternative after
it was added.

When executing, the code automatically updates the python-workflow, removes the old dependecy files and executes the workflow.

The syntax is of the execute command is,
.. code-block::

    >> hera-workflows execute <workflow file>
                         [--projectName <projectName>]
                         [--groupName <groupName>]
                         [--overwrite]
                         [--force]
                         [--assignName]
                         [--execute]


Comparing workflows.
^^^^^^^^^^^^^^^^^^^^

When comparing simulations, the tool lists the differing parameters along with their corresponding values. By default, the simulations are displayed as columns and the parameters are displayed as rows.

.. code-block::

    >> hera-workflows compare <obj1> <obj2> ....
                         [--projectName <projectName>]
                         [--longTable]
                         [--transpose]
                         [--format pandas|json|latex]
                         [--file <outputfileName>]

The input obj can take various forms, such as a simulation name,
a directory path on the disk, a file name on the disk, or a workflow group name.
In the case of a workflow group name, all the simulations within that group will be compared to each other.

* if --projectName is not supplied, the try to read it from the caseConfiguration.json file.

* if --longTable is supplied, then the results are pronted as a long table.
  That is, each parameter (that differs) in each simulation is shown in one line.

* if --transpose is supplied, the the simulations are printed as rows and the parameters are printed as lines.

* The --format prints the comparison in different formats.

* if the --file is supplied, then the output is also printed to a file. If the outputfileName
  does not have extension (i.e it is just the name), the the file name will be appended with

Deleting workflow
^^^^^^^^^^^^^^^^^^
When deleting a workflow from Hera, it's important to note that the deletion process only removes the workflow from the project itself. The files and execution directories associated with the workflow are not automatically deleted, requiring additional action from the user.

When a workflow is deleted from the project, it is exported to a file, and a Python script is generated. This script allows the user to remove all directories associated with the workflow's execution. However, the workflow will not be removed from the project if a file exists in its directory, unless the user explicitly requests overwriting.

It is necessary for the user to manually remove both the workflow file and the execution directories, as these actions need to be performed separately.

To remove the workflow(s) from the project type

.. code-block::

    >> hera-workflows delete <obj1> <obj2> ....
                      [--no-export]
                      [--forceOverwrite]

Where obj<i> can be a simulation name or a workgroup.

* If the --no-export flag is supplied, then the workflow will not be exported to the disk.

* if the --forceOverwrite flag is supplied, then the workflow will be overwrite the currently
 existing workflow  on the disk.

Running this procedure creates a completeRemove.py script that will remove the execution directories.
To remove the execution

.. code-block::

    >> python completeRemove.py

Export workflow
^^^^^^^^^^^^^^^

Exporting workflow saves the workflow in the DB to a file.
If file name is not specified, then the output will be the simulation name

Building/executing a workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Building and running a workflow requires a file on the disk. Hence
this option also include the possibilty to export the file from the DB and then to build it and then
to execute it.

Building and executing a workflow take place similarly to the hemres workflow.


Internals
---------

Hera document structure
^^^^^^^^^^^^^^^^^^^^^^^

The toolkit saves each workflow as a documnet in the project with the following
structure

..  code-block:: javascript

    {
        groupName : <group name>,
        groupID : <group ID>,
        workflowName:  <simulationName>,
        workflow    : workflow JSON,
        parameters: <The parameters of all the nodes>
    }

The resource of the document is the dicrecotry of the simulation, the type is STRING
and the type is the type of the workflow.

Stages in adding a workflow to the project
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adding a workflow to the project using the CLI has  3 stages.

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

#. Perform addition actions that the user requested (using the action flag).




