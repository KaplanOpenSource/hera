.. _HermesWorkflow:

Hermes-workflow toolkit
========================

Overview
--------

The hermes workflow is used to automate the creation of configuration files and running files
for running general application.

Each workflow is stored as a JSON format, and then transformed to python program for execution.
The Hermes-workflow toolkit manages the toolkits in the project. That is, the toolkit allows
the user to add workflows to the project, check if a simulation exists, and compare two simulations.

To simplify the use of the toolkit, the toolkit also supplies command line interface (CLI)
that allows the user to perform all the operations.

General functionality
----------------------

The structure of Hermes-workflow JSON has two parts:

1. workflow     - The workflow (see the Hermes package for documentation).
2. heraMetaData - A description of the metadata that is used to describe.

The structure of the heraMetaData,

.. jsonschema::

    {
        "workflowType": "OF_FlowField" ,
        "projectName": "NTA2022",
        "simulationGroup": "simpleStation",
        "caseExecution": {...}



    }





Usage (CLI)
-----------

In this section we will describe how to use the toolkit with the CLI.

The CLI allows the user to :

1. Add workflows


