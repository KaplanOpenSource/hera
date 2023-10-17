.. _openFOAMToolkit:


openFOAM Toolkit
================

Overview
--------

This toolkit manages the execution and analysis of openFOAM simulations in Hera.

To do so, each openFOAM solver has a dedicated toolkit that allows the users to
setup (using the hermes workflow) and analyze the simulation results.

Generally, the configuration for the case is stored in the file 'caseConfiguration.json'
in the directory that stores all the workflows (the case directories of the simulations are
usually the subdirectories of this directory.

The most basic caseConfiguration.json holds the project name of the simulations:
..  code-block:: javascript

    {
        "projectName": <projectName>,
    }

However, this file might contain more data, which depend on which solver you use.

Concepts
---------

Each simulation is represented by a single hermes workflow file that describes how to generate and entire
openfoam case. The structure of the workflow is described in the hermes documentation.
Note that to execute a workflow, the hermes engine builds (automatically) a python file that executes it. Hence,
if you change the hermes file, you must recreate the python program.

Each simulation is a part of a simulation group. This allows the user to compare a certain set of simulations to each other.
The hera code assumes that the name of each simulation is <simulation group>_<id>.json. If the name of the simulation
is not of that structure, then the user must supply the group name.


List of supported solvers
-------------------------


The current solvers that are supported are:

- simpleFOAM
- :ref:`StochasticLagrangianSolver  <openFOAM_StochasticLagrangianSolver>`


