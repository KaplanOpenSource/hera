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
        "projectName": <projectName>
    }

However, this file might contain more data, which depend on which solver you use.

Concepts
--------

Each simulation is represented by a single hermes workflow file that describes how to generate and entire
openfoam case. The structure of the workflow is described in the hermes documentation.
Note that to execute a workflow, the hermes engine builds (automatically) a python file that executes it. Hence,
if you change the hermes file, you must recreate the python program.

Each simulation is a part of a simulation group. This allows the user to compare a certain set of simulations to each other.
The hera code assumes that the name of each simulation is <simulation group>_<id>.json. If the name of the simulation
is not of that structure, then the user must supply the group name.

Listing all the simulation groups
---------------------------------

Listing the simulation groups is performed by the CLI:

    >>  hera-workflows list groups


Adding simulation to simulation groups (with or without execution)
------------------------------------------------------------------

A workflow can be added to the group with or without executing it.

1. Adding a workflow to the DB without executing it

    >> hera-workflows add <workflow name> [--overwrite]

If it exists, then use the --overwrite flag to update the parameters parameter

Note, that in order to be added to a workflow the workflow name must be 'groupName'_'identifier'.
For example the file LargeRoom_1.json is a simulation of the LargeRoom group.

2. Execute workflow and add it to the DB:

  Either:

	>> hera-workflows add <workflow name>	--execute

  Or:

	>> hera-workflows execute <workflow name> [--overwrite]

If it exists, then use the --overwrite flag to update parameters and rerun the flow.


Listing/comparing all the simulations in the  group
---------------------------------------------------



List of supported solvers
-------------------------

The following list is supported for openFOAM V10:

Flow computation
~~~~~~~~~~~~~~~~

- simpleFOAM            : Steady state, incompressible solver.
- indoorFOAMBoussinesq  : Dynamic, incopressible solver. Following Ogura and Phillips (1962).

Dispersion
~~~~~~~~~~~~~~~~

- :ref:`StochasticLagrangianSolver  <openFOAM_StochasticLagrangianSolver>`  : Stochastic lagrangian solver of passive/inertial particles.
- EvaporationStochasticLagrangianSolver                                     : Stochastic lagrangian solver of passive/inertial evaporating particles (with solution) and
                                                                              the evaported cloud.



