.. _openFOAMToolkit:

openFOAM Toolkit
================

This toolkit manages the execution and analysis of openFOAM simulations in Hera.

Basically, the openFOAM toolkit is an extension of the `hermes workflow <HermesWorkflow>`, and provides some specialized
operations for the openFOAM solver. Therefore, each simulation is represented by a single hermes workflow file
that describes how to generate and entire openfoam case. The structure of the workflow is described in the hermes documentation.
Note that to execute a workflow, the hermes engine builds (automatically) a python file that executes it.
Hence, if you change the hermes file, you must recreate the python program.

As was defined in the hermes workflow, each simulation is a part of a simulation group. This allows the user to compare a certain set of simulations to each other.
The code assumes that the name of each simulation is <simulation group>_<id>.json. If the name of the simulation
is not of that structure, then the user must supply the group name.

Since each openFOAM solver has a specialized characteristics, the openFOAM toolkit
also inludes specialized toolkits for each solver (these are called toolkit extensions).




Supported solvers
*****************

Currently, we support 2 types of solvers. The first are the flow computation
solvers that compute the air flow and the other is the dispersion model that
solve the dispersion using stochastic lagrangian solver.

Flow computation
~~~~~~~~~~~~~~~~

- :ref:`simpleFOAM <openFOAM_simpleFoam>`            : Steady state, incompressible solver.
- indoorFOAMBoussinesq  : Dynamic, incopressible solver. Following Ogura and Phillips (1962).

Dispersion
~~~~~~~~~~

- :ref:`StochasticLagrangianSolver  <openFOAM_StochasticLagrangianSolver>`  : Stochastic lagrangian solver of passive/inertial particles.
- EvaporationStochasticLagrangianSolver                                     : Stochastic lagrangian solver of passive/inertial evaporating particles (with solution) and
                                                                              the evaported cloud.

API
===



