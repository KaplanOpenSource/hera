.. _openFOAM_StochasticLagrangianSolver:

StochasticLagrangianSolver
##########################

General
*******

Tutorial
********


Dispersion workflow
*******************



Using the StochasticLagrangianSolver
************************************

In order to run a stochastic lagrangian solver, we need to perform the following steps:

1. Create a flow field for dispersion from an existing flow field.
2. Create the dispersion case
3. Run the dispersion simulation.

Create a flow field for dispersion from an existing flow field.
===============================================================

The Dispersion Flow Field (DFF) differs from the original Flow Field (OFF) in several aspects.
Firstly, if the OFF is in a steady-state, the DFF will have two time steps: one for the actual time step used in the simulation,
and another time step that is longer than the expected dispersion time.
This is because the stochastic solver interpolates between adjacent time steps, and setting
them as equal would result in a de-facto steady field. If the OFF is dynamic,
it usually contains time steps used to bootstrap the simulation to avoid the effects of initial conditions.
Therefore, the DFF will include only the time after the initialization, and for simplicity, we set that time step to 0.

The DFF also includes fields that are necessary for the dispersion solver but are not part of the solution itself (e.g., ustar).

Input Parameters
^^^^^^^^^^^^^^^^

To initialize a DFF, you start with a ready OFF and use scripts to create it. The creation of the DFF from the OFF requires the following parameters:

* Flow name:
    * case directory
    * simulation name. If using a simulation name, it must be present in your project.
* Flow dynamics:
    * SteadyState: f the flow is in a steady state, specify the time step to use and the duration of the dispersion. In this case, the time in the dispersion field will vary from 0 to the maximum time.
    * Dynamic    :  In this case, the time in the dispersion field will use the time of the flow simulation. The user specifies the first time to be used (to ignore bootstrapping).

* ustar : An estimation of the ustar (friction velocity) in the domain. Currently, we use a constant value, but it can be changed in a later procedure.
* Hmix  : The height of the mixing layer. For indoor simulations, simply type 1000 or another appropriate value.
* CellHeights:  The distance of each cell from the ground. This field is calculated using the buildDistanceFromWalls flag.

How-to create a dispersion flow field (DFF)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The hera-openfoam interface allows the creation of a flow-dispersion case
based on a previously computed OpenFOAM flow, referred to as the original flow.
It is preferable to have the original flow in the project to keep track of its parameters, but it is not mandatory.

To handle the numerous parameters involved, a JSON configuration file is utilized to specify
the dispersion flow field based on the original flow field. The configuration file follows a specific
structure (refer to the documentation on how to create an empty configuration file).

..  code-block:: javascript

    {
        "projectName": <projectName>,
        "Flows": [
                {
                    "originalFlow" : {
                        "source" : <name>,
                        "time" : {
                            type : "steadyState|dynamic",
                            "timestep" : <time>
                        },
                        linkMeshSymbolically : True
                    },
                    dispersionDuration : <duration>
                    dispersionFields : {
                            <see below>
                    }
                },
                .
                .
                .
        ]
    }

* The original flow name (stated in originalFlow.source) can be a workflow name or a directory.

    * The time determines which time steps will be used for the dispersion.
      If the solution is parallel, we use the parallel solution, otherwise,
      we use the single core solution.

        - Steady-state: Use the timestep as time 0 in the dispersion.
                        This timestep is copied to the final dispersion time (and so get a de-facto steady flow).
       -  dynamic: Use all the timesteps from timestep and on. The timestep is time 0 in the
                   dispersion, and that time is subtracted from all the other timesteps.

    * linkMeshSymbolically determines whether the mesh will be copied or just linked symbolically. Symbolic linking
      saves disk space, but makes the present simulation depend on the existence of the original flow.
      the default is True

* dispersionDuration: determines the duration of the dispersion simulation.

* The dispersionFields key determines the fields that will be added to the dispersion flow.
  A field is defined by its dimensions, components (1 for scalar, 3 for vector and 9 for tensor),

  and the values of the boundary fields.

  It is possible to select a predefined field or define the field. The boundary conditions should be stated for either.
  We note, that the boundaries that were not stated are added automatically with the boundary condition zeroGradient.

  When using a predefined field, it is only necessary to state the flow type (compressible, incompressible, dispersion).
  This is becuase sometimes the dimension of the fielding depend on the context (for example pressure has different
  units for copressible and incompressible flows).
  For predefined fields the structure is:

..  code-block:: javascript

                            <FieldName> : {
                                "flowType"   : "compressible|incompressible|dispersion"
                                "boundaryFields" : {
                                        <boundary name> : {
                                                <property 1> : <property value>,
                                                .
                                                .
                                                <property n> : <property value>
                                        },
                                        .
                                        .
                                }
                                "internalField" : <value>|string
                            }

For defined fields it is necessary to define their units and the name of the each component (for example, for velocity it is usually
Ux,Uy and Uz. For scalars it is null.

..  code-block:: javascript

                            <FieldName> : {
                                "dimensions" : {kg=<int>,m=<int>,s=<int>,K=<int>,mol=<int>,A=<int>,cd=<int>}
                                "componentNames" : [<name X>, <name Y>, <name Z>],
                                "boundaryFields" : {
                                        <boundary name> : {
                                                <property 1> : <property value>,
                                                .
                                                .
                                                <property n> : <property value>
                                        },
                                        .
                                        .
                                },
                                "internalField"  : <value>|<list>|string
                            }


* The boundaryField lists the values of the boundaries.
    The struct translates the value to the openfoam dict.
    For example the following translates to drichlet boundary condition.

..  code-block:: javascript

                                        "east" : {
                                                "type" : "fixedValue",
                                                "value" : "uniform 0"
                                        },
                                        .
                                        .
                                }
                            }

* The internalField can be a value (constant for scalar, vector and tensor), list (constant for vector/tensor) or a string
  that can include a parquet file to be red. The structure of the parquet should be similar to a parquet that was loaded
  with the load method in OF objects (see ...).

To create the dispersion field with CLI use
.. code-block::

    >> hera-openfoam stochasticLagrangian createDispersionFlow <configuration file>
                                                                     [--projectName <projectName>]


If <configuration file> is not stated, try to use the caseConfiguration.json file.

**Using the toolkit**

First, lets import the toolkitHome

.. code-block::

    from hera import toolkitHome
    projectName = "my-project"

Then, initializa a SIMULATIONS_OPENFOAM toolkit:

.. code-block::

    dispersionToolkit = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)



Create the dispersion case
==========================

Creating a dispersion case includes creating the directories of the case, linking to a dispersion flow,
creating the configuration files required to run the StochasticLagrangianSolver, and adding the workflow
to the DB.

.. code-block::

    >> hera-openfoam stochasticLagrangian createDispersionFlow <configuration file>
                                                                     [--projectName <projectName>]





File examples
**************



Notes for developers
*********************

The procedures for all the stages are implemented in the StochasticLagrangian module of the openfoam toolkit.
See reference below.

Steps in creating the dispersion flow field (DFF)
=================================================

In order to allow maximal flexibility, it is possible to create the DFF
without accessing to the hera mongodb. This is sometimes required (especially when using hermes)
for simpler settings. We note that there is no DB, then only the case directory of the  original flow field can be supplied.

The steps in the creation of DFF will take the following steps:

1. Check if a DFF with the requested parameters is already in the project.
   If the database is not enabled, the treat as if the flow does not exist.

2. If flow does not exist, (or exists and it is to be overwritten) create it.
   The name of the flow is <case name>_Dispersion_<id>.
   id is the first ID that is available in the directory.
3. Create the case:
    3.1 Copy the system and constant from the original flow.
        If the original flow is parallel, the this processes is repeated for all the processor* sub-dirs.
    3.2 For each timestep:
        3.1.1 Copy the time step from the original. Map the time step in the original to be [timestep-starting time]
              in the dispersion simulation. If the original flow is parallel, the this processes is repeated
              for all the processor* sub-dirs.
        3.1.2 Link (or copy) the mesh. If the original flow is parallel, the this processes is repeated
              for all the processor* sub-dirs.
    3.3 create the new fields with their values in each time step.
        Take the boundary conditions from the existing fields.
    3.4 If the original is parallel, create empty directories with the timesteps to overcome a
        a bug in the stochastic solver that recognizes only the time steps in the main directory
        and not in the parallel case.
4. Add The new DFF to the database.

Steps in creating dispersion case and link a dispersion flow field to it
========================================================================

1. Create a the dispersion directory with the system and control dict.
2. Copy from the dispersion flow field the system and control.



API
===

.. autoclass:: hera.simulations.openFoam.StochasticLagrangian.stochasticLagrangianDataLayer
    :members:
    :undoc-members:
    :inherited-members:
