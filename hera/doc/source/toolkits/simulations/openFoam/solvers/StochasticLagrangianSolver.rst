.. _openFOAM_StochasticLagrangianSolver:

StochasticLagrangianSolver
##########################

The Stochastic lagrangian solver extends the stochastic solver of OpenFOAM by adding a stochastic parameterization
of the turbulence to the lagrangian solver. The solver uses the mean flow field that was prevsiouly computed to calculate
the dispersion of lagrangian particles.

To do so, the solver interpolates the values of the mean fields both temporally and spatially (i.e to the location of the
particle). Because the stochastic parametrization requires additional mean fields that are not solved by the
navier solver (for example :math:`u_*`). However, adding those fields to the original solution of the flow field
can be limiting when it is necessary to examine the effect of the parametrization on the dispersion.

To solve this problem, using the StochasticLagrangianSolver includes an intermediate step that creates a
copy of the original flow field (OFF), but adds the required parametrization fields to it. To save space, the
dispersion flow field (DFF) does not copy the mesh and the results, but rather, creates a symbolic links to them.

Once the DFF was created, it is possible to create the dispersion case. Since the StochasticLagrangianSolver required
the configuration in numerous files (as is usual for openFOAM solvers), we use Hermes workflows to automate this procedure
and convert JSON file to the configurations files.

Comparing dispersion workfows and other management can be achieved through the hera-workflows CLI.

10min tutorial
**************

For a 10min tutorial, copy the example that is located in ... .

1. Create the original flow field (OFF) by running the hermes workflow of a flow field.
2. Create the dispersion flow field (DFF) by describing the flow in the caseConfiguration
   and executing

    >> hera-openfoam stochasticLagrangian dispersionFlow create <OriginalFlowField> [--DFF <dispersion flow field names>]

Where the <OriginalFlowField> is the name of the flow, the hermes workflow file or the directory.
The batch file will create one DFF for each type that is defined in the caseConfiguration file (see below).

3. Create the dispersion case by running the hermes workflow of the dispersion.

    >> hera-openfoam ...


Create a flow field for dispersion
**********************************

In this section we describe how to create a Dispersion Flow Field (DFF) from an original Flow Field (OFF).
The DFF differs from the OFF)in several aspects.
Firstly, if the OFF is in a steady-state, the DFF will have two time steps: one for the actual time step used in the simulation,
and another time step that is longer than the expected dispersion time.
This is because the stochastic solver interpolates between adjacent time steps, and setting
them as equal would result in a de-facto steady field. If the OFF is dynamic,
it usually contains time steps used to bootstrap the simulation to avoid the effects of initial conditions.
Therefore, the DFF will include only the time after the initialization, and for simplicity, we set that time step to 0.

The DFF also includes fields that are necessary for the dispersion solver but are not part of the solution itself (e.g., ustar).

Input Parameters
^^^^^^^^^^^^^^^^

The following parameters are required to create a DFF from an existing OFF:

* Flow name:
    * case directory
    * simulation name. If using a simulation name, it must be present in your project.
* Flow dynamics:
    * SteadyState: f the flow is in a steady state, specify the time step to use and the duration of the dispersion. In this case, the time in the dispersion field will vary from 0 to the maximum time.
    * Dynamic    :  In this case, the time in the dispersion field will use the time of the flow simulation. The user specifies the first time to be used (to ignore bootstrapping).

Additional fields are often necessary, depend on the parametrization chosen.
For the Neutral2018, and Indoor2018 parameterizations the following fields are required:

* ustar : An estimation of the ustar (friction velocity) in the domain. Currently, we use a constant value, but it can be changed in a later procedure.
* Hmix  : The height of the mixing layer. For indoor simulations, simply type 1000 or another appropriate value.
* CellHeights:  The distance of each cell from the ground. This field is calculated using the buildDistanceFromWalls flag.


Defining the dispersion flow field (DFF)
========================================

The creation of the DFF requires as input: (1) the name (or location) of the OFF, (2)
the definition the the time step(s) to use for the dispersion and (3) the definition of the fields that will be added to the
OFF.

Currently, the creation of the DFF is supported from the command line or by calling the createDispersionFlowField in the SIMULATIONS_OPENFOAM
toolkit (see below).

Calling from the CLI
~~~~~~~~~~~~~~~~~~~~

The name of the OFF is supplied in the command line by the

    >> hera-openfoam stochasticLagrangian dispersionFlow create <OriginalFlowField> [--DFF <dispersion flow field names>]

Where the <OriginalFlowField> is the name of the flow, the hermes workflow file or the directory.
The batch file will create one DFF for each type that is defined in the caseConfiguration file (see below).

If the DFF is no specified, then the CLI will create all the defined DFFs in the file.


Calling the toolkit directly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to call the toolkit function directly, use

.. code-block:: python

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)
    tk.stochasticLagrangian.createDispersionFlowField(flowName=flowName,
                                                      flowData=flowdata,
                                                      OriginalFlowField=OriginalFlowField,
                                                      overwrite=<true/false>)

Where the flowName is the name of the new flow, flowdata is the JSON that describes the flow (see below),
OriginalFlowField is the name, directory or workflow file of the original flow and overwrite specifies
whether or not it is will be overwritten if it exists.

Definition of the DFF
~~~~~~~~~~~~~~~~~~~~~

The DFF is defined in a JSON file with the follwing structure:

.. literalinclude:: ./example_caseConfiguration.json
   :language: json
   :linenos:
   :caption: case configuration structure

Where:
- The 'directory' key denotes the location that the dispersion workflows are written to.

- <name> is the name of the DFF. The final name of the DFF will be <originalFlow>_DFF_<DFF name>.
- originalFlow: defines how to get the times.
    - time.type : steadyState/dynamic
    - time.timestep : if steadyState, the time to use as start and end (since the flow is constant).
                      if dynamic, the timestep denotes the time that will be mapped to 0.
    - linkMeshSymbolically : If true, links the mesh to the OFF mesh. Otherwise, just copies it.
- dispersionDuration : The last time step to use. If it is steady-state, the first time step will be
                       copied to this time step.


- The dispersionFields key determines the fields that will be added to the dispersion flow.
  A field is defined by its dimensions, components (1 for scalar, 3 for vector and 9 for tensor),
  and the values of the boundary fields.

  It is possible to select a predefined field or define the field. The boundary conditions should be stated for either.
  We note, that the boundaries that were not stated are added automatically with the boundary condition zeroGradient.

  When using a predefined field, it is only necessary to state the flow type (compressible, incompressible, dispersion).
  This is becuase sometimes the dimension of the fielding depend on the context (for example pressure has different
  units for copressible and incompressible flows).
  For predefined fields the structure is:

.. literalinclude:: ./example_field.json
   :language: json
   :linenos:
   :caption: Example of a field

The list of predifined fields is:

.. csv-table:: Predefined fields
   :file: predefinedFields.csv
   :widths: 30, 70
   :header-rows: 1


For fields that are not predefined, it is necessary to define their units and the name of the each component (for example, for velocity it is usually
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


- The boundaryField lists the values of the boundaries.
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

- The internalField can be a value (constant for scalar, vector and tensor), list (constant for vector/tensor) or a string
  that can include a parquet file to be red. The structure of the parquet should be similar to a parquet that was loaded
  with the load method in OF objects (see ...).

Example
~~~~~~~

Here is an example of a full caseConfiguration:

.. literalinclude:: ./caseConfiguration.json
   :language: json
   :linenos:
   :caption: caseConfiguration

Create the dispersion case
***************************

Creating the dispersion case involves in creating the soft links to the mesh directories
of the DFF in the case directory and in the processor<x> (if exists).

Creating the case and the links is usually part of the hermes dispersion workflow, but it can also be
executed manually.



To create the dispersion director manually,

>> hera-openfoam stochasticLagrangian dispersion create <dispersion case name> <DFF name> [--overwrite]

use overwrite to recreate the directory.


Post processing
***************



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


Using Hera directrly (without CLI) [not written]
================================================

First, lets import the toolkitHome

.. code-block::

    from hera import toolkitHome
    projectName = "my-project"

Then, initializa a SIMULATIONS_OPENFOAM toolkit:

.. code-block::

    dispersionToolkit = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM,
                                                projectName=projectName)
