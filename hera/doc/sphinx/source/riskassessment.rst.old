Risk Assessment
===============

The risk assessment module provides tools that help making risk assessments.
It helps assessing the effects of dispersions of hazardous chemicals,
and the potential benefits of protection policies.


10-minute tutorial
------------------

This tutorial demonstrates the basic applications of the method.

Getting the concentration field
*******************************
The first step for making a new risk assessment is creating a concentration field.
This can be done using hera LSM, for example (see LSM tutorial).

.. code-block:: python

    sims = LSM.getSimulations(projectName='example')
    chosen_sim = sims[0]
    concentration = chosen_sim.getConcentration()

The agent object
****************

Then, one should usualy create an instance of the agent object,
which holds the properties of the requested agent.

.. code-block:: python

    from hera import riskassessment as risk
    AgentName = "ExampleAgent"
    Agent = risk.Agent(AgentName)

Each agent may have several injury types.
For example, The example agent has the injury type "RegularPopulation", which has two
levels of injury, light and severe.
The effected risk areas for the chosen LSM simulation can be calculated this way:

.. code-block:: python

    injurytype = Agent.RegularPopulation
    riskAreas = injurytype.calculate(concentration, "C")

Population
**********

The number of injuries may be assessed by projecting the risk areas on population data.
The population data may be loaded using hera GIS (see GIS tutorial).

.. code-block:: python

    from hera import GIS
    population = GIS.population(projectName="Examples")
    populationAreas = population.projectPolygonOnPopulation(Geometry="TelAviv")
    injuryareas = riskAreas.project(populationAreas, releasePoint=[x_coordinate,y_coordinate],mathematical_angle=some_angle)

Protection Policies
*******************

Adding protection policies should be done upon the concentraion field.
For example, if we would like to use wearing masks and getting indoors for 10 minutes,
it can be done this way:

.. code-block:: python

    from unum.units import *
    indoor = dict(name="indoor",params=dict(turnover=15*min,enter="30s",stay="10min"))
    masks = dict(name="masks",params=dict(protectionFactor=10,wear="0s",duration="10min"))
    newConcentraion = risk.ProtectionPolicy(actionList=[indoor,masks]).compute(concentraion).squeeze().compute()

Usage
-----

.. toctree::
    :maxdepth: 1

    riskassessment/agent
    riskassessment/ProtectionPolicies