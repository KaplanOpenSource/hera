.. _RiskAssessmentPage:

Risk assessment
===============

The risk assessment module provides tools that help making risk assessments.
It helps assessing the effects of dispersions of hazardous chemicals,
and the potential benefits of protection policies.

First Use
---------

Before we use the package for the first time in a project, it is necessary to load the relevant databases.

The two databases that are required to run risk analysis are:

- Demography: The spatial distribution of the population (using the :ref:`Demography toolkit <DemographyPage>`).

- Agents: The toxicity defintions of the different agents.

Existsing databases are available in the hera-data pckage.

Other databases/source that can be useful are the shape definitions (using the :ref:`Shape toolkit <ShapePage>`),
for defining areas of interest and topography (using the :ref:`Topography toolkit <TopographyPage>`)/Buildings
(using the :ref:`Buildings toolkit <BuildingsPage>`) to define dispersion scenarios.

Usage
-----

.. toctree::
    :maxdepth: 3
    :caption: Contents:

    riskAssessment/effectCalculation
    riskAssessment/ProtectionPolicies
    riskAssessment/minimalRiskAssessmentProject
    riskAssessment/ExampleOfRiskAssessment1


API
---

.. autoclass:: hera.riskassessment.riskToolkit.RiskToolkit
    :members:
    :undoc-members:
    :inherited-members:
