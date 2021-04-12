"""
==============================================
Risk Areas plot
==============================================
"""

#######################
# This example shows how to plot the risk areas projected on a population polygon.
# We will project a cloud of our example agent on a polygon of Katsrin, populated using the demography data.
# Here we will load the buildings data, the demography and LSM simulation directly from the files directory;
# the data may be managed using the appropriate tools.
#
# First, we will load the polygon and demography.

from hera import toolkitHome
import geopandas
from unum.units import *
import matplotlib.pyplot as plt

projectName = "documentation"
Demography  = toolkitHome.getToolkit(projectName=projectName,toolkitName="GIS_Demography")
katsrinBuildings = geopandas.read_file("KatsrinBuildings")
KatsrinCityOnly = Demography.analysis.createNewArea(shapeNameOrData=katsrinBuildings,dataSourceOrData="KatsrinDemography/Katsrin.shp").getData()

#######################
# Then, we calculate the concentration field.

LSM  = toolkitHome.getToolkit(projectName=projectName,toolkitName="LSM")
simulation = LSM.singleSimulation(resource="netcdf")
Concentration = simulation.getConcentration(Q=10*kg)

#######################
# Now, we will get an agent object using the example agent description, and calculate the risk areas.

risk = toolkitHome.getToolkit(toolkitName="RiskAssessment", projectName = projectName)
description = {
    "effectParameters" : {
        "tenbergeCoefficient" : 1.5
    },
    "effects" : {
        "RegularPopulation":{
            "type": "Lognormal10",
            "calculator":{
                "TenBerge" : {"breathingRate":10}
            },
            "parameters":{
                "type": "Lognormal10DoseResponse",
                "levels":["Severe","Light"],
                "parameters":{
                "Severe": {
                    "TL_50" : 10,
                    "sigma": 0.5
                },
                "Light": {
                    "TL_50" : 1,
                    "sigma": 0.5
                }
            }
            }
        },
        "PAC1min":{
            "type":"Threshold",
            "calculator":{
                "MaxConcentration" : {"sampling":"1min"}
            },
            "parameters":{
                "type":"Threshold",
                "levels":["1","2"],
                "parameters":{
                    "1":{"threshold":"1*mg/m**3"},
                    "2":{"threshold":"2*mg/m**3"},
                }
            }
        }
    }
}
Agent = risk.getAgent(description)
riskAreas = Agent.RegularPopulation.calculate(Concentration, "C", isel={"datetime":-1})

#######################
# Finally, we will plot the areas in which we get casualties on top of the city polygon.
# The release point is indicated by a red star

windAngle = 0
x_coordinate = 263500
y_coordinate = 766750
ax, retProj = risk.presentation.plotCasualtiesProjection(results=riskAreas,area=KatsrinCityOnly,loc=[x_coordinate,y_coordinate],
                                             mathematical_angle=windAngle, severityList=["Light","Severe"],
                                                            cycler=plt.cycler(fc=plt.rcParams['axes.prop_cycle'].by_key()['color'])*plt.cycler(ec=['black']))
x,y = KatsrinCityOnly.geometry[0].exterior.xy
ax.plot(x,y,color="black")
ax.scatter(x_coordinate,y_coordinate,marker="*",color="red")