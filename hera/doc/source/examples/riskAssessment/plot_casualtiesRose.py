"""
==============================================
Casualties rose plot
==============================================
"""

#######################
# This example shows how to plot a casualties rose.
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
        }
    }
}
Agent = risk.getAgent(description)
riskAreas = Agent.RegularPopulation.calculate(Concentration, "C", isel={"datetime":-1})

#######################
# Finally, we will plot the areas in which we get casualties on top of the city polygon.
# The release point is indicated by a red star

windAngles = [0,90,180,270]
x_coordinate = 264500
y_coordinate = 766750
ax, retProj = risk.presentation.plotCasualtiesRose(results=riskAreas,area=KatsrinCityOnly,loc=[x_coordinate,y_coordinate],
                                                   mathematical_angles=windAngles, severityList=["Light","Severe"])

#######################
# One may create different subplots using the ax parameter. This parameter is a list of integers.
# The first value defines the number of rows, the second the number of column, the third the index of the subplot.

ax, retProj = risk.presentation.plotCasualtiesRose(results=riskAreas,area=KatsrinCityOnly,loc=[x_coordinate,y_coordinate],
                                                   mathematical_angles=windAngles, severityList=["Light"],ax=[1,2,1])

ax2, retProj = risk.presentation.plotCasualtiesRose(results=riskAreas,area=KatsrinCityOnly,loc=[x_coordinate,y_coordinate],
                                                   mathematical_angles=windAngles, severityList=["Severe"],ax=[1,2,2])
