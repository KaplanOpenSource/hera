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
import pandas
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
# Now, we will plot the areas in which we get casualties on top of the city polygon.
# The release point is indicated by a red star.
# The boundaries of the areas of the cloud for each severity may be plotted, using the parameter plumSeverity.
# We will plot the boundaries of the area in which the toxic load may cause light injuries.
#
# In the next examples we work with mathematical angles, and therefore ue the parameter mathematical_angle.
# If we would want to work with meteorological angles, the windAngle value would have to be assigned to the parameter
# meteorological_angle.

windAngle = 0
x_coordinate = 263500
y_coordinate = 766750
ax, retProj = risk.presentation.plotCasualtiesProjection(
    results=riskAreas,area=KatsrinCityOnly,loc=[x_coordinate,y_coordinate],
    mathematical_angle=windAngle, severityList=["Light","Severe"], plumSeverity=["Light"],
    cycler=plt.cycler(fc=plt.rcParams['axes.prop_cycle'].by_key()['color'])*plt.cycler(ec=['black']))
x,y = KatsrinCityOnly.geometry[0].exterior.xy
ax.plot(x,y,color="black")
ax.scatter(x_coordinate,y_coordinate,marker="*",color="red")

#######################
# The retProj object is a geopandas dataframe that holds numbers of injuries for different polygons.
# It may be used to make a table of the results.
# For example, we will use it to add a table at of the number of injuries at the top of the plot.
# We will sum the number of injuries over all polygons.
#
# Additionaly, several wind directions may be presented on the same plot, by delivering the same axis to the ax parameter.

windAngles = [30,-30]
fig, ax = plt.subplots()
calculatedData = []
for windAngle in windAngles:
    ax, retProj = risk.presentation.plotCasualtiesProjection(
        results=riskAreas,area=KatsrinCityOnly,loc=[x_coordinate,y_coordinate],
        mathematical_angle=windAngle, severityList=["Light","Severe"], ax=ax,
        cycler=plt.cycler(fc=plt.rcParams['axes.prop_cycle'].by_key()['color'])*plt.cycler(ec=['black']))
    retProj["windDirection"] = windAngle
    calculatedData.append(retProj)
x,y = KatsrinCityOnly.geometry[0].exterior.xy
ax.plot(x,y,color="black")
ax.scatter(x_coordinate,y_coordinate,marker="*",color="red")
tableData = pandas.concat(calculatedData).groupby(["severity","windDirection"])["effectedtotal_pop"].sum().reset_index()
pandas.plotting.table(ax,tableData,loc="top")
plt.subplots_adjust(top=0.8)

#######################
# The plot results may also be plotted over an image.
# Here we will load an image and plot is unsing the raster tool.

toolkitName = "GIS_Raster"
raster  = toolkitHome.getToolkit(projectName=projectName,toolkitName=toolkitName)

location = "Katsrin"
extents = {"minX":259600, "minY":762000, "maxX":269600, "maxY":772000}
loadedData = raster.loadExperiment(fileNameOrData="Katsrin.png", extents = extents, saveMode="DB_overwrite", regionName=location, additionalData=dict(units="ITM"))


fig, ax = plt.subplots()
for windAngle in windAngles:
    ax, retProj = risk.presentation.plotCasualtiesProjection(
        results=riskAreas,area=KatsrinCityOnly,loc=[x_coordinate,y_coordinate],
        mathematical_angle=windAngle, severityList=["Light","Severe"], ax=ax,
        cycler=plt.cycler(fc=plt.rcParams['axes.prop_cycle'].by_key()['color'])*plt.cycler(ec=['black']))
    retProj["windDirection"] = windAngle
    calculatedData.append(retProj)
raster.presentation.plot(imageNameOrData=location,ax=ax)
ax.plot(x,y,color="black")
ax.scatter(x_coordinate,y_coordinate,marker="*",color="red")
#######################
# Finally, we may save the figure. In our example we will use the used parameters in the name.

figName = f"projectedCasualties_releasePoint{x_coordinate}_{y_coordinate}_windDir_{windAngle}.png"
# plt.savefig(figName)
