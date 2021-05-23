"""
==============================================
Casualties bar plot
==============================================
"""

#######################
# This example shows how to plot a casualties barplot.
# This is very usefull in order to compare the influence of different parameters on the results.
# For example, we will project a cloud of our example agent on a polygon of Katsrin, populated using the demography data.
# We will compare between different masses of agent dispersed, and different location of the populations.
#
# Here we will load the buildings data, the demography and LSM simulation directly from the files directory;
# the data may be managed using the appropriate tools.
#
# First, we will load the polygon and demography.

from hera import toolkitHome
import geopandas
from unum.units import *
import pandas
import seaborn as sns
import matplotlib.pyplot as plt

projectName = "documentation"
Demography  = toolkitHome.getToolkit(projectName=projectName,toolkitName="GIS_Demography")
katsrinBuildings = geopandas.read_file("KatsrinBuildings")
KatsrinCityOnly = Demography.analysis.createNewArea(shapeNameOrData=katsrinBuildings,dataSourceOrData="KatsrinDemography/Katsrin.shp").getData()

#######################
# Then, we load the LSM simulations.

LSM  = toolkitHome.getToolkit(projectName=projectName,toolkitName="LSM")
simulation = LSM.singleSimulation(resource="netcdf")

#######################
# Now, we will get an agent object using the example agent description.

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

######################
# Here we will demonstrate how the final number of casualties in each scenario can be attached to a single dataframe.
# In our example the examined variables are mass and protection policy, but every other variables such as agent type, release point, atmospheric stability and so on may be examined
# in the same manner.

windAngle = 0
x_coordinate = 263500
y_coordinate = 766750
datalist = []
Qs = [10*kg,20*kg]
policyNames = ["NoProtection","OpenRooms"]
indoor = dict(name="indoor",params=dict(turnover=15*min,enter="30s",stay="10min"))
for Q in Qs:
    Concentration = simulation.getConcentration(Q=Q)
    for protectionPolicy, policyName in zip([None, indoor],policyNames):
        if protectionPolicy is not None:
            Concentration = risk.ProtectionPolicy(actionList=protectionPolicy).compute(Concentration).squeeze().compute()
        riskAreas = Agent.RegularPopulation.calculate(Concentration, "C", isel={"datetime":-1})
        data = riskAreas.project(KatsrinCityOnly,[x_coordinate,y_coordinate],mathematical_angle=windAngle)
        data = data.groupby("severity")["effectedtotal_pop"].sum().reset_index()
        data["Q"] = Q.asNumber()
        data["Protection"] = policyName
        datalist.append(data)
data = pandas.concat(datalist)

#######################
# Now, we will plot the bar plot of the data, seperated by severity level, mass and protection policy.

x = "Protection"
hue = "Q"
fig,ax = plt.subplots()
for severities, injury, color in zip([["Severe","Light"],["Light"]],["All injuries","Light injuries"],["Reds_d","Greens_d"]):
    plotdata = data.query("severity in @severities").groupby([hue,x]).sum()["effectedtotal_pop"]
    plotdata = plotdata.to_frame().reset_index()
    ax = sns.barplot(x=x,y="effectedtotal_pop",data=plotdata,hue=hue,palette=color)
    ax.legend_.remove()
for p in ax.patches[0:2]:
    ax.text(p.get_x()+p.get_width()/2-0.055, p.get_height()+6,f"{Qs[0].asNumber()} kg",fontsize=8)
for p in ax.patches[2:4]:
    ax.text(p.get_x()+p.get_width()/2-0.055, p.get_height()+6,f"{Qs[1].asNumber()} kg",fontsize=8)
ax.set_xlabel("Protection")
ax.set_ylabel("Number of casualties")

#######################
# Finally, we may save the figure. In our example we will use the examined variables values in the name.

figName = "compareVariables"
for policyName in policyNames:
    figName += f"_{policyName}"
for Q in Qs:
    figName += f"_{Q.asNumber()}_kg"
figName += ".png"
# plt.savefig(figName)