__version__=(1,1,0)

from unum.units import *
from unum import Unum
from unum import NameConflictError


try:
	dosage = Unum.unit('dosage',mg*min/m**3,"Exposure dosage")
# except NameConflictError:
except NameConflictError:
	pass

tonumber = lambda x,theunit: x.asNumber(theunit) if isinstance(x,Unum) else x
tounum = lambda x,theunit: x.asUnit(theunit) if isinstance(x,Unum) else x*theunit


def toMeteorologicalAngle(mathematical_angle): 
	return (270-mathematical_angle) if ((270-mathematical_angle)   >= 0)   else (630-mathematical_angle)
def toMatematicalAngle(meteorological_angle): 
	return (270-meteorological_angle) if ((270-meteorological_angle) >= 0) else (630 - meteorological_angle)

def toAzimuthAngle(mathematical_angle): 
	return (90-mathematical_angle) if ((90-mathematical_angle) >= 0) else (450 - mathematical_angle)

from .riskDatalayer.riskDatalayer import riskDatalayer
from .agents.Agents import Agent
from .presentation.casualtiesFigs import casualtiesPlot
from .protectionpolicy.ProtectionPolicy import  ProtectionPolicy
from .analysis.riskAreas import getRiskAreaAlgorithm

from .agents.effects.thresholdGeoDataFrame import thresholdGeoDataFrame

casualtiesPlots = casualtiesPlot()

#from pynumericalmodels.LSM import  datalayer as datalayerLSM
from hera import SingleSimulation as datalayerLSM

#from pynumericalmodels.utils.matplotlibCountour import toGeopandas
from hera.simulations.utils.matplotlibCountour import toGeopandas

dispersion = {  "LSM" : datalayerLSM }

#from pymeteoGIS.demography.Population import projectIsoplethOnPopulation
#from pymeteoGIS import GISResources
#from hera.GIS import population

"""
	Version (1,1,0)
	---------------
    - Adding the HERA system.   


	Version (1,0,13)
	---------------
    - Fixed the bug with the dt in the indoor protectionPolicy  
    
	Version (1,0,12)
	---------------
	- Fixed a bug in the maxmial concentration 
	- added Metyl phtalate. 

	Version (1,0,11)
	---------------
	- Fixed the switch backend.   

	Version (1,0,10)
	---------------
	- Fixed the bug in calculation of raw  

	Version (1,0,9)
	---------------
	- Fixed some bug in the Threshold injury. return numpy.array and not a list.  

	Version (1,0,8)
	---------------
	- Fixed some bug in the Threshold injury. Maily relate to units and getPercent. 


	Version (1,0,7)
	---------------
	- ProtectionPolicy no accepts a list of actions as a json dict in its 
	  contructor.     

	Version (1,0,6)
	---------------
	- copied the attrs from the data in the policy.    

	Version (1,0,4)
	---------------
	- Fixed the bug in the breathing rate.   

	Version (1,0,3)
	---------------
	- Fixed the single level bug in the threshold injury  

	Version (1,0,1)
	---------------
	- Added the calculateRaw to calculate effect from pandas. 
	
	Version (1,0,0)
	---------------
	- Must supply the units of the concentration or the  code will break.
	 

	Version (0,3,6)
	---------------
 	- fixing the calculator to get the time column in the calculation and not during the init. 

	Version (0,3,5)
	--------------
	- Fixing the 'calcualte' to use the x,y columns from the user. 


	Version (0,3,4)
	--------------
	- Fixed the estimation of the bounds in the risk area calculations when multiple time 
	  steps are involved. 
	


	Version (0,3,3)
	--------------

	- We add the levels of higher toxicit y to solve the problem
	 of a wide gape between the two severities (which led to double counting the casualties). 

	Version (0,3,2)
	--------------
	- fixing a bug in the presentation of the casualties - the ax selection problem.
	
	Version (0,3,1)
	--------------
	-  Fixed the bug of the returned projected in casualty figs. 
	-  Fixed a bug in the meteorological->mathematical conversion. 

	Version (0,3,0)
	---------------
	- Added the area Risk module to analyze the risk in missiles that fall within a certain area. 

	Version (0,2,0) 
	---------------

	- Added a protection policy module. 
            Computes the concentration indoor and with/without masks. 

	Version (0,1,0) 
	---------------

	- Finished the plotting module .
	-

	Version (0,0,2)
	---------------
	- More debugging. 
	- extended the project to a list of directions. 
	- added rose plots


	Version (0,0,1)
	---------------
	- works with threshold and lognormal dose response 
	- implements Tenberge (classic), Haber and Maximal concentrations. 
	- Assumes that the concentration results are in m. 
	

"""
