__version__ = '1.0.0'

import sys
import os
import json
## Load modules if it is python 3.
version = sys.version_info[0]
if version==3:
    from .simulations import WRF
    #from .measurements import meteorological as meteo
    from .simulations import LSM
    from .simulations.LSM.DataLayer import SingleSimulation
    from .measurements import GIS
    from .simulations.interpolations.interpolations import spatialInterpolate
    #from .risk import riskassessment

from .simulations import openfoam

import logging
import logging.config

with open(os.path.join(os.path.dirname(__file__),'logging','heraLogging.config'),'r') as logconfile:
     log_conf_str = logconfile.read().replace("\n","")
     log_conf = json.loads(log_conf_str.replace("{herapath}",os.path.dirname(__file__)))

EXECUTION = 15
logging.addLevelName(EXECUTION, 'EXECUTION')

def execution(self, message, *args, **kws):
    self.log(EXECUTION, message, *args, **kws)

logging.Logger.execution = execution

logging.config.dictConfig(log_conf)


"""
 1.0.0
 -----
    - Introduced tools. 
        - The GIS tools work  
    - Project classes are equipped with a logger. 
    

 0.7.0
 -----
  - Adding the riskassessmet package. 
  - Adding the simulations/gaussian package. 
  - adding the simulation/evaporation package. 

 0.6.1
 -----
  - Fixing the simulations.interpolations package. 
  - Renanimg interpolation->spatialInterpolations.  
  
 0.6.0
 -----
  - adding tonumber,tounum and tometeorological/to mathematical functions to the utils. 

 0.5.1
 -----
 CampbellBinary parser and datalayer fixed.

 0.5.0
 -----
 Changed the meteorological structure(datalayer and presentaionlayer)

 0.4.1
 -----
 With demography in GIS

 0.4.0
 -----
 Added features to the turbulence calculator.
 Added options to the db documents search.

 0.3.0
 -----
 Changed the datalayer.analysis to datalayer.cache.
 Added more documentation.

 0.2.2
 -----
 Turbulence calculator working with sampling window None.

 0.2.1
 -----
 More turbulence calculator fix

 0.2.0
 -----
 Turbulence calculator fix

 0.1.1
 -----
 Removed the necessity to have a public DB

 0.1.0
 -----
 getData() from datalayer returns list of data.


 0.0.2
 -----
 
 LSM - * Tiding up the datalayer a bit 
       * LagrangianReader - changing the order of the x and y coordinates



"""
