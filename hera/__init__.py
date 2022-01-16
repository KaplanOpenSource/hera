__version__ = '2.4.1'

import os
import json
import sys
version = sys.version_info[0]
if version > 2:
    from .toolkit import ToolkitHome
    toolkitHome =ToolkitHome()

    # Setup the units for the model
    # The project name is not imporant becuase logging is universal for user.
    loggingHome = toolkitHome.getToolkit(toolkitName="Logging",projectName=None)

"""

 GIS_Buildings 
 --------------
    - Adding the i0,j0 index of the bloc to the lambda. 
    - canging names. 
    - refactoring the canopywindProfile. 

 Util 
 ----
    - Adding special units and fixing the multiple load problem under the unum 
    - Adding logging toolkit to handle the local logs. 
             + Changed the default logging directory to .pyhera/log. 
    
 Simulations
 -----------

    - Updating the hermes/hera interafces. 
    - hera workflow - moved to the hermes. 
    - updated the interfaces and the work on the hera-of package. 
    
    OpenFoam  
        * Build a toolkit that unifies the access to NS and LSM  
        * Build a reader for the eulerian data, either decomposed or composed cases.
        * New interface to write and read OF fields (eulerian or lagrangian).
        * Improving the Hera interfaces     
 
    LSM - small fixes to handle parameters without units.
    Gaussian - Fixing the toolkit of large droplets 
             - Adding a constant meteorology
             
 Riskassessment 
 --------------
    Changing the computation of the risk areas back to contour from tricontour. 
 

    - Building new GIS toolkit . 
        - Depracating the cut of a region. 

2.4.1 
------

    LSM
    ----
    
    - Fixed some bugs in getSimulation and getSimulationTable

Datalayer
---------

- Project: Added get<simulation/cache/measurements>DocumentByID to the interface. 


2.4.0
-----

   -measurements.meteorology.highfreqdata:
    - Added an AveragingCalculator class to deal with TRH data
    - Major revision of MeanDataCalculator
    - Added a few new functions to singlePointTurbulenceStatistics
    - Renamed Monin-Obukhov length related quantities in singlePointTurbulenceStatistics to make it apparent that 
        they are calculated only from Sonic raw data

    LSM: 
        - template.run - change directory back to its location.
                       - changes the params of a template back to a separate branch. 
                       - changes the units of the parameters in the template to be more generic. 
                        
        - presentation layer - raise exception/print if there are no casualties in that direction.
                             - Fixed the rose to be with the right directions. 
                                 

2.3.0
-----
    - Added Ofir version
    - LSM:
        * Added simulation name to each run. 
        * template.run returns a single simulation
        * Added TRUE, FALSE constants.   

    - Removed some last toNumber (some remain in the Gaussian toolkit). 
    - Removed the ProjectMultiDB. and MultiDBPublic. 

2.2.0
-----
    - removed the simulation.util.toUnum and toNumber. Use util tounit and tonumber. 
    - Extended the GIS.location.raster to save the image locally only if it is not a local image. 
    - Added the TOOLKIT save modes to the toolkit home. 
    - measurements.demography was extended to allow loading an existing demography shp file. 
    
      

2.1.4
-----
    - Removing the big files from the history. 

2.1.3
-----
    - Few corrections of the documentation and some minor bugs. 

2.1.2
-----
   - Fixing some bugs in the previous version 
   - Added printing to formatted text of an agent. 

2.1.1
-----


buildings:
 - refactoring the code

OF-LSM:
 - fixed small bug in reading points

Topography:
 - Removed extra code    
    
Experiment: 
     -   Fixed the dynamic load of the experiment.  

2.1.0
-----
     
    - GIS
        Topography: 
                - Add height. Adds the topography height to a regular data. 
                              The topography is interpolated from pandas (with xarray).   
              
    - openFOAM : 
            NavierStokes
                - Adding Canopy profile.
                
            LSM:
                - Updating Ustar/Hmix 
                - Fixed some bugs in the OF-LSM reading files.
                - Changed getSource to makeSource
               
    - datalayer:
        - project to specify the logger name.         
        - adding 'dict' data format to store dict as a resource. 
    
2.0.2
-----
    - Examples gallery for the risk assessment  
    - Fixing the pvServerRun 

2.0.1
-----
    - Example for the toolkits raster 
    - fixed bug in the get concentration of the LSM
    - Updating the documentaiton 

2.0.0
-----
    - Fixing the commit
    - Changing the structure of the toolkits. 
    - Fixing the imports to be lighter
    - Some other changes to make the risk assessment procedure work.  
    - Updating the AgentsHome according to the existing agents description
    - changed wind_speed to horizontal_wind_speed in the turbulence calculator. 

 1.1.3
------
   - Updates of the documentation. 
   - Minor refactoring of the utils. 

 1.1.2
------
  - Changes to the intepolations in the simulations module. 
  - updated documentations. 
   

 1.1.1
------ 

    - Changes to the  command lines. 
    - Rearranging the meteorology. 


 1.0.0
 -----
    - Introduced tools. 
        - The GIS tools work  
    - Project classes are equipped with a logger. 

 This version consists a structural change that introduces concept of Tools. 
 
 A tool is a set of library functions that is designed to handle a single type of data. 
 Many tools also include a capability to search for public data automatically. 
 
 The change is cosmetic, but also includes several new concepts in tools such as datasource. 
 A datasource is the name of data that will be used by default by the tool. For examplt the BNTL data 
 is the default datasource of the mesasurements.GIS.topography tool. In this 
 example, the default datasource is stored in the public database. 
 
 0.7.0
 -----
  - Adding the riskassessment package. 
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
 Changed the meteorology structure(datalayer and presentaionlayer)

 0.4.1
 -----
 With demography in GIS

 0.4.0
 -----
 Added features to the turbulence calculator.
 Added options to the db documents search.

 0.3.0
 -----
 Changed the datalayer.datalayer to datalayer.cache.
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
