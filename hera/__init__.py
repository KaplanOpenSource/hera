__version__ = '2.15.5'

import sys
python_version = sys.version_info
if python_version < (3, 8):
    raise NotImplementedError("Hera does not support Python versions lower than 3.8")

# Have some initial default logging configuration in case the user hasn't set any
from hera.utils.logging.helpers import initialize_logging
from hera.utils.logging import with_logger,initialize_logging, get_logger,getClassLogger,get_classMethod_logger
initialize_logging(disable_existing_loggers=False)

from hera.toolkit import ToolkitHome
from hera.datalayer import Project,datatypes
from hera.datalayer.autocache import cacheFunction,clearFunctionCache,clearAllFunctionsCache


toolkitHome = ToolkitHome()

"""
2.15.5
    #272    - Fixed the way we check the dataframe type in the analysis of the high freq toolkit. 
    
    #249    - Fixed the bug with JSON containing keys end with __type. Key must end with __type__. 
            - Fixed a bug in the processToJson: now also support dicts inside maps. 
2.15.4
    #193    - OpenFOAM: Uniting the write new fields and the write of the regular openfoam case. 
                        also include procedures to read the field from the case, only the boundaries and ect.
    #187    - Adding the support for VTk pipelines in the openfoam toolkit.   
    #182    - Added the default tile server in the tile toolkit. 
2.15.3
    #174     - Added the ability to activate ageneral function in the repository 
                TODO : add documentation to it. 
              - Refactoring the internals of the data toolkit. 
2.15.2 
    #168:   - Added the building toolkit. 
            - Still need to fix the lambda calculations. 
2.15.1
    #161:   - Adding the jupyter-lab-script that find the ip and runs the jupyter.
            - Adding toSTL for hera-GIS
            - Adding CRS points transformation.   
2.15
----
    #134:   - ExperimentToolkit: Adding support for storage of device data per device or per device type 
            - Adding the CLI for handling experiment. 
            - Adding documentation
    #122:   - Adding support for repositories. 
            - Repositories can be either datasources of the toolkit, or just regular measurements, cache, or simulations. 
            - The data in the repositories is automatically added when a project is created. 

    #119:   - Adding add_formatter,add_FileHandler that could be used 
              to define local log file handlers and formatters. 
              
            - Adding logging helpers.
            - Adding unum tutorial.  
            - Adding query tutorial 
            

2.14.1
------
    #109 : * Updating the stl creation in the topography and buildings toolkit. 

2.14.0
------
    #117 :  * Completing the management of projects. 
            * Completing the management of databases.
            * Completing the management of repositories.  
            * Adding the documentations.
            * Refactoring the documentation structure. 
    
    #117 : Updating the format of the hera data. 
           Adding a toolkit to handle the data, and use the standard loading to make sure that 
           all the data is loaded through this mechanis. Thus, the hera-data-load is redundent and 
           the project can be initialized with all the data using the CLI
           
           We have changed the format of the datasources. So it should also be updated with this version.  

          * a bug in the objects_createVerticesAndBoundary: 
            the boundaries were tight on the object and therefor it didnt remove the old boundary 
            of the domain which caused problems. Creating now constant 10cm gap in the model. 
            In the case that it is too large, change to be 0.1% of the length or something. 

    #108: * Updating the openFOAM CLI. 
          * Fixing a bug in the creation of the dispersion flow field. It was created with a local 
            directories and not with the global ones. 

    #110:   Updating the datalayer documentation and converting it to Jupyter-lab
            Moving the config from the toolkit to the project. 
2.13.2
--------

    #113:   * Adding a default read-only project that can be used to access general databases. 
            * When projectName is None in the contruction of the Project class, attempt to load the 
              project name from the caseConfiguration.json.  
        

2.13.0
--------
    #95:  Updating the changes for the experiment in 10.2023
    #97:  Updated the interface to the openfoam solver
    #102: 
        - Fixing the LSM to use the configuration to JSON. 

    #37: 
        - Adding CLI for StochasticLagrangian 
        - Changed the interface of ofobjects 
        - updated documentation 

2.10.0
------

    #38: - Rerwiting the hera-workflows to have simpler interface. 
         - Fixing bugs in dataframeutils and json utils. 
         - Changing the nonDBMetadataFrame interface.  

    #77: serialization and de-serialization of JSON that include unum objects.
          These functions handles the dictionary with unum object (to string and from string). 
    
        - changed to convertJSONtoConf to JSONToConfiguration 

2.9.1
------
     Issue#38: Updating the documentation. 
     
2.9.0
------
     #68: fix import of non existing packages like FreeCAD and hermes. 
     #65: Refix logging 
2.8.0
------
    #63: Issue#38- Working with openfoam 
    #61: Simpler logging method. 
    #34: A starting point for the unittest. 
    

2.7.0
------
    Issue#47: Updating the datalayer 
    Issue#48: The new requirement type
    PullRequest#51: Fix bug in the in demography  demography.analysis.createNewArea
    PullRequest#52: Initial setup files 
    Issue#13 (pull request #53): Create log file if it does not already exist
    Issue#55 (pull request #56): Fixed the bug, hope that well.
    PullRequest#58: Minor doc build fixes
    Issue#49/#59: Issues with building 

2.6.0
-----
    Issue#46: Updating the documentation to the datalayer. 
    Issue#44: Improved comparison of n-JSON files. 

2.5.2   
    ______________________
    Issue#40: DS changes to measurments.experiment
    Issue#24: Handling domain size #25
    issue#28: Adding experimentSetupWithData class to the experiment modoule
    Issue#26: RiskAssessment - Enhancement - Adding point wise risk calculations

2.5.1
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
        * Build a toolkit that unifies the access to NS and LSM.old  
        * Build a reader for the eulerian data, either decomposed or composed cases.
        * New interface to write and read OF fields (eulerian or lagrangian).
        * Improving the Hera interfaces     
 
    LSM.old - small fixes to handle parameters without units.
    Gaussian - Fixing the toolkit of large droplets 
             - Adding a constant meteorology
             
 Riskassessment 
 --------------
    Changing the computation of the risk areas back to contour from tricontour. 
 

    - Building new GIS toolkit . 
        - Depracating the cut of a region. 

2.4.1 
------

    LSM.old
    ----
    
    - Fixed some bugs in getSimulation and getSimulationTable

Datalayer
---------

- Project: Added get<simulation/cache/measurements.old>DocumentByID to the interface. 


2.4.0
-----

   -measurements.old.meteorology.highfreqdata:
    - Added an AveragingCalculator class to deal with TRH data
    - Major revision of MeanDataCalculator
    - Added a few new functions to singlePointTurbulenceStatistics
    - Renamed Monin-Obukhov length related quantities in singlePointTurbulenceStatistics to make it apparent that 
        they are calculated only from Sonic raw data

    LSM.old: 
        - template.run - change directory back to its location.
                       - changes the params of a template back to a separate branch. 
                       - changes the units of the parameters in the template to be more generic. 
                        
        - presentation layer - raise exception/print if there are no casualties in that direction.
                             - Fixed the rose to be with the right directions. 
                                 

2.3.0
-----
    - Added Ofir version
    - LSM.old:
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
    - measurements.old.demography was extended to allow loading an existing demography shp file. 
    
      

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

OF-LSM.old:
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
                
            LSM.old:
                - Updating Ustar/Hmix 
                - Fixed some bugs in the OF-LSM.old reading files.
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
    - fixed bug in the get concentration of the LSM.old
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
  - Changes to the intepolations in the simulations.old module. 
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
  - Adding the simulations.old/gaussian package. 
  - adding the simulation/evaporation package. 

 0.6.1
 -----
  - Fixing the simulations.old.interpolations package. 
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
 
 LSM.old - * Tiding up the datalayer a bit 
       * LagrangianReader - changing the order of the x and y coordinates



"""
