__version__ = '0.2.2'






import sys
version = sys.version_info[0]
if version==3:
    from .simulations import WRF
    from .measurements import meteorological as meteo
    from .simulations import LSM
    from .simulations.LSM.DataLayer import SingleSimulation
    from .measurements import GIS
from .simulations import openfoam

"""
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
