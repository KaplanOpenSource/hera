__version__ = (0,0,2)

"""
0.0.2
-----
   * Added the breathing rate considarations. 

0.0.1 
----- 
A package to calculate the doasge and project them on 
geographical area (+ calculating the casualties).

  * Defines the 3 tyes of thresholds: 
     - Threshold 
     - Lognormal dose response 
     - Exponential dose response. 

  * Define 3 types of calculators: 
     - Haber 
     - TenBerge. 
     - MaxConcentration. 
"""

from .Injury import InjuryFactory
injuryfactory = InjuryFactory()
