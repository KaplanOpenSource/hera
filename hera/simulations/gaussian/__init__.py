"""

Gaussian module.
----------------

Calculates the semi-analytical gaussian models.

TODO:
- Add cloud particle ditribution.


"""

from .Sigma import briggsRural
from .MeshUtils import gaussianToMesh,GaussianIntegrationToMesh
from .DropletCloud import FixedPositionDropletsCloud,LinePositionDropletsCloud
from .Meteorology import meteorologyFactory




