import warnings
from .datalayer.OFObjects import OFObjectHome
from .datalayer.hermesWorkflow import Workflow_Eulerian,Workflow_Lagrangian

OFObjectHome = OFObjectHome()


DECOMPOSED_CASE = 'Decomposed Case'
RECONSTRUCTED_CASE ="Reconstructed Case"


TYPE_VTK_FILTER = "vtk_filter"

try:
    from .datalayer.hermesWorkflow import Workflow_Eulerian, Workflow_Lagrangian
except ImportError:
    warnings.warn("hermes is not installed. some features will not work.")
