import sys
version = sys.version_info[0]
if version > 2:
    from .datalayer.OFObjects import ofObjectHome
    from .datalayer.hermesWorkflow import Workflow_Flow

    OFObjectHome = ofObjectHome()


DECOMPOSED_CASE = 'Decomposed Case'
RECONSTRUCTED_CASE ="Reconstructed Case"


TYPE_VTK_FILTER = "vtk_filter"

from .datalayer.hermesWorkflow import Workflow_Flow,Workflow_Dispersion,HERAMETADATA
