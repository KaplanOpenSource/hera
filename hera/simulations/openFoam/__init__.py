import warnings
import sys
version = sys.version_info[0]
if version > 2:
    from .datalayer.OFObjects import ofObjectHome
    try:
        from .datalayer.hermesWorkflow import Workflow_Flow
    except:
        warnings.warn("hermes is not installed. some features will not work.")

    OFObjectHome = ofObjectHome()


DECOMPOSED_CASE = 'Decomposed Case'
RECONSTRUCTED_CASE ="Reconstructed Case"


TYPE_VTK_FILTER = "vtk_filter"

try:
    from .datalayer.hermesWorkflow import Workflow_Flow,Workflow_Dispersion,HERAMETADATA
except:
    warnings.warn("hermes is not installed. some features will not work.")
