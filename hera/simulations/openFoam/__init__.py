import sys
version = sys.version_info[0]
if version > 2:
    from .datalayer.OFObjects import ofObjectHome
    from .datalayer.hermesWorkflow import Workflow_Flow

    OFObjectHome = ofObjectHome()