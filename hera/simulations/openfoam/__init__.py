
import sys
version = sys.version_info[0]
if version==2:
    from .postprocess.extractVTK import VTKpipeline
    from .analysis.analysys import tests
else:
    from .process.process import process
    from .presentationLayer.Plotting import Plotting
    from .postprocess.dataManipulations import dataManipulations