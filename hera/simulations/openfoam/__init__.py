
import sys
version = sys.version_info[0]
if version==2:
    from .postprocess.extractVTK import VTKpipeline
else:
    from .process.process import process
    from .analysis.Plotting import Plotting
    from .postprocess.dataManipulations import dataManipulations