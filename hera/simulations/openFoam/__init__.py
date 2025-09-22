import warnings
from .OFWorkflow import workflow_Eulerian,workflow_Lagrangian

CASETYPE_DECOMPOSED = 'Decomposed Case'
CASETYPE_RECONSTRUCTED = "Reconstructed Case"
TYPE_VTK_FILTER = "vtk_filter"


FLOWTYPE_COMPRESSIBLE = "compressible"
FLOWTYPE_INCOMPRESSIBLE = "incompressible"
FLOWTYPE_DISPERSION = "dispersion"

FIELDTYPE_SCALAR = "scalar"
FIELDTYPE_VECTOR = "vector"
FIELDTYPE_TENSOR = "tensor"

FIELDCOMPUTATION_EULERIAN = "eulerian"
FIELDCOMPUTATION_LAGRANGIAN = "lagrangian"