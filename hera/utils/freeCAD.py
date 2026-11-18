import numpy
import warnings

# The import failure is a warning rather than an error to allow partial usage
# of hera.utils; the names stay defined (as None) so that the functions below
# can say what is missing instead of raising NameError.
try:
    import FreeCAD
except ImportError:
    FreeCAD = None

try:
    import Mesh
except ImportError:
    Mesh = None

if FreeCAD is None or Mesh is None:
    warnings.warn("freecad is not installed. some features will not work.")


def getObjFileBoundaries(fileName):
    """
            Loads the file in FreeCAD and finds its boundaries.
            Can work with anyfile freeCAD recognizes (stl,obj,...).
    Parameters
    ----------
    fileName :  str

    Returns
    -------
        dict
            with keys:
                - 'XMax', 'YMax', 'ZMax'
                - 'XMin', 'YMin', 'ZMin'
    """
    if FreeCAD is None or Mesh is None:
        raise ImportError("freecad is not installed. please install it before trying again.")

    Mesh.open(fileName)
    objFile = FreeCAD.getDocument("Unnamed")

    bboxes = [x.Mesh.BoundBox for x in objFile.findObjects()]
    nameMapping = dict(XMax = 'XMax',
                       YMax = 'YMax',
                       ZMax = 'ZMax',
                       XMin='XMin',
                       YMin='YMin',
                       ZMin='ZMin')

    maxPropList = ['XMax', 'YMax', 'ZMax']
    corners = dict()
    for propName in maxPropList:
        corners[nameMapping[propName]] = numpy.max([getattr(x, propName) for x in bboxes]) / 1000

    minPropList = ['XMin', 'YMin', 'ZMin']
    for propName in minPropList:
        corners[nameMapping[propName]] = numpy.min([getattr(x, propName) for x in bboxes]) / 1000

    return corners