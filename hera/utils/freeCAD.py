import numpy
import warnings

try:
    from freecad import app as FreeCAD
    import Mesh
except ImportError:
    # I changed the raise error to warning to allow partial usage.
    #raise ImportError("freecad is not installed. please install it before trying again.")
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
    Mesh.open(fileName)
    objFile = FreeCAD.getDocument("Unnamed")

    bboxes = [x.Mesh.BoundBox for x in objFile.findObjects()]

    maxPropList = ['XMax', 'YMax', 'ZMax']
    corners = dict()
    for propName in maxPropList:
        corners[propName] = numpy.max([getattr(x, propName) for x in bboxes]) / 1000

    minPropList = ['XMin', 'YMin', 'ZMin']
    for propName in minPropList:
        corners[propName] = numpy.min([getattr(x, propName) for x in bboxes]) / 1000

    return corners