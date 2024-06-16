from hera import toolkitHome
from hera.utils import WSG84,ITM
import logging

def topography_vector_list(arguments):
    """
        Lists the topography datasources.

    Parameters
    ----------
    arguments:
        No arguments needed

    Returns
    -------

    """
    if "projectName" in arguments:
        projectName = arguments.projectName
    else:
        projectName = None

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_VECTOR_TOPOGRAPHY,projectName=projectName)
    print(f"Loaded topography (vector) datasources for the project {tk.projectName}")
    print(tk.getDataSourceTable())

def topography_raster_list(arguments):
    """
        Lists the topography datasources.

    Parameters
    ----------
    arguments:
        No arguments needed

    Returns
    -------

    """
    if "projectName" in arguments:
        projectName = arguments.projectName
    else:
        projectName = None

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_RASTER_TOPOGRAPHY,projectName=projectName)
    print(f"Loaded topography (raster) datasources for the project {tk.projectName}")
    print(tk.getDataSourceTable())

def topography_raster_topography(arguments):
    """
        Gets the coordinates and saves an STL file.

    Parameters
    ----------
    arguments:
        - left
        - right
        - top
        - bottom
        - inputCRS
        - outpuCRS [default : ITM]
        - dataSourceName
        - outputFilename
        - projectName
    Returns
    -------

    """
    logger = logging.getLogger("hera.bin.measuerments.GIS")

    if "projectName" in arguments:
        projectName = arguments.projectName
    else:
        projectName = None

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_RASTER_TOPOGRAPHY,projectName=projectName)
    logger.info(f"Working on project {projectName}")

    dxdy = arguments.dxdy if "dxdy" in arguments else 30
    inputCRS = WSG84 if arguments.inputCRS is None else arguments.inputCRS
    outputCRS = ITM if arguments.outputCRS is None else arguments.outputCRS
    dataSourceName = None if arguments.dataSourceName is None else arguments.dataSourceName

    stlString = tk.getDomainElevation_STL(left=arguments.left,
                              bottom=arguments.bottom,
                              right=arguments.right,
                              top=arguments.top,
                              dxdy=dxdy,
                              inputCRS=inputCRS,
                              outputCRS=outputCRS,
                              dataSourceName=dataSourceName)
    fileName =arguments.fileName
    if '.stl' not in fileName:
        fileName = f"{fileName}.stl"

    logger.info(f"Writing STL to the file:  {fileName}")
    with open(fileName,"w") as outSTLFile:
        outSTLFile.write(stlString)



def buildings_parser_list(arguments):
    """
        Lists the topography datasources.

    Parameters
    ----------
    arguments:
        No arguments needed

    Returns
    -------

    """
    if "projectName" in arguments:
        projectName = arguments.projectName
    else:
        projectName = None

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_BUILDINGS,projectName=projectName)
    print(f"Loaded buildings datasources for the project {tk.projectName}")
    print(tk.getDataSourceTable())