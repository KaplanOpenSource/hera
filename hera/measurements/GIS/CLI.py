import pandas as pd
from hera import toolkitHome
from hera.utils import WSG84, ITM
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

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_VECTOR_TOPOGRAPHY, projectName=projectName)
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

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_RASTER_TOPOGRAPHY, projectName=projectName)
    print(f"Loaded topography (raster) datasources for the project {tk.projectName}")
    print(tk.getDataSourceTable())


def topography_raster_toSTL(arguments):
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

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_RASTER_TOPOGRAPHY, projectName=projectName)
    logger.info(f"Working on project {projectName}")

    dxdy = arguments.dxdy if "dxdy" in arguments else 30
    inputCRS = WSG84 if arguments.inputCRS is None else arguments.inputCRS
    outputCRS = ITM if arguments.outputCRS is None else arguments.outputCRS
    dataSourceName = None if arguments.dataSourceName is None else arguments.dataSourceName

    stlString = tk.getDomainElevation_STL(minx=arguments.minx,
                                          miny=arguments.miny,
                                          maxx=arguments.maxx,
                                          maxy=arguments.maxy,
                                          dxdy=dxdy,
                                          inputCRS=inputCRS,
                                          outputCRS=outputCRS,
                                          dataSourceName=dataSourceName)
    fileName = arguments.fileName
    if '.stl' not in fileName:
        fileName = f"{fileName}.stl"

    logger.info(f"Writing STL to the file:  {fileName}")
    with open(fileName, "w") as outSTLFile:
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

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_BUILDINGS, projectName=projectName)
    print(f"Loaded buildings datasources for the project {tk.projectName}")
    print(tk.getDataSourceTable())


def buildings_raster_toSTL(arguments):
    logger = logging.getLogger("hera.bin.measuerments.GIS")

    if "projectName" in arguments:
        projectName = arguments.projectName
    else:
        projectName = None

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_BUILDINGS, projectName=projectName)
    logger.info(f"Working on project {projectName}")

    dxdy = arguments.dxdy if "dxdy" in arguments else 30
    inputCRS = WSG84 if arguments.inputCRS is None else arguments.inputCRS
    outputCRS = ITM if arguments.outputCRS is None else arguments.outputCRS
    dataSourceName = None if arguments.dataSourceName is None else arguments.dataSourceName

    buildings = tk.getBuildingsFromRectangle(minx=arguments.minx,
                                             miny=arguments.miny,
                                             maxx=arguments.maxx,
                                             maxy=arguments.maxy,
                                             withElevation=True)
    fileName = arguments.fileName

    if '.stl' not in fileName:
        fileName = f"{fileName}.stl"

    tk.buildingsGeopandasToSTLRasterTopography(buildings, "BLDG_HT", "elevation", fileName)
    logger.info(f"Writing STL to the file:  {fileName}")


def get_landocver(arguments):
    logger = logging.getLogger("hera.measuerments.GIS.CLI.get_landocver")
    logger.info(f"Arguments: {arguments}")

    if "projectName" not in arguments:
        arguments.projectName = None

    if "dataSourceName" not in arguments:
        arguments.dataSourceName = None

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_LANDCOVER, projectName=arguments.projectName)

    dxdy = int(arguments.dxdy) if "dxdy" in arguments else 30
    inputCRS = WSG84 if arguments.inputCRS is None else int(arguments.inputCRS)
    # outputCRS = ITM if arguments.outputCRS is None else int(arguments.outputCRS)
    dataSourceName = None if arguments.dataSourceName is None else arguments.dataSourceName
    windDirection = None if arguments.windDirectionis is None else float(arguments.windDirectionis)
    resolution = None if arguments.resolution is None else float(arguments.resolution)


    if arguments.roughness:
        xarray = tk.getRoughness(arguments.minx,
                        arguments.miny,
                        arguments.maxx,
                        arguments.maxy,
                        dxdy=dxdy,
                        inputCRS=inputCRS,
                        dataSourceName=dataSourceName,
                        isBuilding=arguments.isBuilding,
                        windMeteorologicalDirection=windDirection,
                        resolution=resolution)
    else:
        xarray = tk.getLandCover(arguments.minx,
                        arguments.miny,
                        arguments.maxx,
                        arguments.maxy,
                        dxdy=dxdy,
                        inputCRS=inputCRS,
                        dataSourceName=dataSourceName)

    if 'z0' in xarray.coords:
        df = pd.DataFrame({
            'lat': xarray['lat'].values.ravel(),
            'lon': xarray['lon'].values.ravel(),
            'landcover': xarray['landcover'].values.ravel(),
            'z0': xarray['z0'].values.ravel()
        })
    else:
        df = pd.DataFrame({
            'lat': xarray['lat'].values.ravel(),
            'lon': xarray['lon'].values.ravel(),
            'landcover': xarray['landcover'].values.ravel()
        })

    df.to_csv(arguments.filePath,index=False)


