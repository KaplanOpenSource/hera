from shapely import geometry
import geopandas

from hera.utils.unitHandler import  *

def standardize_polygon(poly, units_conversion):
    if isinstance(poly, list):
        xs = [p[0] for p in poly]
        ys = [p[1] for p in poly]
    else:
        xs = poly[:, 0]
        ys = poly[:, 1]
    return [(x * units_conversion, y * units_conversion) for (x, y) in zip(xs, ys)]


def toGeopandas(ContourData, inunits=ureg.m):
    """
        Converts the contours of matplotlib to polygons.


    :param ContourData:
                The output of a matplotlib counrour (maybe also contourf)
    :param inunits:
        The output will be in meters, but input can be in other units.x
        So use this to generate a utilsOld factor.
    :return:
        A geopandas object with the contours as polygons and levels as attributes.

    """
    units_conversion = inunits.asNumber(ureg.m)
    polyList = []
    levelsList = []
    for col, level in zip(ContourData.collections, ContourData.levels):
        # Loop through all polygons that have the same intensity level
        for contour_path in col.get_paths():
            polygons = (standardize_polygon(p,units_conversion) for p in contour_path.to_polygons())
            try:
                # break the list -- "shell" takes the first, "holes" takes the rest
                (shell, *holes) = polygons
            except ValueError:
                # There was nothing in the list
                pass
            else:
                # There was "shell" and maybe some "holes"
                shape = geometry.Polygon(shell, holes)
                levelsList.append(level)
                polyList.append(shape)


    ret = geopandas.GeoDataFrame({"Level": levelsList, "contour": polyList}, geometry="contour")
    return ret


