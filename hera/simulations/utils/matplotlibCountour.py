from shapely import geometry
import geopandas

from unum.units import *


def toGeopandas(ContourData, inunits=m):
    """
        Converts the contours of matplotlib to polygons.


    :param ContourData:
                The output of a matplotlib counrour (maybe also contourf)
    :param inunits:
        The output will be in meters, but input can be in other units.
        So use this to generate a utils factor.
    :return:
        A geopandas object with the contours as polygons and levels as attributes.

    """
    unitsconversion = inunits.asNumber(m)
    polyList = []
    levelsList = []
    for col, level in zip(ContourData.collections, ContourData.levels):
        # Loop through all polygons that have the same intensity level
        for contour_path in col.get_paths():
            # Create the polygon for this intensity level
            # The first polygon in the path is the main one, the following ones are "holes"
            for ncp, cp in enumerate(contour_path.to_polygons()):
                if isinstance(cp,list):
                    x = [x[0] for x in cp]
                    y = [x[1] for x in cp]
                else:
                    x = cp[:, 0]
                    y = cp[:, 1]
                new_shape = geometry.Polygon([(i[0]* unitsconversion, i[1]* unitsconversion) for i in zip(x, y)])
                if ncp == 0:
                    poly = new_shape
                else:
                    # Remove the holes if there are any
                    poly = poly.difference(new_shape)
                    # Can also be left out if you want to include all rings
            levelsList.append(level)
            polyList.append(poly)

    ret = geopandas.GeoDataFrame({"Level": levelsList, "contour": polyList}, geometry="contour")
    return ret
