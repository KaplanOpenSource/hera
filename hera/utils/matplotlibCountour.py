def standardize_polygon(poly, units_conversion):
    """Scale polygon coordinates by a units conversion factor."""
    if isinstance(poly, list):
        xs = [p[0] for p in poly]
        ys = [p[1] for p in poly]
    else:
        xs = poly[:, 0]
        ys = poly[:, 1]
    return [(x * units_conversion, y * units_conversion) for (x, y) in zip(xs, ys)]


def toGeopandas(ContourData, inunits=None):
    
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
    from hera.utils.unitHandler import ureg, unumToPint

    try:
        from shapely import geometry
        import geopandas
    except ImportError:
        print("gis support not installed. ")
    inunits = inunits if inunits is not None else 1*ureg.m

    units_conversion = unumToPint(inunits).m_as(ureg.m)
    polyList = []
    levelsList = []
    if hasattr(ContourData, 'collections'):
        path_level_pairs = [
            (contour_path, level)
            for col, level in zip(ContourData.collections, ContourData.levels)
            for contour_path in col.get_paths()
        ]
    else:
        path_level_pairs = list(zip(ContourData.get_paths(), ContourData.levels))

    for contour_path, level in path_level_pairs:
        polygons = (standardize_polygon(p, units_conversion) for p in contour_path.to_polygons())
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

def plot_polygons(ax, geom, facecolor="none", edgecolor="black", linewidth=1):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MplPolygon
    from matplotlib.collections import PatchCollection

    from shapely.geometry import (
        Polygon, MultiPolygon,
        LineString, MultiLineString,
        Point, MultiPoint,
        GeometryCollection
    )

    if geom.is_empty:
        return

    if isinstance(geom, Polygon):
        # exterior
        patch = MplPolygon(
            list(geom.exterior.coords),
            closed=True,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=linewidth,
        )
        ax.add_patch(patch)

        # holes
        for interior in geom.interiors:
            x, y = interior.xy
            ax.plot(x, y, color=edgecolor, linewidth=linewidth)

    elif isinstance(geom, MultiPolygon):
        for part in geom.geoms:
            plot_polygons(ax, part, facecolor, edgecolor, linewidth)

    elif isinstance(geom, LineString):
        x, y = geom.xy
        ax.plot(x, y, color=edgecolor, linewidth=linewidth)

    elif isinstance(geom, MultiLineString):
        for part in geom.geoms:
            plot_polygons(ax, part, facecolor, edgecolor, linewidth)

    elif isinstance(geom, Point):
        ax.plot(geom.x, geom.y, "o", color=edgecolor)

    elif isinstance(geom, MultiPoint):
        for part in geom.geoms:
            plot_polygons(ax, part, facecolor, edgecolor, linewidth)

    elif isinstance(geom, GeometryCollection):
        for part in geom.geoms:
            plot_polygons(ax, part, facecolor, edgecolor, linewidth)

    else:
        raise TypeError(f"Unsupported geometry type: {geom.geom_type}")
