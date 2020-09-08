import geopandas
import shapely
import pandas

def PolygonDataFrameIntersection(dataframe, polygon):
    """
    Creates a new dataframe based on the intersection of a dataframe and a polygon.

    Parameters:
    ----------
    dataframe: A geopandas dataframe.
    polygon: A shapely polygon

    Returns:
    --------

    geopandas.dataframe.

    A new geopandas dataframe

    """

    newlines = []
    for line in dataframe["geometry"]:
        newline = polygon.intersection(line)
        newlines.append(newline)
    dataframe["geometry"] = newlines
    dataframe = dataframe[~dataframe["geometry"].is_empty]

    return dataframe

def getBoundaries(doc):
    """
    Returns a dict with  the boundaries of the document.

    Parameters:
    -----------

    doc:

    Returns:
    ---------

    dict with the boundaries.
    """

    dataframe = doc.getData()
    points = doc.desc["points"]

    boundaries = dict(xmin=points[0], xmax=points[2], ymin=points[1], ymax=points[3], zmin=dataframe["HEIGHT"].min(), zmax=dataframe["HEIGHT"].max())

    return boundaries

def makePolygonFromEndPoints(points):

    polygon = shapely.geometry.Polygon([[points[0],points[1]],
                                        [points[0],points[3]],
                                        [points[2],points[3]],
                                        [points[2],points[1]]])
    return polygon

def ConvexPolygons(data, buffer=100):
    """
    Returns polygons of groups of buildings.
    """
    data = data.reset_index()
    d = data.buffer(buffer)
    indicelist=[[0]]
    for i in range(1,len(data)):
        found = False
        for g in range(len(indicelist)):
            for n in indicelist[g]:
                if d[i].intersection(d[n]).is_empty:
                    continue
                else:
                    indicelist[g].append(i)
                    found = True
                    break
            if found:
                break
            if g==len(indicelist)-1:
                indicelist.append([i])

    geo = data.loc[indicelist[0]].unary_union.convex_hull
    gpd = geopandas.GeoDataFrame.from_dict([{"geometry":geo,"area":geo.area}])
    for indice in indicelist[1:]:
        geo = data.loc[indice].unary_union.convex_hull
        gpd = pandas.concat([gpd,geopandas.GeoDataFrame.from_dict([{"geometry":geo,"area":geo.area}])])

    gpd = gpd.sort_values(by="area", ascending=False).reset_index()
    found=False
    for i in range(len(gpd)):
        for j in range(i+1,len(gpd)):
            if gpd.loc[i].geometry.intersection(gpd.loc[j].geometry).is_empty:
                continue
            else:
                found = True
                break
        if found:
            break
    if found:
        gpd = ConvexPolygons(gpd,buffer=1)

    return gpd


