Welcome to hera-data!

This project holds the relevant data for the toolkits and a script to load them to your project.

Available data
--------------
    - GIS_Buildings :
            * BNTL   - Building polygons shapefile of Israel

    - GIS_Topography :
            * BNTL    - Elvation shapefile of Israel
            * SRTM    - tif file of Israel
            * USGS

    - GIS_Topography:
            * BNTL - polygons shapefile of Israel


Listing data
------------

./hera-data-load list [data type]

Lists all the data types if [data type] is not supplied, else lists all the available databases.

Loading data
------------

./hera-data-load load [project Name] [data type] [data source]

Loads the data type/data source to the project.