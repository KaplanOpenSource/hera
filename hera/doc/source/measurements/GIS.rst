GIS
===

The GIS method is used to manage GIS data.
Its main use is adding geopandas dataframes, polygons, points and picturs to the database
and easily retrieving them using simple queries.
In addition, it enables common data manipulations as built-in functions, and helps
building stl files based on the stored data.



10-minute tutorial
------------------

This tutorial demonstrate how to store and retrieve data from the database.

Geopandas data
**************

first, import the requested module:

.. code-block:: python

    from hera import GIS

The module used for storing and retrieving data is called GIS_datalayer.

.. code-block:: python

    datalayer = GIS.GIS_datalayer(projectName = "Example", FilesDirectory = "/home/ofir/Projects/2020/GIS")

Storing and retrieving GIS data is done using the same function.
The first time one asks for data of a certain rectangle in Israels' map,
one has to deliver the rectangle's coordinates and a name for the new file that saves it.
Any additional parameters that describe the data may be added, such as the region's name.

.. code-block:: python

    points = [193000, 731000, 218500, 763000] # The rectangle coordinantes
    #example_data = datalayer.getGISDocuments(points=points, CutName="Example", Region="Haifa")

The function returns a list of the documents that fit the parameters.
The data itself can be achieved by choosing a document and using the function getData().

For retrieving data that already exists, one can deliver any parameters that describe the data,
or even no parameters at all for getting all the GIS geopandas data in the project.
for example, for getting a specific data one may use a variation of the following line:

.. code-block:: python

    datalayer.getGISDocuments()[3].getData()


Shapely Data
************

Shapely points and polygons may be added to the database.
Adding shapes is done like this:

.. code-block:: python

    from shapely import geometry
    Haifa_Port = geometry.Point([200000,748000])
    dataHandler.addGeometry(name="Haifa_Port", Geometry=Haifa_Port)

The shape is retrieved this way:

.. code-block:: python

    shape = dataHandler.getGeometry(name="Haifa_Port")

It can also be used in oreder to retrieve a goeopandas dataframe that contains the shape:

.. code-block:: python

    datalayer.getGISDocuments(Geometry="Haifa_Port"])[0].getData()

Demography
**********

Demography data in an area may be retrieved.
If we have a document whose name or CutName is "TelAviv", for example,
which holds coordinates that may define a polygon,
we can find the amount of people who live in that polygon.

.. code-block:: python

    population = GIS.population(projectName = "Example")
    data = population.projectPolygonOnPopulation(Geometry="TelAviv")

Usage
-----

.. toctree::
    :maxdepth: 1

    GIS/examples/ManagingData
    GIS/examples/Plotting
    GIS/examples/Convertions
    GIS/examples/DataManipulations
    GIS/examples/Population
