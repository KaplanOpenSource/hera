import numpy
import geopandas
import pydoc
from shapely import affinity
from shapely.geometry import GeometryCollection,LineString,point
import pandas
from itertools import product,chain
import multiprocessing
from functools import partial

from ...utils import toMeteorologicalAngle, toMathematicalAngle
def getRiskAreaAlgorithm(algorithmName,**kwargs):
	"""
		Return an estimator class. 
		kwargs are the parameters passed  as the constructor. 

		currently only Sweep is implemented
	"""
	estimatorCLS = pydoc.locate("pyriskassessment.datalayer.riskAreas.riskAreaAlgorithm_%s" % algorithmName.title())
	if estimatorCLS is None: 
		mymod = pydoc.locate("pyriskassessment.datalayer.riskAreas")
		raise ValueError("estimator %s not found. Available estimators are: %s" % (algorithmName,",".join([x.split("_")[1] for x in dir(mymod) if x.startswith("riskAreaAlgorithm_")])))

	return estimatorCLS(**kwargs)


class riskAreaAlgorithm_Sweep(object):
	"""Sweep-based risk area algorithm.

	Estimates the number of casualties at each point on a regular grid
	by sweeping effect isopleths over a demographic layer.  Supports
	parallel execution via ``multiprocessing``.

	Attributes
	----------
	_dxdy : float or None
	    Grid spacing (in metres) for the sweep mesh.
	_outlayers : int or None
	    Number of extra grid cells to add around the bounding box.
	_workerCount : int or None
	    Number of parallel workers to use.
	_runParallel : bool or None
	    Whether to execute the sweep in parallel.
	"""

	_dxdy = None # The density of points to plot.
	_outlayers = None # The number of outlayers to do around the object.
	_workerCount = None
	_runParallel = None

	@property
	def workerCount(self):
		"""
		int
		    The number of parallel worker processes used for the sweep.
		"""
		return self._workerCount

	@workerCount.setter
	def workerCount(self,value):
		"""Set the number of parallel worker processes."""
		self._workerCount = int(value)

	@property
	def parallel(self):
		"""
		bool
		    Whether the sweep calculation runs in parallel using
		    ``multiprocessing``.
		"""
		return self._runParallel

	@parallel.setter
	def parallel(self,value):
		"""Set whether to run the sweep in parallel."""
		self._runParallel = bool(value)

	@property
	def outlayers(self):
		"""
		int
		    The number of extra grid cells added around the bounding box
		    in each direction to avoid boundary effects.
		"""
		return self._outlayers

	@property
	def dxdy(self):
		"""
		float
		    The grid spacing in metres between sweep points.
		"""
		return self._dxdy

	@dxdy.setter
	def dxdy(self,value):
		"""Set the grid spacing in metres."""
		self._dxdy = float(value)

	def __init__(self,dxdy=150,outlayers=3,parallel=True):
		"""
			Initialize the plots. 
			Use 50m as a default. 
		"""
		self._dxdy = dxdy
		self._outlayers = outlayers
		self._workerCount = multiprocessing.cpu_count()
		self.parallel = parallel


	def _findBoundingBox(self,effectIsopleths, demog, mathematical_angle, geometryColumn="TotalPolygon",severityColumn="severity"):
		"""Compute a grid of candidate release points covering the area of interest.

		The bounding box is calculated by rotating the demographic region
		so that the wind direction aligns with mathematical angle 0, then
		expanding it by the width/length of the largest effect isopleth
		plus ``outlayers * dxdy`` in each direction.  The resulting grid of
		points is rotated back to the original coordinate frame.

		Parameters
		----------
		effectIsopleths : thresholdGeoDataFrame
		    The effect isopleths containing threshold polygons and severity
		    information.
		demog : geopandas.GeoDataFrame
		    The demographic data layer used to determine the spatial extent.
		mathematical_angle : float
		    Wind direction in mathematical (counter-clockwise from east)
		    degrees.
		geometryColumn : str, optional
		    Name of the geometry column in ``effectIsopleths`` that holds the
		    total polygons.  Default is ``'TotalPolygon'``.
		severityColumn : str, optional
		    Name of the column that identifies severity levels.  Default is
		    ``'severity'``.

		Returns
		-------
		geopandas.GeoDataFrame
		    A GeoDataFrame of grid points covering the bounding area, rotated
		    back to the original coordinate system.
		"""
		united = demog.convex_hull.unary_union
		minX, minY, maxX, maxY = affinity.rotate(united, origin=united.centroid, angle=-mathematical_angle).bounds

		maxdatetime = effectIsopleths.datetime.max()
		bounds_effectIsopleths = effectIsopleths.set_index("datetime").loc[maxdatetime]

		severityForMesh 	 = bounds_effectIsopleths.set_geometry(geometryColumn).dissolve(by=severityColumn).area.sort_values(ascending=False).index[0]
		effectIsopleths_mesh = bounds_effectIsopleths.query("%s=='%s'" % (severityColumn, severityForMesh))

		minXi, minYi, maxXi, maxYi = effectIsopleths_mesh.set_geometry(geometryColumn).unary_union.bounds

		width = maxYi - minYi
		length = maxXi  # Assume the min is 0.

		minX -= (length + self.outlayers * self.dxdy)
		maxX += self.outlayers * self.dxdy
		minY -= (width + self.outlayers * self.dxdy)
		maxY += width + self.outlayers * self.dxdy

		xCoords = numpy.arange(minX, maxX, self.dxdy)
		yCoords = numpy.arange(minY, maxY, self.dxdy)

		L = []
		for xx, yy in product(xCoords, yCoords): L.append(point.Point(xx, yy))
		ret = geopandas.GeoDataFrame({"points": L}, geometry="points")
		ret = ret.rotate(angle=mathematical_angle, origin=united.centroid)

		return ret

	def _doCalculation(self,releaseLoc,params):
		"""Project effect isopleths at a single release location and return casualties."""
		effectIsopleths = params["effectIsopleths"]
		projected = effectIsopleths.datalayer(params["demog"], releaseLoc, mathematical_angle=params["rotate_angle"])
		if projected is not None:
			data = projected.groupby(["severity", "datetime"]).sum().reset_index()[['severity', 'datetime', params["valueColumn"]]]
			data['x'] = releaseLoc[0]
			data['y'] = releaseLoc[1]
		else:
			L = []
			for severity,timeslice in product(effectIsopleths['severity'].unique(),effectIsopleths['datetime'].unique()):
				L.append(dict(x=releaseLoc[0],y=releaseLoc[1],effectedPopulation=0,severity=severity,datetime=timeslice))
			data = pandas.DataFrame(L)
		return [data]

	def calculate(self,effectIsopleths,demog,mathematical_angle=None,meteorological_angle=None,severityColumn="severity",valueColumn="effectedPopulation"):
		"""	
			Calculates the number of casualties for each point in the find bounding box. 
		"""
		rotate_angle 	= mathematical_angle if meteorological_angle is None else toMathematicalAngle(meteorological_angle)
		severity_effectIsopleths = effectIsopleths
		pointsList = self._findBoundingBox(severity_effectIsopleths,demog,
										   mathematical_angle=rotate_angle,
										   geometryColumn="TotalPolygon",
										   severityColumn=severityColumn)


		params = { "valueColumn" : valueColumn, 
		   	   "demog"    : demog, 
			   "effectIsopleths"  : effectIsopleths, 
			   "rotate_angle" : rotate_angle}

		if self.parallel: 
			F = partial(self._doCalculation,params = params)
			pool = multiprocessing.Pool(self._workerCount)
			listofeffectIsopleths = pool.map(F,[(pointIter.x,pointIter.y) for pointIter in pointsList])
			pool.terminate()
			pool.join()
			resList = chain(*listofeffectIsopleths)
		else:
			resList = []
			for I,pointIter in enumerate(pointsList): 
				print("Processing point %d out of %d (%s)" %(I,len(pointsList),pointIter)) #,end="\r",flush=True)
				data  = self._doCalculation((pointIter.x,pointIter.y),params)
				resList.append(data[0])
			print("")

		return pandas.concat(resList,ignore_index=True,sort=False) 
	
		
		


