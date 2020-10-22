import numpy
import geopandas
import pydoc
from shapely import affinity
from shapely.geometry import GeometryCollection,LineString,point
import pandas
from itertools import product,chain
import multiprocessing
from functools import partial

from ....utils import toMeteorologicalAngle,toMatematicalAngle,toAzimuthAngle

def getRiskAreaAlgorithm(algorithmName,**kwargs):
	"""
		Return an estimator class. 
		kwargs are the parameters passed  as the constructor. 

		currently only Sweep is implemented
	"""
	estimatorCLS = pydoc.locate("pyriskassessment.analysis.riskAreas.riskAreaAlgorithm_%s" % algorithmName.title())
	if estimatorCLS is None: 
		mymod = pydoc.locate("pyriskassessment.analysis.riskAreas")
		raise ValueError("estimator %s not found. Available estimators are: %s" % (algorithmName,",".join([x.split("_")[1] for x in dir(mymod) if x.startswith("riskAreaAlgorithm_")])))

	return estimatorCLS(**kwargs)


class riskAreaAlgorithm_Sweep(object): 

	_dxdy = None # The density of points to plot. 
	_outlayers = None # The number of outlayers to do around the object. 
	_workerCount = None 
	_runParallel = None 

	@property 
	def workerCount(self):
		return self._workerCount 

	@workerCount.setter
	def workerCount(self,value):
		self._workerCount = int(value)

	@property
	def parallel(self):
		return self._runParallel

	@parallel.setter
	def parallel(self,value):
		self._runParallel = bool(value)

	@property 
	def outlayers(self):
		return self._outlayers

	@property
	def dxdy(self): 
		return self._dxdy

	@dxdy.setter
	def dxdy(self,value): 
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
		"""
			Returns a polygon of the area to search.
			The polygon will be rotated such that the wind is in mathematical angle of 0.


			Then it will encompass the demog_polygon + the width/length of the effectPolygons + outlayers*dxdy.

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
		effectIsopleths = params["effectIsopleths"]
		projected = effectIsopleths.project(params["demog"],releaseLoc,mathematical_angle=params["rotate_angle"])
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
		rotate_angle 	= mathematical_angle if meteorological_angle is None else toMatematicalAngle(meteorological_angle)
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
	
		
		


