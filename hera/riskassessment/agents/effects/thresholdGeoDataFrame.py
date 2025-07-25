import numpy
import collections
import pandas
import geopandas

from ....utils import toMeteorologicalAngle,toMathematicalAngle


from hera.measurements.GIS.vector.demography import DemographyToolkit as demoDatalayer

pop = demoDatalayer(projectName="Demography") # used only for the analysis layer, which is static.



class thresholdGeoDataFrame(geopandas.GeoDataFrame): 

	def __init__(self,*args,**kwargs): 
	
		super(thresholdGeoDataFrame,self).__init__(*args,**kwargs)

	def shiftLocationAndAngle(self,loc,meteorological_angle=None,mathematical_angle=None,geometry="ThresholdPolygon"):
		"""
			returns a new thresholdGeoDataFrame with a shifted geometry polygons. 
		"""
		ret = self.copy()
		ret[geometry] = self._shiftPolygons(loc=loc,meteorological_angle=meteorological_angle,mathematical_angle=mathematical_angle,geometry=geometry)
		ret = ret.set_geometry(geometry)
		return ret 

	def _shiftPolygons(self,loc,meteorological_angle=None,mathematical_angle=None,geometry="ThresholdPolygon"):
		self = self.set_geometry(geometry)

		rotate_angle = mathematical_angle if meteorological_angle is None else toMathematicalAngle(meteorological_angle)

		if (rotate_angle is None):
			raise ValueError("either met_angle or math_angle should be provided") 			

		shiftedPolygons = self.rotate(rotate_angle,origin=(0,0,0)).translate(*loc)   
		return shiftedPolygons


	def project(self,demographic,loc,meteorological_angle=None,mathematical_angle=None,geometry="ThresholdPolygon",population="total_pop"):
		"""
			Projects the results on a demographic data set. 
			The rotation and translation are in respect to the (0,0) coordinate of the concentration field. 

			:param: demographic - either a geopandas dataframe or a dictionary 
                        that represetents the parameters that will pass to the GISResource. 

			:param: loc - the location of the fall. 
			:param: meteorological_angle - wind  in meteorology angle.
							either this or the mathematical_angle
			:param: mathematical_angle   - wind  in meteorology angle.
							either this or the meteorological_angle 
		"""
		if isinstance(meteorological_angle,collections.Iterable): 
			retList = []
			radList = []
			for metangle in meteorological_angle: 
				projectedValue = self._project(demographic=demographic,loc=loc,meteorological_angle=metangle,geometry=geometry)
				if projectedValue is None: 
					continue
				projectedValue["meteorological_angle"] 		= metangle
				projectedValue["mathematical_angle_rad"] 	= numpy.deg2rad(toMathematicalAngle(metangle))
				projectedValue["mathematical_angle"] 		= toMathematicalAngle(metangle)
				radList.append(projectedValue)
			ret = pandas.concat(retList) 

		elif isinstance(mathematical_angle,collections.Iterable): 
			retList = []
			for mathangle in mathematical_angle: 
				projectedValue=self._project(demographic=demographic,loc=loc,mathematical_angle=mathangle,geometry=geometry)
				if projectedValue is None: 
					continue
				projectedValue["meteorological_angle"] 		= toMeteorologicalAngle(mathangle)
				projectedValue["mathematical_angle_rad"] 	= numpy.deg2rad(mathangle)
				projectedValue["mathematical_angle"] 		= mathangle
				retList.append(projectedValue)

			ret = pandas.concat(retList) 
		else:
			ret = self._project(demographic=demographic,loc=loc,meteorological_angle=meteorological_angle,
								mathematical_angle=mathematical_angle,geometry=geometry, population=population)
		return ret 

	def _project(self,demographic,loc,meteorological_angle=None,mathematical_angle=None,geometry="ThresholdPolygon",population="total_pop"):
		"""
			Projects the polygons of the thresholds on the demography.
			Shifts and rotates the polygons according to loc and the angle of the wind.

		Parameters
		----------
		demographic : geopandnaas
				Holds the deomgraphy
		loc : tuple
			The location of the thresholds.
		meteorological_angle : float
			The angle of the wind in meteorological angles.
			Only meteorological_angle or mathematical_angle should be supplied.

		mathematical_angle
			Only meteorological_angle or mathematical_angle should be supplied.
		geometry : string
			The name of the column to use for the thresholds.

		population : string
			The name of the column to use for the population.

		Returns
		-------

		"""
		localcrs = {"init":"epsg:2039"} # itm

		demog_data = demographic
		demog_data = demog_data.to_crs(localcrs) # convert to itm. It is in m**2. 

		shiftedPolygons = self._shiftPolygons(loc=loc,meteorological_angle=meteorological_angle,mathematical_angle=mathematical_angle,geometry=geometry)

		retList = []
		population = [population] if type(population)==str else population
		for ((severity,timestamp),data) in self.groupby(["severity","datetime"]):
			for indx,row in data.iterrows():
				curpoly = shiftedPolygons.loc[indx]
				if curpoly.is_valid:
					res     = pop.analysis.calculatePopulationInPolygon(dataSourceOrData=demog_data,shapeNameOrData=curpoly, populationTypes=population)
					if len(res) > 0:
						for popu in population:
							res['effected%s' % popu] = res[popu]*row.percentEffected
						res['percentEffected']  = row.percentEffected
						res['ToxicLoad']        = row.ToxicLoad
						res['severity'] = severity
						res['datetime'] = timestamp
						retList.append(res)

		return None if len(retList) ==0 else pandas.concat(retList)
