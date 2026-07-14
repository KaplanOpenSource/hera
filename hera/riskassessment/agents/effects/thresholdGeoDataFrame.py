import numpy
import collections
import pandas
import geopandas

from hera.utils import toMeteorologicalAngle, toMathematicalAngle
from hera.measurements.GIS import ITM


import geopandas as gpd



class thresholdGeoDataFrame(geopandas.GeoDataFrame):
	"""A GeoDataFrame that holds threshold-based injury isopleths.

	Extends ``geopandas.GeoDataFrame`` with convenience methods for
	shifting, rotating, and projecting threshold polygons onto
	demographic datasets.
	"""

	def __init__(self, *args, **kwargs):
		"""Initialize a thresholdGeoDataFrame.

		Parameters
		----------
		*args
		    Positional arguments forwarded to ``geopandas.GeoDataFrame``.
		**kwargs
		    Keyword arguments forwarded to ``geopandas.GeoDataFrame``.
		"""
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
		"""Rotate and translate polygons to a given location and wind angle."""
		self = self.set_geometry(geometry)

		rotate_angle = mathematical_angle if meteorological_angle is None else toMathematicalAngle(meteorological_angle)

		if (rotate_angle is None):
			raise ValueError("either met_angle or math_angle should be provided") 			

		shiftedPolygons = self.rotate(rotate_angle,origin=(0,0,0)).translate(*loc)   
		return shiftedPolygons


	def project(self,demographic,loc,meteorological_angle=None,mathematical_angle=None,geometry="ThresholdPolygon",population="total_pop",projectName=None):
		"""
		Projects the results on a demographic data set.
		The rotation and translation are in respect to the (0,0) coordinate of the concentration field.

		Parameters
		----------
		demographic : geopandas.GeoDataFrame
			The demographic data.
		loc : tuple
			The (x, y) location of the release.
		meteorological_angle : float or list of float, optional
			Wind angle(s) in meteorological convention.
		mathematical_angle : float or list of float, optional
			Wind angle(s) in mathematical convention.
		geometry : str
			The geometry column name.
		population : str or list of str
			The population column(s) to use.
		projectName : str, optional
			The project name for the DemographyToolkit. If None, uses the
			project associated with the demographic data source.
		"""
		if isinstance(meteorological_angle,collections.Iterable):
			retList = []
			for metangle in meteorological_angle:
				projectedValue = self._project(demographic=demographic,loc=loc,meteorological_angle=metangle,geometry=geometry,projectName=projectName)
				if projectedValue is None:
					continue
				projectedValue["meteorological_angle"] 		= metangle
				projectedValue["mathematical_angle_rad"] 	= numpy.deg2rad(toMathematicalAngle(metangle))
				projectedValue["mathematical_angle"] 		= toMathematicalAngle(metangle)
				retList.append(projectedValue)
			ret = pandas.concat(retList)

		elif isinstance(mathematical_angle,collections.Iterable):
			retList = []
			for mathangle in mathematical_angle:
				projectedValue=self._project(demographic=demographic,loc=loc,mathematical_angle=mathangle,geometry=geometry,projectName=projectName)
				if projectedValue is None:
					continue
				projectedValue["meteorological_angle"] 		= toMeteorologicalAngle(mathangle)
				projectedValue["mathematical_angle_rad"] 	= numpy.deg2rad(mathangle)
				projectedValue["mathematical_angle"] 		= mathangle
				retList.append(projectedValue)

			ret = pandas.concat(retList)
		else:
			ret = self._project(demographic=demographic,loc=loc,meteorological_angle=meteorological_angle,
								mathematical_angle=mathematical_angle,geometry=geometry, population=population,projectName=projectName)
		return ret 

	@staticmethod
	def _calculatePopulationInPolygon(demography, poly, populationTypes):
		"""
		Calculate the population within a polygon by spatial intersection.

		Parameters
		----------
		demography : geopandas.GeoDataFrame
			The demographic data with geometry and population columns.
		poly : shapely.geometry
			The polygon to intersect with.
		populationTypes : list of str
			Column names for population counts.

		Returns
		-------
		geopandas.GeoDataFrame
			The intersection with area-weighted population.
		"""
		res_intersect_poly = demography.loc[demography["geometry"].intersection(poly).is_empty == False]
		intersection_poly = res_intersect_poly["geometry"].intersection(poly)

		res_intersection = geopandas.GeoDataFrame.from_dict(
			{"geometry": intersection_poly.geometry,
			 "areaFraction": intersection_poly.area / res_intersect_poly.area})

		for populationType in populationTypes:
			if populationType in res_intersect_poly:
				res_intersection[populationType] = intersection_poly.area / res_intersect_poly.area * res_intersect_poly[populationType]

		return res_intersection

	def _project(self,demographic,loc,meteorological_angle=None,mathematical_angle=None,geometry="ThresholdPolygon",population="total_pop",projectName=None):
		"""
		Projects the polygons of the thresholds on the demography.
		Shifts and rotates the polygons according to loc and the angle of the wind.

		Parameters
		----------
		demographic : geopandas.GeoDataFrame
			The demographic data.
		loc : tuple
			The (x, y) location of the release.
		meteorological_angle : float
			The angle of the wind in meteorological angles.
			Only meteorological_angle or mathematical_angle should be supplied.
		mathematical_angle : float
			Only meteorological_angle or mathematical_angle should be supplied.
		geometry : str
			The name of the column to use for the thresholds.
		population : str or list of str
			The name of the column(s) to use for the population.
		projectName : str, optional
			Unused. Kept for backward compatibility.

		Returns
		-------
		pandas.DataFrame or None
			Casualty estimates per severity and time step, or None if no
			population is affected.
		"""
		localcrs = ITM

		demog_data = demographic
		demog_data = demog_data.to_crs(localcrs) # convert to itm. It is in m**2.

		shiftedPolygons = self._shiftPolygons(loc=loc,meteorological_angle=meteorological_angle,mathematical_angle=mathematical_angle,geometry=geometry)

		retList = []
		population = [population] if type(population)==str else population
		for ((severity,timestamp),data) in self.groupby(["severity","datetime"]):
			for indx,row in data.iterrows():
				curpoly = shiftedPolygons.loc[indx]
				if curpoly.is_valid:
					res = self._calculatePopulationInPolygon(demog_data, curpoly, population)
					if len(res) > 0:
						for popu in population:
							res['effected%s' % popu] = res[popu]*row.percentEffected
						res['percentEffected']  = row.percentEffected
						res['ToxicLoad']        = row.ToxicLoad
						res['severity'] = severity
						res['datetime'] = timestamp
						retList.append(res)

		return None if len(retList) ==0 else pandas.concat(retList)
