import pandas
import json
import pydoc
from unum.units import *
from .thresholdGeoDataFrame import thresholdGeoDataFrame 

class InjuryFactory(object): 

	_name = None 

	@property 
	def name(self): 
		return self._name 

	def getInjury(self,name,cfgJSON,**additionalparameters): 
		"""
			Creates the appropriate injury 

			Structure: 
			{
				"type" : <name> . Appends 'pyriskassessment.agents.effects.InjuryLevel.Injury' to the name.
				"calculator" :{ <calculatorName> : {"time" : "datetime"} }  Appends 'pyriskassessment.agents.effects.Calculator.Calculator' to the name.
				"parameters": { 
					... Injury level json. 
				}
			}

			additionalparameters: 
			=====================
				These are parameters that might be used by the calculator

				For example the tenberge coefficient. 
		"""		
		try:
	                injuryType = cfgJSON["type"]
		except KeyError:
                        raise ValueError("Injury type is not defined")
                
		try:
			calculatorTypeAndParams = cfgJSON["calculator"]
			calcType,calcParam       = [x for x in calculatorTypeAndParams.items()][0]
		except KeyError: 
                        raise ValueError("Calculator not defined")

		calculatorCLS = pydoc.locate("hera.riskassessment.agents.effects.Calculator.Calculator%s" % calcType)
		calculator    = calculatorCLS(**calcParam,**additionalparameters)

		injuryCLS = pydoc.locate("hera.riskassessment.agents.effects.Injury.Injury%s" % injuryType)
		
		if injuryCLS is None: 
			
			injuryExists  = ",".join([x[6:] for x in dir(pydoc.locate("hera.riskassessment.agents.effects.Injury")) if x.startswith("Injury")])

			raise NotImplementedError("The injury %s is not defined. Injuries found: %s  " % (injuryType,injuryExists))

		return injuryCLS(name,cfgJSON["parameters"],calculator,units=cfgJSON.get("units"))

class Injury(object): 
	"""
                Holds a list of injury levels 
	"""

	_name = None
	_levelsmap = None   # just a map of name->injury.
	_levels = None 
	_calculator = None

	@property
	def calculator(self):
	    return self._calculator

	@property 
	def levels(self): 
		return self._levels

	@property
	def levelNames(self): 
		return [x.name for x in self._levels]

	def __getitem__(self,name):
		return self._levelsmap[name]

	def __init__(self,name,cfgJSON,calculator,units=None):
		"""
			Loads the relevant injury levels from the JSONcfg

			{

				"type" : <injuryLevel type>.  Appends 'pyriskassessment.casualties.InjuryLevel.InjuryLevel' to the name. 
				"levels" : [....],
			        parameters: { 
					     "InjuryName" : { injury parameters },
						.
						.			
				
			}

		"""
		self._name = name
		self._calculator = calculator

		injuryType = cfgJSON.get("type",None)
		if (injuryType is None):
			raise ValueError("InjuryLevel type is nor defined")

		injuryCLS = pydoc.locate("hera.riskassessment.agents.effects.InjuryLevel.InjuryLevel%s" % injuryType)
		self._levelsmap = {}
		self._levels = []
		levelNames = cfgJSON.get("levels")

		if units is not None:
			units = eval(units)

		for lvl in levelNames:
			try:
				lvlparams = cfgJSON["parameters"][lvl].copy()
			except KeyError:
				raise ValueError("parameters for level %s not found !" % lvl)

			curindex = levelNames.index(lvl)
			lvlparams["higher_severity"] = None if curindex == 0 else self.levels[curindex-1]

			injry = injuryCLS(lvl, units=units, **lvlparams)
			self._levels.append(injry)
			self._levelsmap[lvl] = injry



	def getPercent(self,level, ToxicLoad):
		"""
			calculate the toxic load and subtract the percent from the percent above. 
		"""
		severityList = self.levelNames
		curindex     = severityList.index(level)

		if (curindex == 0): 
			val = self.levels[curindex].getPercent(ToxicLoad) 
		else: 
			val = self.levels[curindex].getPercent(ToxicLoad)  - self.levels[curindex-1].getPercent(ToxicLoad)
		return val

	def _postCalculate(self,retList,time):
		"""
			Apply some post calculations on the results.

			The post calculation depends on the injury type, and therefore implemented in the inherited classes.

			Parameters
			----------
				retList : list of pandas.
						The values of the calculation in different levels.

				time : str
						The name of the time columns.
			
		"""
		raise NotImplementedError("Abstract class")

	def _postCalculatePointWise(self, retList):
		"""
			apply some post calculations on the results

		"""
		pass

	def calculate(self, concentrationField, field, time="datetime", x="x", y="y", breathingRate=10 * L / min,
				  **parameters):
		import warnings
		warnings.warn("This function is obselete. Use calculateRegionOfInjured instead")
		return self.calculateRegionOfInjured(concentrationField=concentrationField,
											  field=field,
											  time=time,
											  x=x,
											  y=y,
											  breathingRate=breathingRate,
											  **parameters)

	def calculateRegionOfInjured(self,concentrationField,field,time="datetime",x="x",y="y",breathingRate=10*L/min,**parameters):
		"""
			Calculates the fraction of the population that was effected in each point, and returns its contour.
			The levels are determined by the injury type.

			Parameters
			----------
			concentrationField:
				an xarray with the concentrations at each time.

			time: str
				The name time dimension. if None then don't look up for time (assume it is without time).

			x:  str
				The name x dimension.

			y:  str
				The name of the y dimension.

			breathingRate:
					The breathing rate used. The default is man at Rest.

			parameters: kwargs
					Additional parameters that are needed for the calculations.


			Returns
			--------

            	A thresholdGeoDataFrame with the columns:

				time |  injury name | total polygon  | diff polygon                             |       addtional fields
					 |              | the total area | a differance from the level above.       |       if necessary.
		"""
                
		retList = [] 

		toxicLoads = self.calculateToxicLoads(concentrationField=concentrationField,
											  time=time,
											  breathingRate=breathingRate,
											  field=field,
											  **parameters).compute()
		for lvl in self.levels:
			data = lvl.calculateContours(toxicLoads=toxicLoads,time=time,x=x,y=y)
			if data is not None: 
				retList.append(data)

		ret = self._postCalculate(retList,time)
		return thresholdGeoDataFrame(ret)

	def calculateToxicLoads(self,concentrationField,time="datetime",breathingRate=10*L/min,field=None):
		"""
			Calculates the toxic loads of the concetration field.

			Paramters
			---------
					concentrationField : pandas or xarray.
							The concentrations to calculate.

							If xarray, then we assume that it is a 2D map with an additional 'time' coordinate.
							If pandas, then we assume that it is a 'time' in the index and each other column is a point to calculate.

					time : str
							The name of the time column (or the name of the index column if it is pandas).

					breathingRate: unum [L/min].
							The breathing rate used. The default is man at Rest is 10L/min.

					field : str
							If xarray, the name of the field to use.
							ignored if pandas

		"""
		return self.calculator.calculate(concentrationField, field, breathingRate=breathingRate, time=time)


	def calculatePointWiseFractionInjured(self,timeConcentration,time="datetime",breathingRate=10*L/min,field=None):
		"""
			Calculates the fraction of injury over time in each point.

			Can be used with pandas.DataFrame or  xarray.Dataset [not implemented yet].

			Parameters
			----------

				timeConcentration : pandas.DataFrame, or xarray.DataFrame.
							Holds the concentration in time.

							If xarray, also has a 'time' coordinate that will be calculated. [not implemented yet]

							If DataFrame:
								Each point is represeneted by a column, and the time is an index (with the name 'datetime').

								So the structure of the input is :
									  P1   P2   P3
								time
								00:00 0     0   0
								00:01 0.1   0.1 0.01
								00:02 0.4   0.1 0.01
								00:05 0.5   0.1 0.01

		    time  : str
		    			The name of the time column (or the name of the index).

		    breathingRate : unum, L/min
		    			The breathing rate of the population.

		    parameters: kwargs.
		    			Additional parameter to the calculator (for example ten-berge coefficient).

						selection parameters (xarray only):

							sel - select according to the coordinates (see sel funcion of the xarray).
							isel - select according to the coordinate index (see isel function of the xarray).

		"""




		# 2. For each injury level:
		#      Create a dataframe with the fields:
		#
		#					P1        injury
		#		datetime
		#		  0  	[fraction]        level 1 name
		#		  1  	[fraction]        level 1 name
		#					...

		if not isinstance(timeConcentration,pandas.DataFrame):
			raise ValueError("Still not implemented....")

		# 1. Calculate the toxic load for each point.
		toxicLoads = self.calculateToxicLoads(concentrationField=timeConcentration,
											  time=time,
											  breathingRate=breathingRate,
											  field=field)
		retList = []
		for lvl in self.levels:
			data = pandas.DataFrame()

			for device in toxicLoads:
				prct = toxicLoads[device].apply(lambda x: self.getPercent(lvl.name,x))
				data = prct.to_frame("injuryPercent").assign(deviceName=device,level=lvl.name)
				if data is not None:
					retList.append(data)

		ret = pandas.concat(retList)
		return ret


	def calculateThresholdPolygon(self,data,time):
		"""
			Calculates the diff of the polygon based on the toxic load. 
		"""
		ret = data.sort_values("ToxicLoad")
		polyList = []
		indexList = []
		for timeseries in ret.groupby(time):
			timedata = timeseries[1].sort_values("ToxicLoad",ascending=False)
			polyList.append(timedata.iloc[0]["TotalPolygon"])
			indexList.append(timedata.index[0]) 
			for curPolyIndex,prevPolyIndex in zip(timedata.index[1:],timedata.index[:-1]): 
				curPoly = ret.loc[curPolyIndex,"TotalPolygon"] 
				prevPoly = ret.loc[prevPolyIndex,"TotalPolygon"] 
				polyList.append(curPoly.difference(prevPoly))
				indexList.append(curPolyIndex)

		diff = pandas.DataFrame({"ThresholdPolygon" : polyList},index=indexList) 
		ret = data.merge(diff,left_index=True,right_index=True)		
				
		return ret 


	def toJSON(self):
		ret = dict()
		#ret['name']  = self._name
		ret['levels'] = dict([(lvl.name,lvl.toJSON()) for lvl in self.levels])
		ret['calculator'] = self._calculator.toJSON()
		return ret

	def __str__(self):
		json.dumps(self.toJSON(),indent=4)




class InjuryLognormal10(Injury): 
	def _postCalculate(self,retList,time):
		"""
                        Fill in the actual % that was effected. 

		"""
		modList = []
		for vals in retList: 
			vals["percentEffected"] = vals.apply(lambda x: self.getPercent(x['severity'],x['ToxicLoad']),axis=1)
			vals = self.calculateThresholdPolygon(vals,time)
			modList.append(vals)

		return pandas.concat(modList,ignore_index=True)


class InjuryThreshold(Injury): 
	def _postCalculate(self,retList,time):
		"""
                        Calculate the percent Effected and calculate the differential polygons. 
		"""
		if len(retList) > 0: 
			ret = pandas.concat(retList,ignore_index=True)
			ret = self.calculateThresholdPolygon(ret,time)
			ret['percentEffected'] = 1.
		else: 
			ret = None
		return ret 


class InjuryExponential(Injury):
	def _postCalculate(self,retList,time):
		"""
                        Calculate the percent Effected and calculate the differential polygons. 
		"""
		modList = []
		for vals in retList: 
			vals["percentEffected"] = vals.apply(lambda x: self.getPercent(x['severity'],x['ToxicLoad']),axis=1)
			vals = self.calculateThresholdPolygon(vals,time)
			modList.append(vals)

		return pandas.concat(modList,ignore_index=True)


