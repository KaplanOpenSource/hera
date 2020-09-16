import pandas
import geopandas 
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

		calculatorCLS = pydoc.locate("hera.risk.riskassessment.agents.effects.Calculator.Calculator%s" % calcType)
		calculator    = calculatorCLS(**calcParam,**additionalparameters)

		injuryCLS = pydoc.locate("hera.risk.riskassessment.agents.effects.Injury.Injury%s" % injuryType)
		
		if injuryCLS is None: 
			
			injuryExists  = ",".join([x[6:] for x in dir(pydoc.locate("hera.risk.riskassessment.agents.effects.Injury")) if x.startswith("Injury")])

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

		injuryCLS = pydoc.locate("hera.risk.riskassessment.agents.effects.InjuryLevel.InjuryLevel%s" % injuryType)
		self._levelsmap = {}
		self._levels = []
		levelNames = cfgJSON.get("levels")

		if units is not None:
			units = eval(units)

		for lvl in levelNames:
			try:
				lvlparams = cfgJSON["parameters"][lvl]
			except KeyError:
				raise ValueError("parameters for level %s not found !" % lvl)

			curindex = levelNames.index(lvl)
			lvlparams["higher_severity"] = None if curindex == 0 else self.levels[curindex-1]
			injry = injuryCLS(lvl, calculator=calculator, units=units, **lvlparams)
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

		

	def _postCalculate(self,retList): 
		"""
			apply some post calculations on the results 
			
		"""
		pass

	def calculate(self,concentrationField,field,time="datetime",x="x",y="y",breathingRate=10*L/min,**parameters):
		"""
                        Calculates the number of people in the isopleths of the concentrationField. 
                        The object itself determines which levels to draw. 

                        :param: concentrationField: an xarray with the concentrations at each time. 

                        :param: time: The time dimension. if None then don't look up for time (assume it is without time). 
			:param: x: the x dimension. 
			:param: y: the y dimension. 
			:param: parameters: Additional parameters that are needed for the calculations. 
			:param: breathingRate: The breathing rate used. The default is man at Rest. 
                        :return 
                                A pandas object with the columns: 
                                time |  injury name | total polygon  | diff polygon                             |       addtional fields 
                                     |              | the total area | a differance from the level above.       |       if necessary. 
		"""
                
		retList = [] 

		for lvl in self.levels:
			data = lvl.calculate(concentrationField=concentrationField,
					     field=field,
					     time=time,
					     x=x,
					     y=y,
					     breathingRate=breathingRate,
 					     **parameters)
			if data is not None: 
				retList.append(data)

		ret = self._postCalculate(retList,time) 
		return thresholdGeoDataFrame(ret)

	def calculateRaw(self,concentrationField,field,time="datetime",x="x",y="y",breathingRate=10*L/min,**parameters):
		data = self.levels[0].calculateRaw(concentrationField=concentrationField,
											field=field,
											time=time,
											x=x,
											y=y,
											breathingRate=breathingRate,
											**parameters)
		return data


		
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


class ExponentialThreshold(Injury): 

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


