import geopandas
import pandas
import numpy
import json

from hera.utils.matplotlibCountour import toGeopandas
from ....utils import tounit

from unum.units import *
from unum import Unum

import matplotlib.pyplot as plt 
from scipy.stats import lognorm

dosage = mg/m**3

class InjuryLevel(object): 
	"""
		Represents an injury level 
	"""
	_name =None

	_parameters = None 

	_calculator = None

	_units 	    = None 
	
	@property 
	def name(self): 
		return self._name 


	@property 
	def units(self):
		return self._units

	def __init__(self,name,units=None,**parameters):
		self._name = name
		self._parameters = parameters
		self._units = units


	def calculateContours(self, toxicLoads, time="datetime", x="x", y="y"):
		"""
			Calculates the contours of the current injury from the toxic loads.

			Parameters
			----------
				toxicLoads : xarray.Dataset
						The 2D map of the toxicloads in time.


				time : str
						The name of the time column

				x : str
						The name of the x column.

				y : str
						The name of the y column.



			Return
			-------
			GeoPandas.GeoDataFrame with the fields:

				time |  injury name | total polygon  | diff polygon 				| percentInjury
           			 |       	    | the total area | a differance from the level above. 	| if necessary.
			
	
			The data returns in km because itm is in km. 

		"""
		curBackend = plt.get_backend()	
		plt.switch_backend("pdf")

		if time in toxicLoads.dims:
			if toxicLoads[time].size == 1:
				timeList = [pandas.to_datetime(toxicLoads[time].item())]
				timesel  = False 
			else: 
				timeList = toxicLoads[time]
				timesel = True
		elif hasattr(toxicLoads,time):
			timeList = [pandas.to_datetime(getattr(toxicLoads,time).item())]
			timesel  = False
		else: 
			timeList= [0] 
			timesel = False 

		retList = []

		for tt in timeList: 
			concField = toxicLoads.sel(**{time:tt}) if timesel else toxicLoads

			ret = self._getGeopandas(concField,x,y) #,**parameters)
			if not ret.empty:

				ret = ret[ret.area > 1]   # 1 m**2. 
				ret['contour'] = ret['contour'].buffer(0) # sometimes fixes invalid polygons.

				ret = ret.dissolve(by="Level").reset_index()

				if "cloud_math_direction" in toxicLoads.attrs:
					ret['contour'] = ret['contour'].rotate(-toxicLoads.attrs["cloud_math_direction"], origin=(0, 0))

				ret["severity"] = self.name
				ret = ret.rename(columns={"Level" : "ToxicLoad","contour":"TotalPolygon",time : "datetime"}) 
				ret = ret.set_geometry("TotalPolygon")
				ret[time] = pandas.to_datetime(tt.item()) if timesel else tt
				retList.append(ret)
		
		ret = None if ret.empty else pandas.concat(retList,ignore_index=True)
		plt.switch_backend(curBackend)
		
		return ret


	def _getGeopandas(self,concentrationField,x,y,**parameters):
		"""
			Return the correct geopandas of the  Injury level 
		"""
		raise NotImplementedError("abstract function")



	def toJSON(self):
		return dict(name=self.name,
					units=str(self.units))


	def __str__(self):
		json.dumps(self.toJSON(),indent=4)




# ******************************************************

class InjuryLevelLognormal10DoseResponse(InjuryLevel): 
	"""
		Calculates the % of the population that will be affected at a given ToxicLoad 
		
		This class also return the toxic load of a requested fraction. 
	"""

	@property 
	def TL_50(self): 	
		return self._parameters["TL_50"]
	
	@property 
	def sigma(self):
		return self._parameters["sigma"]
	
	def __init__(self,name,**parameters):
		"""
			Initiates the injury level. 
			
			parameters are 
				- TL_50: The ToxicLoad in which 50% are effected. 
				- sigma: The probit inverse. 
		"""
		if "TL_50" not in parameters: 
			raise ValueError("Must supply TL_50") 
		else: 
			parameters["TL_50"] = tounit(parameters["TL_50"],dosage)
		if "sigma" not in parameters: 
			raise ValueError("Must supply sigma (1/probit)") 

		super(InjuryLevelLognormal10DoseResponse,self).__init__(name,**parameters)

	def getPercent(self,ToxicLoad):
		"""
	        Using lognormal distribution with base 10.

	        In that case define a new probable variable:
	
	        Y such that Y = X/a. (X is the original dose in base 10, Y is in natural base).
	        Then, X_avg/a =  Y_avg and  Sig_x/a = Sig_y.
	
	        Where a is the transformation between the bases:
	        10^(x) = exp(y) --> x = log10(e)*y --> a = log10(e).

	        Therefore, in order to use the lognormal with the natural base we divide
	        all the data (and the mean,std) with the 'a' coefficient.

	        Important Notes:
	        ---------------
	        1. See that the a factor in v/a and Effect[0]/a actually cancel.
	        2. lognorm calculates log(x/scale)/sigma --> log(x/scale)=  log10(x/scale)/a
	           But we use sigma/a ==> and we get that lognorm.cdf(x/a,sigma/a,scale/a) = lognorm.cdf(x,sigma/a,scale) = Phi(log10(x/scale)/a / (sigma/a) )
	                                  = Phi(log10(x/scale)/sigma ) which is what we got from the manuals (and see Ziv).

		"""
		a = numpy.log10(numpy.e)
		prob = lognorm.cdf(ToxicLoad/a,self.sigma/a,scale=self.TL_50.asNumber(dosage)/a)
		return prob


	def getToxicLoad(self,Percent): 
		a = numpy.log10(numpy.e)
		percent = lognorm.ppf(Percent,self.sigma/a,scale=self.TL_50.asNumber(dosage)/a)*a
		return percent

	def _getGeopandas(self,concentrationField,x,y,**parameters):
		"""
			Return the correct geopandas of the  Injury level
			higher_severity
		"""
		defaultLevels = numpy.arange(0.05, 1, 0.05)
		percentList = parameters.get("percentInjury",defaultLevels)
		ToxicLoads = self.getToxicLoad(percentList)

		## We add the levels of higher toxicit y to solve the problem
		# of a wide gape between the two severities (which led to double counting the casualties). 
		if self._parameters.get("higher_severity",None) is not None:
			HigherToxicLoads = self._parameters["higher_severity"].getToxicLoad(defaultLevels)
			ToxicLoads = numpy.unique(numpy.sort(numpy.concatenate([ToxicLoads,HigherToxicLoads])))

		CS =  plt.contour(concentrationField[x],concentrationField[y],concentrationField.squeeze(),levels=ToxicLoads)

		# Changes by ofir... not clear yet when the original fails....
		#concentrationField = concentrationField.to_dataframe().reset_index()
		#CS = plt.tricontour(concentrationField[x], concentrationField[y], concentrationField, levels=ToxicLoads) This line has bug because
		#	concentrationField is panas and we have to select the right column.

		if numpy.max(CS.levels) < numpy.min(ToxicLoads): 
			ret = geopandas.GeoDataFrame()
		else:
			ret = toGeopandas(CS)
		return ret

	def toJSON(self):
		ret = super().toJSON()
		ret['type'] = 'logNormalBase10'
		ret['TL_50'] = str(self.TL_50)
		ret['sigma'] = str(self.sigma)
		return ret


################################################################################################

class InjuryLevelThreshold(InjuryLevel): 

	@property 
	def threshold(self): 
		return self._parameters["threshold"]

	def __init__(self,name,units=None,**parameters):

		actualunit = mg/m**3 if units is None else units
		if "threshold" not in parameters: 
			raise ValueError("Cannot find the threshold")

		thr = eval(parameters["threshold"])

		parameters["threshold"] = thr if isinstance(thr,Unum) else tounit(eval(parameters["threshold"]),actualunit)

		super(InjuryLevelThreshold,self).__init__(name,units=actualunit,**parameters)

	def getPercent(self,ToxicLoad):
		threshold = self.threshold
		ToxicLoad = numpy.atleast_1d(ToxicLoad)
		ret = numpy.array([ 1 if tounit(x,self.units) > threshold else 0 for x in ToxicLoad])
		return ret[0] if len(ToxicLoad)==1 else ret

	def _getGeopandas(self,concentrationField,x,y,**parameters):
		"""
			Return the correct geopandas of the  Injury level 
		"""
		level = self.threshold.asNumber()
		CS =  plt.contour(concentrationField[x],concentrationField[y],concentrationField.squeeze(),levels=numpy.atleast_1d(level))
		if numpy.max(CS.levels) < level: 
			ret = geopandas.GeoDataFrame()
		else:
			ret = toGeopandas(CS)
		return ret

	def toJSON(self):
		ret = super().toJSON()
		ret['type'] = 'threshold'
		ret['threshold'] = str(self.threshold)
		return ret

################################################################################################

class InjuryLevelExponential(InjuryLevel): 

	@property 
	def k(self): 
		return self._parameters["k"]

	def __init__(self,name,**parameters):


		parameters["k"] = float(parameters["k"])

		super().__init__(name,**parameters)

	def getPercent(self,ToxicLoad):
		
		ret = 1-numpy.exp(-self.k*numpy.array(ToxicLoad))
		return ret 

	def _getGeopandas(self,concentrationField,x,y,**parameters):
		"""
			Return the correct geopandas of the  Injury level 
		"""
		level = self.k
		CS =  plt.contour(concentrationField[x],concentrationField[y],concentrationField.squeeze(),levels=level)
		if numpy.max(CS.levels) < level: 
			ret = geopandas.GeoDataFrame()
		else:
			ret = toGeopandas(CS)

		return ret


	def toJSON(self):
		ret = super().toJSON()
		ret['type'] = 'exponential'
		ret['k'] = str(self.k)
		return ret




