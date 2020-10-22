import pandas
import xarray 
from unum.units import *

class AbstractCalculator(object): 
	"""
		Holds the abstract calculator. 
	"""

	# The breathing rate with which the injury level was determined. 
	_injury_breathingRate = None 
	@property 
	def injuryBreathingRate(self): 
		return self._injury_breathingRate

	def __init__(self,breathingRate=10*L/min): 
		self._injury_breathingRate=breathingRate

class CalculatorHaber(AbstractCalculator): 
	"""
		Implements the haber law (simple integral) on the concentrations. 
	"""

	def __init__(self,breathingRate=10*L/min,**kwargs): 
		"""
			Calculates the Haber dosage. The integration is performed on the time axis. 

			:param: time: a string with the name of the axis to perform the integration on. 
		"""
		super().__init__(breathingRate=breathingRate)

	def calculate(self,concentrationField,field,breathingRate=10*L/min,time="datetime"): 
		breathingRatio = (breathingRate/self.injuryBreathingRate).asNumber()
		CunitConversion = concentrationField.attrs[field].asNumber(mg/m**3)
		return concentrationField[field].cumsum(dim=time)*concentrationField.dt.asNumber(min)*breathingRatio*CunitConversion


class CalculatorTenBerge(AbstractCalculator): 
	"""
		Implements the ten=berge law (int C^n dt) on the concentrations.
	"""

	n = None

	def __init__(self,tenbergeCoefficient,breathingRate=10*L/min,**kwargs): 
		"""
			:param: n: the ten-berge coefficient. 
			:param: time: the time variable. 
		"""
		super().__init__(breathingRate=breathingRate)
		self.n 	   = tenbergeCoefficient

	def calculate(self,concentrationField,field,breathingRate=10*L/min,time="datetime"): 
		breathingRatio = (breathingRate/self.injuryBreathingRate).asNumber()
		CunitConversion = concentrationField.attrs[field].asNumber(mg/m**3)
		return ((concentrationField[field]*CunitConversion)**self.n).cumsum(dim=time)*concentrationField.dt.asNumber(min)*breathingRatio


	def calculateRaw(self,pandasField,dt,time="datetime",C="C",breathingRate=10*L/min,inUnit=mg/m**3):
		breathingRatio = (breathingRate/self.injuryBreathingRate).asNumber()
		CunitConversion = inUnit.asNumber(mg / m ** 3)
		pfield = pandasField.set_index(time).sort_index()
		TL     = ((pfield[C]* CunitConversion) ** self.n).cumsum() * dt.asNumber(min) * breathingRatio
		return pandas.merge(pandasField,TL.to_frame("ToxicLoad"),left_on=time,right_index=True)



class CalculatorMaxConcentration(AbstractCalculator): 
	
	def __init__(self,sampling,breathingRate=10*L/min,**kwargs): 
		"""
			Implements the maximal concentration on the time axis.
		"""
		super().__init__(breathingRate=breathingRate)
		self._sampling = sampling

	def calculate(self,concentrationField,field,breathingRate=10*L/min,time="datetime"): 
		breathingRatio = (breathingRate/self.injuryBreathingRate).asNumber()
		itemstep = pandas.to_timedelta(self._sampling)/pandas.to_timedelta(str(concentrationField.attrs["dt"]).replace("[","").replace("]",""))
		samplingparam = {time : int(itemstep)}
		CunitConversion = concentrationField.attrs[field].asNumber(mg/m**3)
		return concentrationField[field].chunk(chunks={time: int(3*itemstep)}).rolling(**samplingparam).mean().max(dim=time)*breathingRatio*CunitConversion

