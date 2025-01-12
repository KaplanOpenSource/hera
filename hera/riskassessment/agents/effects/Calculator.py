import json
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


	def toJSON(self):
		return dict(breathingRate=str(self.injuryBreathingRate))

	def __str__(self):
		return json.dumps(self.toJSON())

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

	def calculate(self,concentrationField,breathingRate=10*L/min,time="datetime",field=None,inUnits=None):
		"""
			Calculates the dosage from a concentration field.

			Concentration field can be pandas.DataFrame or xarray.Dataset.
			If it is xarray.Dataset, then use the field parameter to calculate the toxic from the requested field.

			Parameters
			----------
				concentrationField : pandas.DataFrame or xarray dataframe.
						The concentrations. Must have 'time' index.

						For pandas, we do not assume that the time is equispaced.

				breathingRate : unum, L/min
		    			The breathing rate of the population.

	  			time : str
	  					The name of the time coordinates. In pandas it is the name of the index.
	  					That is, the index will be accessed after a reset_index() function, and that field must be present.

	  			field : str
	  					Only in xarray. The name of the field to use from the xarray (default None).
	  					Ignored in pandas.
	  			inUnits: unum
	  					The units of the concentration.
	  					If None, use the default units.

	  					The default units of pandas are mg/m**3. For xarray it the units in the attributes.
	  		Return
	  		-------
	  				An xarray or pandas (depend on the input) with the toxic loads
	  				in each time step.


		"""
		breathingRatio = (breathingRate/self.injuryBreathingRate).asNumber()
		if inUnits is None:
			if hasattr(concentrationField, "attrs"):
				attrs_units = concentrationField.attrs.get("field", None)
			inUnits = concentrationField.attrs[field] if attrs_units is not None  else mg / m ** 3

		CunitConversion = inUnits.asNumber(mg / m ** 3)

		if isinstance(concentrationField,xarray.Dataset):
			return concentrationField[field].cumsum(dim=time)*concentrationField.dt.asNumber(min)*breathingRatio*CunitConversion
		elif isinstance(concentrationField,pandas.DataFrame):
			dt_min = concentrationField.reset_index()[time].diff().apply(lambda x: x.seconds)/60.
			return (concentrationField[:-1].fillna(0)*dt_min[1:]).cumsum()*breathingRatio * CunitConversion
		else:
			raise ValueError("concentrationField is not a pandas.DataFrame or xarray.Dataset")

	def toJSON(self):
		ret = super().toJSON()
		ret['type'] = 'haber'
		return ret


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

	def calculate(self,concentrationField,field,breathingRate=10*L/min,time="datetime",inUnits=None):
		"""
            Calculates the toxic load  from a concentration field.
            \begin{equation}
            		D(T) = \int_0^T C^n dt
            \end{equation}

            Concentration field can be pandas.DataFrame or xarray.Dataset.
            If it is xarray.Dataset, then use the field parameter to calculate the toxic from the requested field.

            Parameters
            ----------
                concentrationField : pandas.DataFrame or xarray dataframe.
                        The concentrations. Must have 'time' index.

                        For pandas, we do not assume that the time is equispaced.

                breathingRate : unum, L/min
                        The breathing rate of the population.

                  time : str
                          The name of the time coordinates. In pandas it is the name of the index.
                          That is, the index will be accessed after a reset_index() function, and that field must be present.

                  field : str
                          Only in xarray. The name of the field to use from the xarray (default None).
                          Ignored in pandas.
                  inUnits: unum
                          The units of the concentration.
                          If None, use the default units.

                          The default units of pandas are mg/m**3. For xarray it the units in the attributes.
              Return
              -------
                      An xarray or pandas (depend on the input) with the toxic loads
                      in each time step.


        """
		breathingRatio = (breathingRate/self.injuryBreathingRate).asNumber()

		if inUnits is None:
			if hasattr(concentrationField, "attrs"):
				attrs_units = concentrationField.attrs.get("field", None)
			inUnits = concentrationField.attrs[field] if attrs_units is not None  else mg / m ** 3
		CunitConversion = inUnits.asNumber(mg / m ** 3)

		if isinstance(concentrationField, xarray.Dataset):
			return ((concentrationField[field]*CunitConversion)**self.n).cumsum(dim=time)*concentrationField.dt.asNumber(min)*breathingRatio
		elif isinstance(concentrationField,pandas.DataFrame):
			dt_min = concentrationField.reset_index()[time].diff().apply(lambda x: x.seconds) / 60.
			C_to_n = ((concentrationField[:-1].fillna(0)* CunitConversion )**self.n)
			# we have to ignore the index, so we transform one of the vectors to numpy.
			return (C_to_n * dt_min[1:].values.reshape(C_to_n.shape[0],1)).cumsum() * breathingRatio
		else:
			raise ValueError("concentrationField is not a pandas.DataFrame or xarray.Dataset")

	def toJSON(self):
		ret = super().toJSON()
		ret['type'] = 'tenBerge'
		ret['n'] = str(self.n)
		return ret

class CalculatorMaxConcentration(AbstractCalculator): 
	
	def __init__(self,sampling,breathingRate=10*L/min,**kwargs): 
		"""
			Implements the maximal concentration on the time axis.
		"""
		super().__init__(breathingRate=breathingRate)
		self._sampling = sampling

	def calculate(self,concentrationField,field,breathingRate=10*L/min,time="datetime",inUnits=None):
		"""
            Calculate the maximal concentrations

            Concentration field can be pandas.DataFrame or xarray.Dataset.
            If it is xarray.Dataset, then use the field parameter to calculate the toxic from the requested field.

            Parameters
            ----------
                concentrationField : pandas.DataFrame or xarray dataframe.
                        The concentrations. Must have 'time' index.

                        For pandas, we do not assume that the time is equispaced.

                breathingRate : unum, L/min
                        The breathing rate of the population.

                  time : str
                          The name of the time coordinates. In pandas it is the name of the index.
                          That is, the index will be accessed after a reset_index() function, and that field must be present.

                  field : str
                          Only in xarray. The name of the field to use from the xarray (default None).
                          Ignored in pandas.
                  inUnits: unum
                          The units of the concentration.
                          If None, use the default units.

                          The default units of pandas are mg/m**3. For xarray it the units in the attributes.
              Return
              -------
                      An xarray or pandas (depend on the input) with the toxic loads
                      in each time step.


        """
		breathingRatio = (breathingRate / self.injuryBreathingRate).asNumber()

		if inUnits is None:
			if hasattr(concentrationField, "attrs"):
				attrs_units = concentrationField.attrs.get("field", None)
			inUnits = attrs_units[field] if attrs_units is not None  else mg / m ** 3
		CunitConversion = inUnits.asNumber(mg / m ** 3)

		if isinstance(concentrationField, xarray.Dataset):
			itemstep = pandas.to_timedelta(self._sampling) / pandas.to_timedelta(
				str(concentrationField.attrs["dt"]).replace("[", "").replace("]", ""))
			samplingparam = {time : int(itemstep)}
			return concentrationField[field].chunk(chunks={time: int(3*itemstep)}).rolling(**samplingparam).mean().fillna(0).max(dim=time)*breathingRatio*CunitConversion
		elif isinstance(concentrationField,pandas.DataFrame):
			return concentrationField.rolling(window=self._sampling).mean().fillna(0).max()*breathingRatio*CunitConversion
		else:
			raise ValueError("concentrationField is not a pandas.DataFrame or xarray.Dataset")


	def toJSON(self):
		ret = super().toJSON()
		ret['type'] = 'maxConcentrations'
		return ret