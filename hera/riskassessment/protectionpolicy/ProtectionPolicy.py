import numpy 
import json
import xarray 
import pandas
import pydoc 
from unum.units import * 
from ....utils import dosage,tounum,tonumber

class ProtectionPolicy(object): 
	"""
		Calculates the expected concentration of the defined protection policy. 
		The policy should include: 
			- indoor time. 
			- masking. 
			- evacuation. 

		So a policy is a list of actions that take place: 
		each action modifies the concentration field and is defined by a begining and end. 
		The end can be infite (=None), like in evacuation. 

		The object should support chaining of actions: 

		ProtectionPolicy.indoor(alpha=0.2,begin=0,end=10).masking(factor=1000,begin=0,end=1000).compute(C,concentrationField="C")
				  			         .masking(distribution="beta",mean=1000,sigma=10,begin=0,end=1000,bins=100).compute(C,concentrationField="C")
	
		This will create attributes with the details of the execution pipeline. 
		in each attribute we will have: 
				actionID (sequential number) =>  { "params" : { paramname => value }
                                                                   "outputfields" : ["field"] 
        							 } 

		outputfields is used for the subsequent protection (if there is multiple output). 

		Currently, only a single indoor time is defined (without exit time) 

		In the future, we should calculate the concentrations field for each action imposed on the concentrations of the previous 
 		one (i.e add a field). 
	"""
	_data = None
	_xname = None
	_yname = None
	_datetimename = None

	_actionList = None

	_finalname = None

	_params = None

	@property 
	def params(self): 
		return self._params

	@property 
	def finalname(self): 
		return self._finalname

	@property
	def xname(self):
		return self._xname
	@property 
	def yname(self): 
		return self._yname

	@property
	def datetimename(self): 
		return self._datetimename

	@property 
	def data(self):
		return self._data

	def __init__(self,actionList=[],x="x",y="y",datetime="datetime"): 
		"""
			A basic action list. 

			[
				{
					name: "indoor|masking|evacuation",
					params: { action params }
				}
			]
		"""
		self._xname 		= x
		self._yname 		= y 
		self._datetimename 	= datetime
		self._actionList 	= []
		self._finalname  	= "C"
		self.addActions(dict(actions=numpy.atleast_1d(actionList)))


	def addActions(self,jsonStrOrFile): 
		"""
			add actions from a file. 

			The actions are under the "actions" key. 
		"""
		if isinstance(jsonStrOrFile,str): 
			if os.path.exists(jsonStrOrFile):
				with open(jsonStrOrFile,"r") as jsonFile: 
					jsonFile = json.load(jsonFile)
			else: 
				jsonFile = json.load(jsonStrOrFile)
		else:
			jsonFile = jsonStrOrFile

		for actionJSON in jsonFile["actions"]:
			self.addAction(actionJSON["name"],actionJSON["params"])


	def addAction(self,name,params): 
		"""
			Adds json with the params
		"""
		newaction = abstractAction.getAction(len(self._actionList)+1,self,name,params)
		self._actionList.append(newaction)
		return newaction

	def indoor(self,**kwargs): 
		"""
			Adds an indoor protection to the list of actions. 

			Must supply alpha or turnover. If indoor is restricted by time, 
			then must supply (begin,end) or (enter,stay). 

			:params: kwargs: 
			
			* alpha (1/[time])  - room alpha = 1/[turn over rate].
			* turnover ([time]) - room turn over rate. 
			* begin - the time in which the population enters indoor
			* end   - the time in which the population leaves indoor. 
			* enter - timedelta str after the begining of the simulation that the population enters indoor  
			* stay  - timedelta str for staying indoor. 
		"""
		if "alpha" not in kwargs and "turnover" not in kwargs: 
			raise ValueError("must supply either alpha or turnover")

		self.addAction("Indoor",kwargs)
		return self

	def masks(self, **kwargs):
		"""
			Adds a masking protection. 

			:params: kwargs: 
			
			* protectionFactor  - The protectionFactor of the masks. 
			* begin - the time in which the population enters indoor.
			* end   - the time in which the population leaves indoor. 
			* enter - timedelta str after the begining of the simulation that the population enters indoor.
			* stay  - timedelta str for staying indoor. 
		"""

		if "protectionFactor" not in kwargs: 
			raise ValueError("must supply protectionFactor")

		self.addAction("Masks",kwargs)

		return self
		
	def compute(self,data,C="C"): 
		"""
			Executes the pipeline. 
		"""
		self._data = xarray.Dataset()
		self._data.attrs = data.attrs
		self._data["outdoor_0"] = data[C]
		self._data[self.finalname] = data[C]
		self._data.attrs[0] = { "type" : "outdoor"}
		for action in self._actionList: 
			action.compute()

		self._data.compute()
		return self.data

	@property
	def hdfkey(self):
		return "/".join([action.hdfkey for action in self._actionList])



class abstractAction(object): 
	_actionID = None
	_actiontype = None
	_policy = None
	_params = None

	@property 
	def policy(self): 
		return self._policy

	@property 
	def actionid(self): 	
		return self._actionID

	@property 
	def actiontype(self): 	
		return self._actiontype

	@property 
	def params(self): 	
		return self._params

	@property 	
	def name(self,previdef=None): 
		return "%s_%s" % (self._actiontype,self._actionID) 

	def __init__(self,actionID,actiontype,policy,**kwargs): 
		"""
		actionID is a global identifier of the 
		"""
		self._actionID 	 = actionID
		self._actiontype = actiontype
		self._policy 	 = policy
		self._params     = kwargs

	@classmethod 
	def getAction(cls,actionID,policy,name,params): 
		actionCLS = pydoc.locate("pyriskassessment.protectionpolicy.ProtectionPolicy.Action%s" % name.title())
		return actionCLS(actionID,policy,**params)


	@property
	def hdfkey(self):
		"""
			Returns a str that can be used as an HDF key for the current
			policy.

		:return:
		"""
		pass

class ActionIndoor(abstractAction): 

	_alpha = None
	
	@property 
	def turnover(self):
		return 1/self._alpha 

	@property 
	def alpha(self): 
		return self._alpha
		
	def __init__(self,actionID,policy,**kwargs): 
		"""
			the parameters are: 

			alpha - the 1/turnover rate. default unit [1/h] 
			turnover rate -  turn over rate default unit [h]

			begin - the time in which the population enters indoor
			end   - the time in which the population leaves indoor. 

			enter - timedelta str after the begining of the simulation that the population enters indoor  
			stay  - timedelta str for staying indoor. 
		"""
		super().__init__(actionID,"indoor",policy,**kwargs) 
		if "turnover" in kwargs: 
			self._alpha = tounum(1/kwargs["turnover"],1/h) 
		elif "alpha" in kwargs: 
			self._alpha = tounum(kwargs["alpha"],1/h) 
		else: 
			raise ValueError("Must supply either alpha or turnover rate to indoor calculation (ID %s) " % self.actionid)
		


	def compute(self):
		"""
			Compute the indoor concentrations (Cin). 
		"""
		data     = self.policy.data
		Cin      = xarray.zeros_like(data[self.policy.finalname]).squeeze().compute()
		alphanum = tonumber(self.alpha,1/s) 
		Cout     = data[self.policy.finalname].squeeze() 
		
		for I in range(1,len(self.policy.data[self.policy.datetimename])): 
			curstep = {self.policy.datetimename : I}
			prevstep = {self.policy.datetimename : I-1}
			dt =  (pandas.to_timedelta((data.datetime[curstep] - data.datetime[prevstep]).values)).total_seconds()
			Cin[curstep] = (Cin[prevstep] + alphanum*dt*Cout[curstep])/(1+alphanum*dt)

		data[self.name] = Cin

		## Setting the final values as Cout when outside and Cin between begin and end. 
		if "enter" in self.params: 
			if "stay" not in self.params: 
				raise ValueError("Must supply both stay and enter (as timedelta str)")

			abegin = data[self.policy.datetimename][0].values + pandas.to_timedelta(self.params["enter"])
			aend   = abegin + pandas.to_timedelta(self.params["stay"])

		else:
			abegin = self.params.get("begin",None) 
			aend   = self.params.get("end"  ,None) 

		abegin = data[self.policy.datetimename].to_series()[0]  if abegin is None else abegin
		aend   = data[self.policy.datetimename].to_series()[-1] if aend is None else aend

		actionTimeList = data.datetime.to_series()[data[self.policy.datetimename].to_series().between(abegin,aend)]
		data[self.policy.finalname] = Cin.where(data[self.policy.datetimename].isin(actionTimeList),Cout)
		data.attrs[self.actionid] = { "type" : self.actiontype,"actionid": self.actionid,"name" : self.name,\
						      "params" : {
								"alpha" : self.alpha,
								"turnover":self.turnover,
								"begin" : abegin,
								"end"   : aend 
								},"outputs" : [self.name]
						       }

	@property
	def hdfkey(self):
		return "indoorT%dmin%s" % (tonumber(self.turnover,min),self._timekey())

	def _timekey(self):
		data = self.policy.data
		if "enter" in self.params:
			enter = self.params["enter"]
			stay  = self.params["stay"]
		else:
			abegin = self.params.get("begin",None)
			aend   = self.params.get("end"  ,None)
			enter = "0min" if abegin is None else pandas.to_timedelta(abegin - data[self.policy.datetimename][0].values)
			stay  = pandas.to_timedelta(data[self.policy.datetimename][-1].values - data[self.policy.datetimename][0].values) if aend is None else\
					pandas.to_timedelta(aend-abegin)

		return "Enter%sStay%s" % (enter,stay)



class ActionMasks(abstractAction): 

	_protectionFactor = None
	
	@property 
	def protectionFactor(self):
		return self._protectionFactor 
	
	def __init__(self,actionID,policy,**kwargs): 
		"""
			the parameters are: 

			alpha - the 1/turnover rate. default unit [1/h] 
			turnover rate -  turn over rate default unit [h]

			begin - the time in which the population wears masks
			end   - the time in which the population removes masks. 

			wear  - timedelta str after the begining of the simulation that the population wear masks  
			duration  - timedelta str for the duration of wearing masks. 
		"""
		super().__init__(actionID,"masks",policy,**kwargs) 
		try:
			self._protectionFactor = float(kwargs["protectionFactor"])
		except KeyError:
			raise ValueError("Must supply protectionFactor to masks calculation (ID %s) " % self.actionid)


	def compute(self):
		"""
			Compute the indoor concentrations (Cin). 
		"""
		data     	= self.policy.data
		Cno_mask     	= data[self.policy.finalname].squeeze()
		data[self.name] = Cno_mask/self._protectionFactor

		## Setting the final values as Cout when outside and Cin between begin and end. 
		if "wear" in self.params: 
			if "duration" not in self.params: 
				raise ValueError("Must supply both wear and duration (as timedelta str)")

			abegin = data[self.policy.datetimename][0].values + pandas.to_timedelta(self.params["wear"])
			aend   = abegin + pandas.to_timedelta(self.params["duration"])

		else:
			abegin = self.params.get("begin",None) 
			aend   = self.params.get("end"  ,None) 

		abegin = data[self.policy.datetimename].to_series()[0]  if abegin is None else abegin
		aend   = data[self.policy.datetimename].to_series()[-1] if aend is None else aend

		actionTimeList = data.datetime.to_series()[data[self.policy.datetimename].to_series().between(abegin,aend)]
		data[self.policy.finalname] = data[self.name].where(data[self.policy.datetimename].isin(actionTimeList),Cno_mask)
		data.attrs[self.actionid] = { "type" : self.actiontype,"actionid": self.actionid,"name" : self.name,\
						      "params" : {
								"protectionFactor" : self._protectionFactor,
								"begin" : abegin,
								"end"   : aend 
								},"outputs" : [self.name]
						       }
	

	@property
	def hdfkey(self):

		return "maskPF%d%s" % (self.protectionFactor,self._timekey())

	def _timekey(self):
		data = self.policy.data
		if "wear" in self.params:
			wear  = self.params["wear"]
			duration  = self.params["duration"]
		else:
			abegin = self.params.get("begin",None)
			aend   = self.params.get("end"  ,None)
			wear = "0min" if abegin is None else pandas.to_timedelta(abegin - data[self.policy.datetimename][0].values)
			duration  = pandas.to_timedelta(data[self.policy.datetimename][-1].values - data[self.policy.datetimename][0].values) if aend is None else\
						pandas.to_timedelta(aend-abegin)

		return "Wear%sDuration%s" % (wear,duration)

