import numpy 
import matplotlib.pyplot as plt 
from ....utils import toMeteorologicalAngle,toMatematicalAngle,toAzimuthAngle
from descartes import PolygonPatch 

class casualtiesPlot(object): 

	"""
		A class for plotting the different plots related to the casualties. 
	"""
	def plotCasualtiesRose(self,
				projectedData,
				severityList,
				ax=None,
				angleColumn="mathematical_angle_rad",
				legend=True,
				weights=None,
				cycler=None,
				coordsTickConvertor=toAzimuthAngle):
		"""
			plots the total valueColumn in a radial bars according to the severity.  

			*param: :data: a pandas like with the columns :severity and [valueColumn],[angleColumn]. 
				:severityList:   The list of severity values to plot. 
				:valueColumn: the value to plot. 
				:angleColumn: the column that holds the angle. 
					      the angle is a mathematical angle. 	
				:weights:     a list or a scalar to determine the width of the columns. 
				:cycler: a plt cycler for plotting properties per severity
		"""
		
		pivotedData = projectedData.dissolve(by=["severity",angleColumn],aggfunc='sum').reset_index().pivot(angleColumn,'severity','effectedPopulation').reset_index().fillna(0)
		if (ax is None): 
			fig = plt.gcf()
			ax = fig.add_subplot(111,polar=True) 
		elif isinstance(ax,list): 
			fig = plt.gcf()
			ax  = fig.add_subplot(*ax,polar=True) 


		if cycler is None: 
			if weights is None: 
				cycler  = plt.cycler(width=[0.18]*len(severityList))
			else: 
				cycler  = plt.cycler(width=[weights]*len(severityList))
		else: 			
			if "width" not in cycler.keys:
				if weights is None: 
					cycler  += plt.cycler(width=[0.18]*len(severityList))
				else: 
					cycler  += plt.cycler(width=[weights]*len(severityList))

		bottom = numpy.zeros(pivotedData.shape[0])
		for severity,plotprops in zip(severityList,cycler):
			if severity not in pivotedData.columns:
				continue

			ax.bar(pivotedData[angleColumn],pivotedData[severity],label=severity,bottom=bottom,**plotprops)
			bottom += pivotedData[severity]
		
		# setting meteorological angles. 
		metlist = ["$%d^o$" % coordsTickConvertor(x) for x in numpy.linspace(0,360,9)]
		ax.set_xticklabels(metlist)

		if legend: 
			plt.legend()
		return ax

	def plotCasualtiesProjection(self,
				     results,
				     area,
				     severityList,
				     loc,
				     meteorological_angle=None,
				     mathematical_angle=None,
				     plumSeverity=[],
				     ax=None,
				     cycler=None,
				     boundarycycler=None): 
		"""
			Plots the projected data isolpeths of the effected population on the map. 

			:results: The concentration/dosage data and polygons. 
			:area:    A dict with the parameters of the loadResource or a geopandas with the contours. 
			:loc:     The location of emission.
			:severityList: List of severities to draw
			:meteorological_angle:  wind direction
			:mathematical_angle:	wind direction
			:valueColumn: The column to paint. 
			:plumSeverity: a list of severties to plot the overall area. 
			:ax: the fig (if does not exist, create). 
			:cycler: a property cycler for the polygons. 
			:boundarycycler: a property cycler for the overall polygons.  
		"""
		if (ax is None): 
			fig = plt.gcf()
			ax = fig.add_subplot(111) 
		elif isinstance(ax,list): 
			fig = plt.gcf()
			ax  = fig.add_subplot(*ax) 
		else: 
			plt.sca(ax)

		if ((meteorological_angle is None) and (mathematical_angle is None)): 
			raise ValueError("Must supply meteorological or mathematical angle") 

		rotate_angle 	= mathematical_angle if meteorological_angle is None else toMatematicalAngle(meteorological_angle)
		retProj		= results.project(area,loc,mathematical_angle=rotate_angle)
		projected  	= retProj.dissolve("severity")

		boundarycycler = plt.cycler(color=plt.rcParams['axes.prop_cycle'].by_key()['color']) if boundarycycler is None else boundarycycler		 		
		cycler = plt.cycler(fc=plt.rcParams['axes.prop_cycle'].by_key()['color'])*plt.cycler(ec=['None']) if cycler is None else cycler		

		patchList = []
		for severity,prop,lineprop in zip(severityList,cycler,boundarycycler):
			if severity not in projected.index:
				continue
			if projected.loc[severity].geometry.type == 'GeometryCollection': 
				for pol in projected.loc[severity].geometry:
					if pol.type == 'LineString':
						ax.plot(*pol.xy,**lineprop)
					else:
						ax.add_patch(PolygonPatch(pol,**prop) )
			else:
				ax.add_patch(PolygonPatch(projected.loc[severity].geometry,**prop) )
		
		for ((severity,severitydata),prop) in zip(results.shiftLocationAndAngle(loc,mathematical_angle=rotate_angle,geometry="TotalPolygon")
								 .query("severity in %s" % numpy.atleast_1d(plumSeverity))
								 .groupby("severity"),
							  boundarycycler):
			plt.plot(*severitydata.TotalPolygon.convex_hull.unary_union.exterior.xy,**prop)

		return ax,retProj

