import numpy 
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist.floating_axes as floating_axes
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import (FixedLocator, MaxNLocator,
                                                 DictFormatter)

import pandas
from ...utils import toMeteorologicalAngle,toMathematicalAngle, toAzimuthAngle

from descartes import PolygonPatch 

class casualtiesPlot(object): 

	"""
		A class for plotting the different plots related to the casualties. 
	"""
	def plotCasualtiesRose(self,
						   results,
						   area,
						   severityList,
						   loc,
						   meteorological_angles=None,
						   mathematical_angles=None,
						   effectedPopulation="effectedtotal_pop",
						   ax=None,
						   legend=True,
						   weights=None,
						   cycler=None,
						   coordsTickConvertor=toAzimuthAngle,
						   windDistribution=None,
						   windTicks=None,
						   plotType="plot"):
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
		if ((meteorological_angles is None) and (mathematical_angles is None)):
			raise ValueError("Must supply meteorology or mathematical angles")

		rotate_angles 	= mathematical_angles if meteorological_angles is None else [toMathematicalAngle(meteorological_angle) for meteorological_angle in meteorological_angles]
		projectedData = []
		for angle in rotate_angles:
			injuryareas = results.project(area, loc=loc, mathematical_angle=angle)
			injuryareas = injuryareas.groupby("severity")[effectedPopulation].sum().reset_index()
			injuryareas["angle"] = angle
			projectedData.append(injuryareas)
		projectedData = pandas.concat(projectedData)

		pivotedData = projectedData.pivot("angle",'severity',effectedPopulation).reset_index().fillna(0)
		if (ax is None):
			fig = plt.gcf()
			ax = fig.add_subplot(111,polar=True)
			axloc = [111]
		elif isinstance(ax,list):
			fig = plt.gcf()
			axloc = [*ax]
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

			ax.bar(pivotedData["angle"],pivotedData[severity],label=severity,bottom=bottom,**plotprops)
			bottom += pivotedData[severity]
		pivotedData["total"] = 0
		for severity in severityList:
			pivotedData["total"] += pivotedData[severity]
		if windDistribution is not None:
			windDist = windDistribution.copy()
			windTicks = [25,50,75,100] if windTicks is None else sorted(windTicks)
			ax.set_ylim(0,pivotedData["total"].max()*2)
			ax.set_yticks([int(i*pivotedData["total"].max()/3) for i in range(4)])
			maxDist = windTicks[-1]
			windDist["distribution"] = windDist["distribution"] / maxDist + 2.
			windDist["angle"] = 2.5*numpy.pi- windDist["angle"]*numpy.pi/180

			tr = PolarAxes.PolarTransform()

			angle_ticks = [(0, ""),
						   (.25 * numpy.pi, ""),
						   (.5 * numpy.pi, ""),
						   (.75 * numpy.pi, ""),
						   (1. * numpy.pi, ""),
						   (1.25 * numpy.pi, ""),
						   (1.5 * numpy.pi, ""),
						   (1.75 * numpy.pi, "")]

			grid_locator1 = FixedLocator([v for v, s in angle_ticks])
			tick_formatter1 = DictFormatter(dict(angle_ticks))

			radius_ticks = [(2., '')]
			for tick in windTicks[:-1]:
				radius_ticks.append((2+tick/maxDist, '%i %%' % (tick)))


			grid_locator2 = FixedLocator([v for v, s in radius_ticks])
			tick_formatter2 = DictFormatter(dict(radius_ticks))

			grid_helper = floating_axes.GridHelperCurveLinear(tr,
															  extremes=(2. * numpy.pi, 0, 3, 2),
															  grid_locator1=grid_locator1,
															  grid_locator2=grid_locator2,
															  tick_formatter1=tick_formatter1,
															  tick_formatter2=tick_formatter2)

			ax1 = floating_axes.FloatingSubplot(fig, *axloc, grid_helper=grid_helper)
			fig.add_subplot(ax1)

			aux_ax = ax1.get_aux_axes(tr)

			aux_ax.patch = ax1.patch
			ax1.patch.zorder = 0.9

			getattr(aux_ax,plotType)(windDist["angle"], windDist["distribution"])

		# setting meteorology angles.
		metlist = ["$%d^o$" % coordsTickConvertor(x) for x in numpy.linspace(0,360,9)]
		ax.set_xticklabels(metlist)

		if legend:
			plt.legend()
		return ax, pivotedData

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
			raise ValueError("Must supply meteorology or mathematical angle")

		rotate_angle 	= mathematical_angle if meteorological_angle is None else toMathematicalAngle(meteorological_angle)
		retProj		= results.project(area, loc, mathematical_angle=rotate_angle)
		projected  	= retProj.dissolve("severity")

		boundarycycler = plt.cycler(color=plt.rcParams['axes.prop_cycle'].by_key()['color']) if boundarycycler is None else boundarycycler		 		
		cycler = plt.cycler(fc=plt.rcParams['axes.prop_cycle'].by_key()['color'])*plt.cycler(ec=['None']) if cycler is None else cycler		

		patchList = []
		for severity,prop,lineprop in zip(severityList,cycler,boundarycycler):
			if severity not in projected.index:
				continue
			if projected.loc[severity].geometry.type == 'GeometryCollection' or projected.loc[severity].geometry.type == 'MultiPolygon':
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

