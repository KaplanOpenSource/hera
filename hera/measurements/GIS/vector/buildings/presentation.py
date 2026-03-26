"""Presentation layer for the Buildings toolkit.

Provides visualisation methods for building footprints, heights,
and urban morphology parameters (lambda_P, lambda_F, hc).
"""

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

from ...utils import WGS84, ITM


class presentation:
    """Presentation layer for the Buildings toolkit.

    Access via ``buildings.presentation``.

    All plot methods accept common styling parameters (``ax``, ``figsize``,
    ``cmap``, ``alpha``, ``edgecolor``, ``linewidth``, ``colorbar``, ``title``,
    ``xlim``, ``ylim``, ``inputCRS``, ``outputCRS``) and return a
    ``matplotlib.axes.Axes`` object.
    """

    def __init__(self, dataLayer, analysisLayer=None):
        """Initialise with references to the parent toolkit and analysis layer.

        Parameters
        ----------
        dataLayer : BuildingsToolkit
            The parent buildings toolkit instance.
        analysisLayer : analysis, optional
            The buildings analysis layer instance.
        """
        self._datalayer = dataLayer
        self._analysis = analysisLayer

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _reproject(data, inputCRS, outputCRS):
        """Reproject a GeoDataFrame if needed."""
        if data.crs is None and inputCRS is not None:
            data = data.set_crs(epsg=inputCRS)
        if outputCRS is not None and data.crs is not None:
            data = data.to_crs(epsg=outputCRS)
        return data

    @staticmethod
    def _setup_axes(ax, figsize, title, xlim, ylim):
        """Create or configure axes."""
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=figsize)
        if title:
            ax.set_title(title)
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        ax.set_aspect("equal")
        return ax

    @staticmethod
    def _add_colorbar(ax, mappable, label):
        """Add a colorbar to the axes."""
        cb = plt.colorbar(mappable, ax=ax, shrink=0.7, pad=0.02)
        if label:
            cb.set_label(label)
        return cb

    # ------------------------------------------------------------------
    # Plot methods
    # ------------------------------------------------------------------

    def plotBuildings(self, data,
                      ax=None, figsize=(12, 10),
                      facecolor="steelblue", edgecolor="0.3", linewidth=0.3,
                      alpha=0.8,
                      title=None, xlim=None, ylim=None,
                      inputCRS=None, outputCRS=ITM):
        """Plot building footprints as filled polygons.

        Example
        -------
        >>> buildings.presentation.plotBuildings(gdf)

        Parameters
        ----------
        data : geopandas.GeoDataFrame
            Building footprints with a ``geometry`` column.
        ax : matplotlib.axes.Axes, optional
            Axes to plot on. If ``None``, creates a new figure.
        figsize : tuple of float
            Figure size ``(width, height)`` in inches. Default: ``(12, 10)``.
        facecolor : str
            Fill color for buildings. Default: ``"steelblue"``.
        edgecolor : str
            Edge color. Default: ``"0.3"``.
        linewidth : float
            Edge width. Default: ``0.3``.
        alpha : float
            Transparency. Default: ``0.8``.
        title : str, optional
            Plot title.
        xlim : tuple of float, optional
            ``(xmin, xmax)`` axis limits.
        ylim : tuple of float, optional
            ``(ymin, ymax)`` axis limits.
        inputCRS : int, optional
            EPSG code of input data CRS.
        outputCRS : int
            EPSG code to reproject to. Default: ITM (2039).

        Returns
        -------
        matplotlib.axes.Axes
        """
        data = self._reproject(data.copy(), inputCRS, outputCRS)
        ax = self._setup_axes(ax, figsize, title, xlim, ylim)
        data.plot(ax=ax, facecolor=facecolor, edgecolor=edgecolor,
                  linewidth=linewidth, alpha=alpha)
        return ax

    def plotBuildingHeights(self, data, heightColumn="BLDG_HT",
                            ax=None, figsize=(12, 10),
                            cmap="YlOrBr", vmin=None, vmax=None,
                            alpha=0.9, edgecolor="0.4", linewidth=0.2,
                            colorbar=True, colorbar_label="Building height (m)",
                            title=None, xlim=None, ylim=None,
                            inputCRS=None, outputCRS=ITM):
        """Plot building footprints colored by height.

        Example
        -------
        >>> buildings.presentation.plotBuildingHeights(gdf, heightColumn="BLDG_HT")

        Parameters
        ----------
        data : geopandas.GeoDataFrame
            Building footprints with a height column.
        heightColumn : str
            Column name containing building height. Default: ``"BLDG_HT"``.
        ax : matplotlib.axes.Axes, optional
            Axes to plot on.
        figsize : tuple of float
            Figure size. Default: ``(12, 10)``.
        cmap : str or Colormap
            Colormap. Default: ``"YlOrBr"``.
        vmin : float, optional
            Min color scale value.
        vmax : float, optional
            Max color scale value.
        alpha : float
            Transparency. Default: ``0.9``.
        edgecolor : str
            Edge color. Default: ``"0.4"``.
        linewidth : float
            Edge width. Default: ``0.2``.
        colorbar : bool
            Show colorbar. Default: ``True``.
        colorbar_label : str
            Colorbar label. Default: ``"Building height (m)"``.
        title : str, optional
            Plot title.
        xlim : tuple of float, optional
            ``(xmin, xmax)`` axis limits.
        ylim : tuple of float, optional
            ``(ymin, ymax)`` axis limits.
        inputCRS : int, optional
            EPSG code of input data CRS.
        outputCRS : int
            EPSG code to reproject to. Default: ITM (2039).

        Returns
        -------
        matplotlib.axes.Axes
        """
        data = self._reproject(data.copy(), inputCRS, outputCRS)
        ax = self._setup_axes(ax, figsize, title, xlim, ylim)
        mappable = data.plot(ax=ax, column=heightColumn, cmap=cmap,
                             vmin=vmin, vmax=vmax, alpha=alpha,
                             edgecolor=edgecolor, linewidth=linewidth,
                             legend=False)
        if colorbar:
            sm = plt.cm.ScalarMappable(
                cmap=cmap,
                norm=mcolors.Normalize(
                    vmin=vmin if vmin is not None else data[heightColumn].min(),
                    vmax=vmax if vmax is not None else data[heightColumn].max(),
                ),
            )
            sm._A = []
            self._add_colorbar(ax, sm, colorbar_label)
        return ax

    def plotMorphology(self, morphologyData, column="lambdaP",
                       ax=None, figsize=(12, 10),
                       cmap="RdYlGn_r", vmin=None, vmax=None,
                       alpha=0.8, edgecolor="0.3", linewidth=0.5,
                       colorbar=True, colorbar_label=None,
                       title=None, xlim=None, ylim=None,
                       inputCRS=None, outputCRS=ITM):
        """Plot urban morphology parameters as a block grid choropleth.

        Use this to visualise the output of
        ``buildings.analysis.LambdaFromBuildingData()``.

        Example
        -------
        >>> morph = buildings.analysis.LambdaFromBuildingData(270, 250, gdf)
        >>> buildings.presentation.plotMorphology(morph, column="lambdaP")
        >>> buildings.presentation.plotMorphology(morph, column="lambdaF")
        >>> buildings.presentation.plotMorphology(morph, column="hc")

        Parameters
        ----------
        morphologyData : geopandas.GeoDataFrame
            Output from ``LambdaFromBuildingData`` with columns
            ``lambdaP``, ``lambdaF``, ``hc``, ``geometry``.
        column : str
            Column to plot. One of ``"lambdaP"``, ``"lambdaF"``, ``"hc"``.
            Default: ``"lambdaP"``.
        ax : matplotlib.axes.Axes, optional
            Axes to plot on.
        figsize : tuple of float
            Figure size. Default: ``(12, 10)``.
        cmap : str or Colormap
            Colormap. Default: ``"RdYlGn_r"`` (red=high density, green=low).
        vmin : float, optional
            Min color scale value.
        vmax : float, optional
            Max color scale value.
        alpha : float
            Transparency. Default: ``0.8``.
        edgecolor : str
            Block edge color. Default: ``"0.3"``.
        linewidth : float
            Block edge width. Default: ``0.5``.
        colorbar : bool
            Show colorbar. Default: ``True``.
        colorbar_label : str, optional
            Colorbar label. If ``None``, uses the column name.
        title : str, optional
            Plot title.
        xlim : tuple of float, optional
            ``(xmin, xmax)`` axis limits.
        ylim : tuple of float, optional
            ``(ymin, ymax)`` axis limits.
        inputCRS : int, optional
            EPSG code of input data CRS.
        outputCRS : int
            EPSG code to reproject to. Default: ITM (2039).

        Returns
        -------
        matplotlib.axes.Axes
        """
        data = self._reproject(morphologyData.copy(), inputCRS, outputCRS)
        ax = self._setup_axes(ax, figsize, title, xlim, ylim)

        label = colorbar_label if colorbar_label else column
        data.plot(ax=ax, column=column, cmap=cmap,
                  vmin=vmin, vmax=vmax, alpha=alpha,
                  edgecolor=edgecolor, linewidth=linewidth,
                  legend=False)
        if colorbar:
            sm = plt.cm.ScalarMappable(
                cmap=cmap,
                norm=mcolors.Normalize(
                    vmin=vmin if vmin is not None else data[column].min(),
                    vmax=vmax if vmax is not None else data[column].max(),
                ),
            )
            sm._A = []
            self._add_colorbar(ax, sm, label)
        return ax

    def plotBuildingsWithMorphology(self, buildingsData, morphologyData,
                                    morphColumn="lambdaP",
                                    ax=None, figsize=(14, 10),
                                    buildingColor="0.2", buildingAlpha=0.6,
                                    buildingEdgecolor="0.1", buildingLinewidth=0.2,
                                    cmap="RdYlGn_r", vmin=None, vmax=None,
                                    morphAlpha=0.5, morphEdgecolor="0.3",
                                    morphLinewidth=0.5,
                                    colorbar=True, colorbar_label=None,
                                    title=None, xlim=None, ylim=None,
                                    inputCRS=None, outputCRS=ITM):
        """Overlay building footprints on a morphology block grid.

        Example
        -------
        >>> buildings.presentation.plotBuildingsWithMorphology(
        ...     gdf, morph, morphColumn="lambdaP")

        Parameters
        ----------
        buildingsData : geopandas.GeoDataFrame
            Building footprints.
        morphologyData : geopandas.GeoDataFrame
            Output from ``LambdaFromBuildingData``.
        morphColumn : str
            Morphology column to color blocks by. Default: ``"lambdaP"``.
        ax : matplotlib.axes.Axes, optional
            Axes to plot on.
        figsize : tuple of float
            Figure size. Default: ``(14, 10)``.
        buildingColor : str
            Fill color for buildings. Default: ``"0.2"`` (dark gray).
        buildingAlpha : float
            Building transparency. Default: ``0.6``.
        buildingEdgecolor : str
            Building edge color. Default: ``"0.1"``.
        buildingLinewidth : float
            Building edge width. Default: ``0.2``.
        cmap : str or Colormap
            Colormap for morphology blocks. Default: ``"RdYlGn_r"``.
        vmin : float, optional
            Min color scale value.
        vmax : float, optional
            Max color scale value.
        morphAlpha : float
            Block transparency. Default: ``0.5``.
        morphEdgecolor : str
            Block edge color. Default: ``"0.3"``.
        morphLinewidth : float
            Block edge width. Default: ``0.5``.
        colorbar : bool
            Show colorbar. Default: ``True``.
        colorbar_label : str, optional
            Colorbar label.
        title : str, optional
            Plot title.
        xlim : tuple of float, optional
            ``(xmin, xmax)`` axis limits.
        ylim : tuple of float, optional
            ``(ymin, ymax)`` axis limits.
        inputCRS : int, optional
            EPSG code of input data CRS.
        outputCRS : int
            EPSG code to reproject to. Default: ITM (2039).

        Returns
        -------
        matplotlib.axes.Axes
        """
        # Plot morphology blocks first (background)
        ax = self.plotMorphology(
            morphologyData, column=morphColumn,
            ax=ax, figsize=figsize,
            cmap=cmap, vmin=vmin, vmax=vmax,
            alpha=morphAlpha, edgecolor=morphEdgecolor,
            linewidth=morphLinewidth,
            colorbar=colorbar, colorbar_label=colorbar_label,
            title=title, xlim=xlim, ylim=ylim,
            inputCRS=inputCRS, outputCRS=outputCRS,
        )
        # Overlay buildings
        bld = self._reproject(buildingsData.copy(), inputCRS, outputCRS)
        bld.plot(ax=ax, facecolor=buildingColor, edgecolor=buildingEdgecolor,
                 linewidth=buildingLinewidth, alpha=buildingAlpha)
        return ax

    def plotHeightHistogram(self, data, heightColumn="BLDG_HT",
                            ax=None, figsize=(10, 6),
                            bins=30, color="steelblue", edgecolor="white",
                            title=None, xlabel="Building height (m)",
                            ylabel="Count"):
        """Plot a histogram of building heights.

        Example
        -------
        >>> buildings.presentation.plotHeightHistogram(gdf)

        Parameters
        ----------
        data : geopandas.GeoDataFrame
            Building data with a height column.
        heightColumn : str
            Column name for height. Default: ``"BLDG_HT"``.
        ax : matplotlib.axes.Axes, optional
            Axes to plot on.
        figsize : tuple of float
            Figure size. Default: ``(10, 6)``.
        bins : int
            Number of histogram bins. Default: ``30``.
        color : str
            Bar color. Default: ``"steelblue"``.
        edgecolor : str
            Bar edge color. Default: ``"white"``.
        title : str, optional
            Plot title.
        xlabel : str
            X-axis label. Default: ``"Building height (m)"``.
        ylabel : str
            Y-axis label. Default: ``"Count"``.

        Returns
        -------
        matplotlib.axes.Axes
        """
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=figsize)

        heights = data[heightColumn].dropna()
        heights = heights[heights > 0]
        ax.hist(heights, bins=bins, color=color, edgecolor=edgecolor)
        if title:
            ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        return ax

    def plotBuildingsOnMap(self, data, tilesToolkit,
                           heightColumn=None,
                           ax=None, figsize=(14, 10),
                           cmap="YlOrBr", vmin=None, vmax=None,
                           alpha=0.7, edgecolor="0.3", linewidth=0.2,
                           colorbar=True, colorbar_label="Building height (m)",
                           zoomlevel=15, tileServer=None,
                           title=None, xlim=None, ylim=None,
                           inputCRS=None, outputCRS=ITM):
        """Overlay building footprints on a tile server map.

        Example
        -------
        >>> tiles = toolkitHome.getToolkit(toolkitHome.GIS_TILES, projectName="MY_PROJECT")
        >>> buildings.presentation.plotBuildingsOnMap(gdf, tiles, heightColumn="BLDG_HT")

        Parameters
        ----------
        data : geopandas.GeoDataFrame
            Building footprints.
        tilesToolkit : TilesToolkit
            Initialised Tiles toolkit for fetching map images.
        heightColumn : str, optional
            If provided, color buildings by this column. If ``None``,
            plots all buildings in a uniform color.
        ax : matplotlib.axes.Axes, optional
            Axes to plot on.
        figsize : tuple of float
            Figure size. Default: ``(14, 10)``.
        cmap : str or Colormap
            Colormap (used when *heightColumn* is set). Default: ``"YlOrBr"``.
        vmin : float, optional
            Min color scale value.
        vmax : float, optional
            Max color scale value.
        alpha : float
            Building transparency. Default: ``0.7``.
        edgecolor : str
            Edge color. Default: ``"0.3"``.
        linewidth : float
            Edge width. Default: ``0.2``.
        colorbar : bool
            Show colorbar (only when *heightColumn* is set). Default: ``True``.
        colorbar_label : str
            Colorbar label. Default: ``"Building height (m)"``.
        zoomlevel : int
            Tile server zoom level. Default: ``15``.
        tileServer : str, optional
            Tile server URL or datasource name.
        title : str, optional
            Plot title.
        xlim : tuple of float, optional
            ``(xmin, xmax)`` axis limits.
        ylim : tuple of float, optional
            ``(ymin, ymax)`` axis limits.
        inputCRS : int, optional
            EPSG code of input data CRS.
        outputCRS : int
            EPSG code to reproject to. Default: ITM (2039).

        Returns
        -------
        matplotlib.axes.Axes
        """
        data = self._reproject(data.copy(), inputCRS, outputCRS)

        # Compute extent for tile fetching
        bounds = data.total_bounds
        extent = [bounds[0], bounds[1], bounds[2], bounds[3]]

        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=figsize)

        # Render map tiles as background
        tilesToolkit.presentation.plot(
            extent=extent, ax=ax,
            zoomlevel=zoomlevel, tileServer=tileServer,
        )

        # Overlay buildings
        if heightColumn and heightColumn in data.columns:
            data.plot(ax=ax, column=heightColumn, cmap=cmap,
                      vmin=vmin, vmax=vmax, alpha=alpha,
                      edgecolor=edgecolor, linewidth=linewidth,
                      legend=False)
            if colorbar:
                sm = plt.cm.ScalarMappable(
                    cmap=cmap,
                    norm=mcolors.Normalize(
                        vmin=vmin if vmin is not None else data[heightColumn].min(),
                        vmax=vmax if vmax is not None else data[heightColumn].max(),
                    ),
                )
                sm._A = []
                self._add_colorbar(ax, sm, colorbar_label)
        else:
            data.plot(ax=ax, facecolor="steelblue", edgecolor=edgecolor,
                      linewidth=linewidth, alpha=alpha)

        if title:
            ax.set_title(title)
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        ax.set_aspect("equal")
        return ax
