"""hera/measurements/GIS/vector/topography.py -- the *vector* (height-contour)
topography toolkit and its analysis layer.

Covered here
------------
``TopographyToolkit.__init__``, the ``stlFactory`` property,
``cutRegionFromSource``, ``geoPandasToSTL``, ``regionToSTL``, ``toDEM``, and on
the ``analysis`` class: ``__init__``, the ``datalayer`` property and
``addHeight``.

How the inputs were produced
----------------------------
The toolkit's native input is a *contour* layer: LineString geometries plus a
height column named by ``stlFactory.heightColumnsNames`` (default ``HEIGHT``).
``_contours()`` below builds three horizontal, parallel iso-height lines in ITM
metres, so the interpolated surface between them is a linear ramp in y and the
resulting raster/DEM values are known in closed form rather than from a golden
file.  Nothing here needs a real shapefile, GDAL beyond geopandas, or FreeCAD.

Bugs pinned (each has an xfail(strict) for the intended behaviour plus a
passing characterisation of what actually happens)
--------------------------------------------------------------------------
* **B227** -- ``cutRegionFromSource`` calls ``super().cutRegionFromSource``
  with keyword names that do not exist on the parent
  (``shapeDataOrName=`` / ``datasourceName=`` / ``crs=`` against the parent's
  ``datasourceDocument`` / ``shape`` / ``inputCRS``).  Every call raises
  ``TypeError`` on the very first statement, so the method is unreachable in
  its entirety -- and with it ``regionToSTL``, whose only job is to call it.
* **B228** -- ``geoPandasToSTL`` accepts a ``solidName`` argument and then
  forwards the hard-coded literal ``"Topography"`` instead of it, so the
  caller's solid name is silently discarded.
* **B229** -- ``toDEM`` calls ``self.getDatasourceData`` (lower-case ``s``);
  the toolkit only defines ``getDataSourceData``, so passing a datasource name
  -- the documented first parameter -- raises ``AttributeError``.  The same
  typo makes the geoJSON-string fallback on the next line dead code.
* **B230** -- ``toDEM`` reads ``Nx`` and ``Ny`` from ``grid_x.shape[0]`` /
  ``[1]``, then emits ``rasterized['height'][j, i]`` with ``i`` ranging over
  ``Nx`` and ``j`` over ``Ny`` -- the two subscripts are swapped.  On a
  non-square grid this raises ``IndexError``; on a square one it silently
  transposes the DEM.
* **B231** -- ``analysis.addHeight`` calls
  ``coordinateHandler.regularizeTimeSteps(...)``, but the module-level
  ``from hera.simulations.utils import coordinateHandler`` binds the *module*
  of that name, not the class inside it (``hera/simulations/utils/__init__.py``
  is empty).  Every call raises ``AttributeError``.
* **B232** -- a second, independent defect further down ``addHeight``: the
  save-mode comparisons read ``toolkit.TOOLKIT_SAVEMODE_FILEANDDB`` etc., but
  the module-level ``toolkit`` name was rebound by
  ``from hera.measurements.GIS.vector import toolkit`` to the *vector* toolkit
  module, which defines no ``TOOLKIT_SAVEMODE_*`` constants.  The list
  literal is built before the membership test, so *every* save mode fails --
  ``TOOLKIT_SAVEMODE_NOSAVE`` included -- and ``addHeight`` can never return.
  Reached in the test by neutralising B231 only (the coordinate handler), so
  the second fault is pinned independently of the first.

Deliberately not tested, and why
--------------------------------
* The numeric facet geometry of the STL string: that belongs to
  ``stlFactory.rasterToSTL``, which has its own tests in
  ``test_gis_utils_stlfactory.py``.  Here only the delegation, the solid-name
  handling and the "empty region" branch are checked.
* ``cutRegionFromSource``'s internals past line 1 (the ``getRegionData`` /
  ``getDatasourceDocument`` calls, both also non-existent methods, and the
  CRS-defaulting logic): unreachable while B227 stands, so asserting on them
  would be asserting on dead code.  They are noted rather than pinned.
"""
import os

import geopandas
import numpy
import pandas
import pytest
from shapely.geometry import LineString, Polygon, box

from hera.measurements.GIS.utils import ITM, stlFactory
from hera.measurements.GIS.vector import topography as vector_topography_module
from hera.measurements.GIS.vector.topography import TopographyToolkit, analysis

# --------------------------------------------------------------------------
# synthetic contour input
# --------------------------------------------------------------------------

_X0, _X1 = 200000.0, 200400.0   # ITM easting  span 400 m
_YS = (660000.0, 660200.0, 660400.0)  # three contour lines, 200 m apart
_HEIGHTS = (100.0, 120.0, 140.0)      # -> linear ramp of 0.1 m per m of y


def _contours(x0=_X0, x1=_X1, ys=_YS, heights=_HEIGHTS, crs=ITM):
    """Three horizontal iso-height LineStrings: a linear ramp in y."""
    return geopandas.GeoDataFrame(
        {
            "HEIGHT": list(heights),
            "geometry": [LineString([(x0, y), (x1, y)]) for y in ys],
        },
        crs=crs,
    )


@pytest.fixture()
def topo(unit_toolkit_factory):
    from hera import toolkitHome

    return unit_toolkit_factory(toolkitHome.GIS_VECTOR_TOPOGRAPHY)


@pytest.fixture()
def contours():
    return _contours()


# --------------------------------------------------------------------------
# construction
# --------------------------------------------------------------------------

@pytest.mark.unit
class TestConstruction:
    def test_the_toolkit_name_is_topography(self, topo):
        assert topo.toolkitName == "Topography"

    def test_it_is_a_vector_toolkit(self, topo):
        from hera.measurements.GIS.vector.toolkit import VectorToolkit

        assert isinstance(topo, VectorToolkit)

    def test_the_factory_returns_the_expected_class(self, topo):
        assert isinstance(topo, TopographyToolkit)

    def test_the_analysis_layer_is_built_and_points_back_at_the_toolkit(self, topo):
        assert isinstance(topo._analysis, analysis)
        assert topo._analysis.datalayer is topo

    def test_the_files_directory_is_the_one_that_was_passed(self, topo, unit_files_directory):
        assert os.path.abspath(topo.filesDirectory) == os.path.abspath(unit_files_directory)


@pytest.mark.unit
class TestStlFactoryProperty:
    def test_it_exposes_an_stlfactory_instance(self, topo):
        assert isinstance(topo.stlFactory, stlFactory)

    def test_it_is_the_same_object_on_every_access(self, topo):
        """A property, not a factory function -- state set on it must stick."""
        assert topo.stlFactory is topo.stlFactory

    def test_its_height_column_default_is_visible_through_the_toolkit(self, topo):
        assert topo.stlFactory.heightColumnsNames == "HEIGHT"

    def test_the_height_column_can_be_retargeted_through_the_property(self, topo):
        topo.stlFactory.heightColumnsNames = "Z"
        assert topo.stlFactory.heightColumnsNames == "Z"


@pytest.mark.unit
class TestAnalysisLayer:
    def test_it_stores_the_datalayer_it_was_given(self, topo):
        layer = analysis(projectName="IGNORED", dataLayer=topo)
        assert layer.datalayer is topo

    def test_the_project_name_argument_is_accepted_but_not_stored(self, topo):
        """__init__ takes projectName and does nothing with it."""
        layer = analysis(projectName="SOMETHING_ELSE", dataLayer=topo)
        assert not hasattr(layer, "_projectName")

    def test_the_datalayer_property_is_read_only(self, topo):
        layer = analysis(projectName="IGNORED", dataLayer=topo)
        with pytest.raises(AttributeError):
            layer.datalayer = topo


# --------------------------------------------------------------------------
# geoPandasToSTL
# --------------------------------------------------------------------------

@pytest.mark.unit
class TestGeoPandasToSTL:
    def test_it_returns_an_ascii_stl_string(self, topo, contours):
        stl = topo.geoPandasToSTL(contours, dxdy=100)
        assert stl.startswith("solid ")
        assert stl.rstrip().endswith("endsolid Topography")

    def test_every_triangle_is_a_well_formed_facet(self, topo, contours):
        stl = topo.geoPandasToSTL(contours, dxdy=100)
        assert stl.count("facet normal") == stl.count("endfacet")
        assert stl.count("vertex") == 3 * stl.count("endfacet")

    def test_the_vertices_span_the_gridded_part_of_the_contour_ramp(self, topo, contours):
        """The ramp is 0.1 m of height per metre of northing.

        ``numpy.mgrid[ymin:ymax:dxdy]`` excludes the endpoint, so at
        ``dxdy=100`` the northmost gridded row is y=660300, i.e. 300 m up the
        ramp from the 100 m contour -> 130 m.  ``rasterToSTL`` also emits a
        base skirt 10 m below the minimum, hence the -10.
        """
        stl = topo.geoPandasToSTL(contours, dxdy=100)
        zs = [
            float(line.split()[3])
            for line in stl.splitlines()
            if line.strip().startswith("vertex")
        ]
        assert min(zs) == pytest.approx(min(_HEIGHTS) - 10.0)
        assert max(zs) == pytest.approx(130.0)

    def test_a_finer_resolution_produces_more_facets(self, topo, contours):
        coarse = topo.geoPandasToSTL(contours, dxdy=200)
        fine = topo.geoPandasToSTL(contours, dxdy=100)
        assert fine.count("endfacet") > coarse.count("endfacet")

    def test_the_default_resolution_is_fifty_metres(self, topo, contours):
        assert topo.geoPandasToSTL(contours) == topo.geoPandasToSTL(contours, dxdy=50)

    @pytest.mark.xfail(
        strict=True,
        reason="B228: geoPandasToSTL accepts solidName but forwards the literal "
               "'Topography' to stlFactory.vectorToSTL, so the caller's solid "
               "name is discarded. See the consolidated findings issue.",
    )
    def test_the_requested_solid_name_reaches_the_stl(self, topo, contours):
        stl = topo.geoPandasToSTL(contours, dxdy=200, solidName="MyHill")
        assert stl.startswith("solid MyHill")

    def test_the_solid_name_argument_is_ignored(self, topo, contours):
        """Characterisation of B228."""
        stl = topo.geoPandasToSTL(contours, dxdy=200, solidName="MyHill")
        assert stl.startswith("solid Topography")
        assert "MyHill" not in stl


# --------------------------------------------------------------------------
# cutRegionFromSource / regionToSTL
# --------------------------------------------------------------------------

@pytest.mark.unit
class TestCutRegionFromSource:
    @pytest.mark.xfail(
        strict=True,
        reason="B227: TopographyToolkit.cutRegionFromSource calls "
               "super().cutRegionFromSource(shapeDataOrName=..., "
               "datasourceName=..., crs=...) but VectorToolkit's signature is "
               "(datasourceDocument, shape, isBounds, inputCRS); none of the "
               "three keyword names exist there, so the first statement always "
               "raises TypeError. See the consolidated findings issue.",
    )
    def test_it_cuts_the_requested_shape_out_of_the_datasource(self, topo, contours, tmp_path):
        path = tmp_path / "contours.geojson"
        contours.to_file(str(path), driver="GeoJSON")
        topo.addDataSource("CONTOURS", str(path), "geopandas", version=(0, 0, 1), crs=ITM)
        topo.setDataSourceDefaultVersion("CONTOURS", (0, 0, 1))

        result = topo.cutRegionFromSource(
            box(_X0, _YS[0], _X1, _YS[1]), datasourceName="CONTOURS", crs=ITM
        )
        assert isinstance(result, geopandas.GeoDataFrame)

    def test_any_call_dies_on_the_parent_keyword_names(self, topo):
        """Characterisation of B227."""
        with pytest.raises(TypeError) as excinfo:
            topo.cutRegionFromSource(
                box(0, 0, 1, 1), datasourceName="ANY", crs=ITM
            )
        assert "shapeDataOrName" in str(excinfo.value)

    def test_it_fails_before_looking_the_datasource_up(self, topo):
        """Characterisation of B227: a nonexistent datasource never matters."""
        with pytest.raises(TypeError):
            topo.cutRegionFromSource("NO_SUCH_SHAPE", datasourceName="NO_SUCH_SOURCE")


@pytest.mark.unit
class TestRegionToSTL:
    @pytest.mark.xfail(
        strict=True,
        reason="B227: regionToSTL's only substantive statement is a call to "
               "cutRegionFromSource, which always raises TypeError (see B227), "
               "so regionToSTL is unreachable too. See the consolidated "
               "findings issue.",
    )
    def test_it_returns_an_stl_string_for_a_populated_region(self, topo, contours, tmp_path):
        path = tmp_path / "contours.geojson"
        contours.to_file(str(path), driver="GeoJSON")
        topo.addDataSource("CONTOURS", str(path), "geopandas", version=(0, 0, 1), crs=ITM)
        topo.setDataSourceDefaultVersion("CONTOURS", (0, 0, 1))

        stl = topo.regionToSTL(
            box(_X0, _YS[0], _X1, _YS[-1]), dxdy=200, datasourceName="CONTOURS", crs=ITM
        )
        assert stl.startswith("solid ")

    def test_it_propagates_the_typeerror_from_cutregionfromsource(self, topo):
        """Characterisation of B227, seen through regionToSTL."""
        with pytest.raises(TypeError) as excinfo:
            topo.regionToSTL(box(0, 0, 1, 1), dxdy=50, datasourceName="ANY")
        assert "shapeDataOrName" in str(excinfo.value)


# --------------------------------------------------------------------------
# toDEM
# --------------------------------------------------------------------------

@pytest.mark.unit
class TestToDEMInputDispatch:
    def test_something_that_is_neither_a_string_nor_a_geodataframe_is_rejected(self, topo):
        with pytest.raises(ValueError, match="regionNameOrData must be"):
            topo.toDEM(12345)

    def test_a_bare_shapely_polygon_is_rejected(self, topo):
        with pytest.raises(ValueError, match="regionNameOrData must be"):
            topo.toDEM(Polygon([(0, 0), (1, 0), (1, 1)]))

    @pytest.mark.xfail(
        strict=True,
        reason="B229: toDEM resolves a datasource name through "
               "self.getDatasourceData (lower-case 's'); the toolkit defines "
               "only getDataSourceData, so the documented 'region name in the "
               "DB' input raises AttributeError -- which also makes the "
               "readGeoJSONString fallback on the next line dead code. "
               "See the consolidated findings issue.",
    )
    def test_a_registered_datasource_name_is_resolved(self, topo, contours, tmp_path):
        path = tmp_path / "contours.geojson"
        contours.to_file(str(path), driver="GeoJSON")
        topo.addDataSource("CONTOURS", str(path), "geopandas", version=(0, 0, 1), crs=ITM)
        topo.setDataSourceDefaultVersion("CONTOURS", (0, 0, 1))

        assert topo.toDEM("CONTOURS", dxdy=200).startswith(str(_X0))

    def test_a_string_input_dies_on_the_misspelled_accessor(self, topo):
        """Characterisation of B229."""
        with pytest.raises(AttributeError) as excinfo:
            topo.toDEM("CONTOURS")
        assert "getDatasourceData" in str(excinfo.value)

    def test_even_a_valid_geojson_string_never_reaches_the_fallback(self, topo, contours):
        """Characterisation of B229: the fallback branch is dead code."""
        with pytest.raises(AttributeError):
            topo.toDEM(contours.to_json())


@pytest.mark.unit
class TestToDEMSquareGrid:
    """A square grid is the only shape B230 does not turn into an IndexError."""

    @staticmethod
    def _square_contours():
        # 400 m x 400 m at dxdy=200 -> numpy.mgrid gives a 2x2 grid.
        return _contours(x0=_X0, x1=_X0 + 400.0,
                         ys=(660000.0, 660200.0, 660400.0))

    def test_the_header_carries_the_bounding_box_then_the_grid_shape(self, topo):
        dem = topo.toDEM(self._square_contours(), dxdy=200)
        bbox, shape = dem.splitlines()[0], dem.splitlines()[1]
        xmin, xmax, ymin, ymax = (float(v) for v in bbox.split())
        assert (xmin, xmax) == (_X0, _X0 + 200.0)
        assert (ymin, ymax) == (660000.0, 660200.0)
        assert shape.split() == ["2", "2"]

    def test_the_body_has_one_line_per_column_of_the_grid(self, topo):
        dem = topo.toDEM(self._square_contours(), dxdy=200)
        body = [line for line in dem.splitlines()[2:] if line.strip()]
        assert len(body) == 2
        assert all(len(line.split()) == 2 for line in body)

    def test_every_emitted_value_is_a_finite_height_in_the_contour_range(self, topo):
        dem = topo.toDEM(self._square_contours(), dxdy=200)
        values = [float(v) for line in dem.splitlines()[2:] for v in line.split()]
        assert values
        assert all(min(_HEIGHTS) <= v <= max(_HEIGHTS) for v in values)

    def test_the_heights_follow_the_linear_ramp_of_the_contours(self, topo):
        """0.1 m of height per metre of northing, from the 100/120/140 lines."""
        dem = topo.toDEM(self._square_contours(), dxdy=200)
        rows = [line.split() for line in dem.splitlines()[2:] if line.strip()]
        # y = 660000 -> 100 m ; y = 660200 -> 120 m, for both x columns.
        assert [float(v) for v in rows[0]] == pytest.approx([100.0, 100.0])
        assert [float(v) for v in rows[1]] == pytest.approx([120.0, 120.0])


@pytest.mark.unit
class TestToDEMNonSquareGrid:
    @staticmethod
    def _wide_contours():
        # 800 m of easting, 400 m of northing -> mgrid gives a 4x2 grid.
        return _contours(x0=_X0, x1=_X0 + 800.0,
                         ys=(660000.0, 660200.0, 660400.0))

    @pytest.mark.xfail(
        strict=True,
        reason="B230: toDEM sets Nx=grid_x.shape[0] and Ny=grid_x.shape[1], "
               "then indexes rasterized['height'][j, i] with i over Nx and j "
               "over Ny -- the subscripts are swapped, so any grid whose two "
               "dimensions differ raises IndexError (and a square one is "
               "silently transposed). See the consolidated findings issue.",
    )
    def test_a_grid_wider_than_it_is_tall_still_produces_a_dem(self, topo):
        dem = topo.toDEM(self._wide_contours(), dxdy=200)
        assert dem.splitlines()[1].split() == ["4", "2"]

    def test_a_non_square_grid_raises_indexerror(self, topo):
        """Characterisation of B230."""
        with pytest.raises(IndexError):
            topo.toDEM(self._wide_contours(), dxdy=200)

    def test_the_header_is_correct_even_though_the_body_cannot_be_written(self, topo):
        """Characterisation of B230: only the body loop is wrong."""
        rasterized = topo.stlFactory.rasterizeGeopandas(self._wide_contours(), dxdy=200)
        assert rasterized["x"].shape == (4, 2)
        with pytest.raises(IndexError):
            topo.toDEM(self._wide_contours(), dxdy=200)


@pytest.mark.unit
class TestToDEMMaskedInput:
    def test_a_masked_height_array_takes_the_mask_branch(self, topo, monkeypatch):
        """The `numpy.ma.is_masked` branch fills from the mask, not isnan."""
        grid_x, grid_y = numpy.mgrid[0.0:2.0:1.0, 0.0:2.0:1.0]
        height = numpy.ma.masked_array(
            numpy.array([[1.0, 2.0], [3.0, 0.0]]),
            mask=numpy.array([[False, False], [False, True]]),
        )

        def _fake_rasterize(gpandas, dxdy=50):
            return dict(x=grid_x, y=grid_y, height=height)

        monkeypatch.setattr(topo.stlFactory, "rasterizeGeopandas", _fake_rasterize)
        dem = topo.toDEM(_contours(), dxdy=1)
        values = [float(v) for line in dem.splitlines()[2:] for v in line.split()]
        # the masked cell is filled by nearest-neighbour, so nothing is NaN.
        assert not any(numpy.isnan(values))
        assert len(values) == 4


# --------------------------------------------------------------------------
# analysis.addHeight
# --------------------------------------------------------------------------

def _mesh():
    """A tiny cell-centre table: four columns of a mesh with a z coordinate."""
    return pandas.DataFrame(
        {
            "x": [0.0, 0.0, 10.0, 10.0],
            "y": [0.0, 10.0, 0.0, 10.0],
            "z": [5.0, 15.0, 5.0, 15.0],
        }
    )


def _ground(n=6):
    """A ground surface table: z rises 1 m per metre of x."""
    xs = numpy.linspace(0.0, 10.0, n)
    ys = numpy.linspace(0.0, 10.0, n)
    xx, yy = numpy.meshgrid(xs, ys)
    return pandas.DataFrame(
        {"x": xx.ravel(), "y": yy.ravel(), "z": xx.ravel()}
    )


@pytest.mark.unit
class TestAddHeight:
    @pytest.mark.xfail(
        strict=True,
        reason="B231: addHeight calls coordinateHandler.regularizeTimeSteps, "
               "but 'from hera.simulations.utils import coordinateHandler' "
               "binds the submodule of that name -- the class lives inside it "
               "and hera/simulations/utils/__init__.py is empty -- so the "
               "call raises AttributeError on the module. See the "
               "consolidated findings issue.",
    )
    def test_it_adds_ground_and_height_columns_to_the_mesh(self, topo):
        result = topo._analysis.addHeight(
            _mesh(), _ground(), resolution=2,
            saveMode="nosave",
        )
        assert "ground" in result.columns
        assert "height" in result.columns

    def test_it_raises_attributeerror_on_the_coordinate_handler_module(self, topo):
        """Characterisation of B231."""
        with pytest.raises(AttributeError) as excinfo:
            topo._analysis.addHeight(_mesh(), _ground(), resolution=2)
        message = str(excinfo.value)
        assert "regularizeTimeSteps" in message
        assert "coordinateHandler" in message

    def test_the_imported_name_is_a_module_not_the_class(self):
        """Characterisation of B231: the root cause, stated directly."""
        import types

        assert isinstance(vector_topography_module.coordinateHandler, types.ModuleType)
        assert not hasattr(vector_topography_module.coordinateHandler, "regularizeTimeSteps")
        assert hasattr(vector_topography_module.coordinateHandler, "coordinateHandler")


@pytest.mark.unit
class TestAddHeightSaveModeConstants:
    """B232, reached by neutralising B231 and nothing else."""

    @staticmethod
    def _working_handler(monkeypatch):
        """Bind the *class* the module meant to import, leaving B232 exposed."""
        from hera.simulations.utils.coordinateHandler import coordinateHandler

        monkeypatch.setattr(
            vector_topography_module, "coordinateHandler", coordinateHandler()
        )

    @pytest.mark.xfail(
        strict=True,
        reason="B232: addHeight compares saveMode against "
               "toolkit.TOOLKIT_SAVEMODE_FILEANDDB and friends, but the "
               "module-level name 'toolkit' was rebound by 'from "
               "hera.measurements.GIS.vector import toolkit' to the vector "
               "toolkit module, which defines no TOOLKIT_SAVEMODE_* "
               "constants -- so the save block raises AttributeError. "
               "See the consolidated findings issue.",
    )
    def test_the_default_save_mode_writes_the_cell_data_parquet(self, topo, monkeypatch):
        self._working_handler(monkeypatch)
        topo._analysis.addHeight(_mesh(), _ground(), resolution=2)
        assert os.path.exists(os.path.join(topo.filesDirectory, "cellData.parquet"))

    def test_the_save_block_raises_attributeerror_on_the_rebound_module(self, topo, monkeypatch):
        """Characterisation of B232."""
        self._working_handler(monkeypatch)
        with pytest.raises(AttributeError) as excinfo:
            topo._analysis.addHeight(_mesh(), _ground(), resolution=2)
        assert "TOOLKIT_SAVEMODE_FILEANDDB" in str(excinfo.value)

    def test_the_rebound_toolkit_module_has_none_of_the_save_constants(self):
        """Characterisation of B232: the root cause, stated directly."""
        module = vector_topography_module.toolkit
        assert module.__name__ == "hera.measurements.GIS.vector.toolkit"
        for name in (
            "TOOLKIT_SAVEMODE_ONLYFILE",
            "TOOLKIT_SAVEMODE_ONLYFILE_REPLACE",
            "TOOLKIT_SAVEMODE_FILEANDDB",
            "TOOLKIT_SAVEMODE_FILEANDDB_REPLACE",
        ):
            assert not hasattr(module, name)

    def test_the_columns_are_computed_before_the_save_block_fails(self, topo, monkeypatch):
        """Characterisation of B232: the arithmetic itself is sound.

        addHeight mutates the frame it is given, so the ground/height columns
        it computed survive the AttributeError and can be inspected.  Ground
        rises 1 m per metre of x, so at x=0 the ground is 0 and at x=10 it is
        10; height is clamped at 0 from below.
        """
        self._working_handler(monkeypatch)
        mesh = _mesh()
        with pytest.raises(AttributeError):
            topo._analysis.addHeight(mesh, _ground(), resolution=2)

        assert list(mesh["ground"]) == pytest.approx([0.0, 0.0, 10.0, 10.0], abs=1e-6)
        assert list(mesh["height"]) == pytest.approx([5.0, 15.0, 0.0, 5.0], abs=1e-6)

    @pytest.mark.parametrize(
        "saveMode",
        ["nosave", "onlyFile", "onlyFileReplace", "fileAndDB", "anything at all"],
    )
    def test_no_save_mode_escapes_the_broken_constant_lookup(self, topo, monkeypatch, saveMode):
        """Characterisation of B232: even NOSAVE fails.

        The fault is in the ``if saveMode in [...]`` *list literal*, which is
        built before the membership test runs, so no value of ``saveMode`` --
        not even one meant to skip saving -- can get past it.  addHeight can
        therefore never return at all.
        """
        self._working_handler(monkeypatch)
        with pytest.raises(AttributeError, match="TOOLKIT_SAVEMODE_FILEANDDB"):
            topo._analysis.addHeight(
                _mesh(), _ground(), resolution=2, saveMode=saveMode
            )

    def test_the_height_is_clamped_at_zero_below_the_ground(self, topo, monkeypatch):
        """Characterisation of B232: cells below ground get height 0, not a negative.

        Ground rises 1 m per metre of x, so at x=10 the ground is 10 m and a
        cell at z=1 sits 9 m below it.  The clamp runs before the save block,
        so the mutated frame still shows the result.
        """
        self._working_handler(monkeypatch)
        mesh = pandas.DataFrame({"x": [10.0], "y": [0.0], "z": [1.0]})
        with pytest.raises(AttributeError):
            topo._analysis.addHeight(mesh, _ground(), resolution=2)
        assert mesh["ground"].iloc[0] == pytest.approx(10.0, abs=1e-6)
        assert mesh["height"].iloc[0] == 0.0

    def test_the_computed_columns_are_appended_in_order(self, topo, monkeypatch):
        """Characterisation of B232: 'ground' then 'height', after the input columns."""
        self._working_handler(monkeypatch)
        mesh = _mesh()
        with pytest.raises(AttributeError):
            topo._analysis.addHeight(mesh, _ground(), resolution=2)
        assert list(mesh.columns) == ["x", "y", "z", "ground", "height"]
