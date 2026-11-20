"""hera/measurements/GIS/vector/demography.py -- the analysis and
presentation layers that ``test_gis_demography_toolkit.py`` leaves uncovered.

Covered here
------------
``analysis.createNewArea``, ``analysis.calculatePopulationInPolygon``, and the
whole ``presentation`` class: ``plotPopulationDensity``, ``plotPopulation``,
``plotPopulationByType``, ``plotPopulationInPolygon``, ``plotArea`` and
``plotPopulationOnMap``.

How the inputs were produced
----------------------------
``_census()`` builds five 100 m x 100 m squares **in ITM** (EPSG:2039) with
known population counts.  Working natively in ITM is what makes the
area-fraction arithmetic exact rather than approximate: both
``calculatePopulationInPolygon`` and the presentation layer reproject to ITM
before measuring, so ``to_crs(ITM)`` is the identity here and no projection
distortion enters the expected numbers.  Four of the squares form a 2x2 block
around the origin; the fifth sits 1 km away so that the "does not intersect
the query polygon" branch is exercised.  The query polygon
``box(50, 50, 150, 150)`` overlaps exactly one quarter of each of the four,
so every expected value below is a hand-computable quarter of a known count.

``plotPopulationOnMap`` needs a live ``TilesToolkit`` (it fetches map tiles
over the network), so it is driven through ``_FakeTiles``, a recorder that
captures the corner/zoom/CRS arguments and draws nothing.  That keeps the
assertion on what the method is responsible for -- the extent it asks for, in
which CRS, and the choropleth it overlays -- rather than on imagery.

Bugs pinned (each has an xfail(strict) for the intended behaviour plus a
passing characterisation of what actually happens)
--------------------------------------------------------------------------
* **B233** -- ``createNewArea``'s ``TOOLKIT_SAVEMODE_FILEANDDB`` branch reads
  ``toolkit.TOOLKIT_DATASOURCE_NAME`` / ``toolkit.TOOLKIT_TOOLKITNAME_FIELD``
  / ``toolkit.TOOLKIT_DATASOURCE_TYPE``.  The module imports ``from hera
  import toolkit`` on line 9 and then *rebinds the same name* on line 10 with
  ``from hera.measurements.GIS.vector import toolkit``; the vector toolkit
  module defines none of those constants, so asking to save to the DB raises
  ``AttributeError`` after the file has already been written.
* **B234** -- the same branch is guarded by ``== TOOLKIT_SAVEMODE_FILEANDDB``
  only, so ``TOOLKIT_SAVEMODE_FILEANDDB_REPLACE`` ("save to a file and store
  to the DB as a source, replacing the entry if it exists", per the
  docstring) writes the file and silently never touches the database.
* **B235** -- ``TOOLKIT_SAVEMODE_ONLYFILE`` is documented to "raise exception
  if file exists" and ``..._REPLACE`` to replace it, but the two are handled
  identically: nothing checks for an existing file, so ONLYFILE overwrites
  instead of raising.

Deliberately not tested, and why
--------------------------------
* ``createNewArea``'s ``shapeNameOrData`` string branch: it goes through
  ``self.datalayer.shapes``, which is B81 (already pinned in
  ``test_gis_demography_toolkit.py``) -- the ``_shapes`` attribute is never
  assigned anywhere in the class, so the branch raises ``AttributeError``
  before any of this file's logic runs.  Not re-pinned here.
* The ``TopologicalError`` recovery path in ``createNewArea``: it needs a
  geometry whose ``intersection`` raises inside GEOS, which shapely 2 no
  longer does for the usual self-intersecting inputs (it returns a result
  instead).  Left for an integration test rather than faked.
* Pixel-level appearance of the plots: the tests assert on the artists,
  labels, limits and colour-bar text, not on rendered images.
"""
import os

import geopandas
import matplotlib.pyplot as plt
import pytest
from shapely.geometry import box

from hera import toolkitHome
from hera.measurements.GIS.utils import ITM, WSG84
from hera.toolkit import (
    TOOLKIT_SAVEMODE_FILEANDDB,
    TOOLKIT_SAVEMODE_FILEANDDB_REPLACE,
    TOOLKIT_SAVEMODE_NOSAVE,
    TOOLKIT_SAVEMODE_ONLYFILE,
    TOOLKIT_SAVEMODE_ONLYFILE_REPLACE,
)

# --------------------------------------------------------------------------
# synthetic census input
# --------------------------------------------------------------------------

# Four adjoining 100 m squares plus one far away.  total_pop doubles clockwise
# so that no two subsets can share a sum by accident.
_CELLS = [
    dict(name="A", geometry=box(0, 0, 100, 100), total_pop=1000, age_0_14=200),
    dict(name="B", geometry=box(100, 0, 200, 100), total_pop=2000, age_0_14=400),
    dict(name="C", geometry=box(0, 100, 100, 200), total_pop=4000, age_0_14=800),
    dict(name="D", geometry=box(100, 100, 200, 200), total_pop=8000, age_0_14=1600),
    dict(name="FAR", geometry=box(1000, 1000, 1100, 1100), total_pop=99, age_0_14=9),
]

_QUERY = box(50, 50, 150, 150)   # one quarter of each of A, B, C and D
_NEIGHBOURS = ("A", "B", "C", "D")
_TOTAL_IN_QUERY = 0.25 * (1000 + 2000 + 4000 + 8000)      # 3750
_CHILDREN_IN_QUERY = 0.25 * (200 + 400 + 800 + 1600)      # 750


def _census(crs=ITM):
    return geopandas.GeoDataFrame(_CELLS, crs=crs)


class _FakeTiles:
    """Stand-in for TilesToolkit: records the request, draws nothing."""

    def __init__(self):
        self.corner_calls = []
        self.plot_calls = []
        self.presentation = self

    def getImageFromCorners(self, **kwargs):
        self.corner_calls.append(kwargs)
        return "TILE-IMAGE"

    def plot(self, image, **kwargs):
        self.plot_calls.append((image, kwargs))
        return kwargs.get("ax")


@pytest.fixture()
def demography(unit_toolkit_factory):
    return unit_toolkit_factory(toolkitHome.GIS_DEMOGRAPHY)


@pytest.fixture()
def census():
    return _census()


@pytest.fixture(autouse=True)
def _fresh_figures():
    plt.close("all")
    yield
    plt.close("all")


# ==========================================================================
# analysis.calculatePopulationInPolygon
# ==========================================================================

@pytest.mark.unit
class TestCalculatePopulationInPolygon:
    def test_only_the_intersecting_census_cells_survive(self, demography, census):
        result = demography.analysis.calculatePopulationInPolygon(_QUERY, census)
        assert len(result) == len(_NEIGHBOURS)

    def test_each_intersection_covers_a_quarter_of_its_cell(self, demography, census):
        """box(50,50,150,150) takes a 50 m x 50 m bite out of each 100 m cell."""
        result = demography.analysis.calculatePopulationInPolygon(_QUERY, census)
        assert list(result["areaFraction"]) == pytest.approx([0.25] * 4)

    def test_the_population_is_scaled_by_the_area_fraction(self, demography, census):
        result = demography.analysis.calculatePopulationInPolygon(_QUERY, census)
        assert list(result["total_pop"]) == pytest.approx([250.0, 500.0, 1000.0, 2000.0])

    def test_the_total_is_a_quarter_of_the_four_neighbouring_cells(self, demography, census):
        result = demography.analysis.calculatePopulationInPolygon(_QUERY, census)
        assert result["total_pop"].sum() == pytest.approx(_TOTAL_IN_QUERY)

    def test_every_population_column_present_in_the_data_is_scaled(self, demography, census):
        result = demography.analysis.calculatePopulationInPolygon(_QUERY, census)
        assert result["age_0_14"].sum() == pytest.approx(_CHILDREN_IN_QUERY)

    def test_population_columns_absent_from_the_data_are_not_invented(self, demography, census):
        """The toolkit knows six groups; the data carries two."""
        result = demography.analysis.calculatePopulationInPolygon(_QUERY, census)
        assert "age_65_up" not in result.columns

    def test_the_geometry_returned_is_the_clipped_intersection(self, demography, census):
        result = demography.analysis.calculatePopulationInPolygon(_QUERY, census)
        assert list(result.geometry.area) == pytest.approx([2500.0] * 4)

    def test_a_single_population_type_can_be_named_by_its_group_label(self, demography, census):
        result = demography.analysis.calculatePopulationInPolygon(
            _QUERY, census, populationTypes="Children"
        )
        assert result["age_0_14"].sum() == pytest.approx(_CHILDREN_IN_QUERY)
        assert "total_pop" not in result.columns

    def test_a_population_type_can_also_be_named_by_its_column(self, demography, census):
        """``populationTypes.get(x, x)`` lets the raw column name through."""
        result = demography.analysis.calculatePopulationInPolygon(
            _QUERY, census, populationTypes=["total_pop"]
        )
        assert "total_pop" in result.columns
        assert "age_0_14" not in result.columns

    def test_a_polygon_that_misses_every_cell_gives_an_empty_result(self, demography, census):
        result = demography.analysis.calculatePopulationInPolygon(
            box(5000, 5000, 5100, 5100), census
        )
        assert len(result) == 0

    def test_a_polygon_containing_a_whole_cell_gives_it_a_fraction_of_one(self, demography, census):
        result = demography.analysis.calculatePopulationInPolygon(
            box(-10, -10, 110, 110), census
        )
        assert result.loc[0, "areaFraction"] == pytest.approx(1.0)
        assert result.loc[0, "total_pop"] == pytest.approx(1000.0)


@pytest.mark.unit
class TestCalculatePopulationInPolygonDataSources:
    def test_a_registered_datasource_name_is_resolved(self, demography, census, tmp_path):
        path = tmp_path / "census.geojson"
        census.to_file(str(path), driver="GeoJSON")
        demography.addDataSource("CENSUS", str(path), "geopandas", version=(0, 0, 1), crs=ITM)
        demography.setDataSourceDefaultVersion("CENSUS", (0, 0, 1))

        result = demography.analysis.calculatePopulationInPolygon(_QUERY, "CENSUS")
        assert result["total_pop"].sum() == pytest.approx(_TOTAL_IN_QUERY)

    def test_a_bare_filesystem_path_is_read_directly(self, demography, census, tmp_path):
        """The `os.path.exists(dataSourceOrData)` fallback branch."""
        path = tmp_path / "census.geojson"
        census.to_file(str(path), driver="GeoJSON")

        result = demography.analysis.calculatePopulationInPolygon(_QUERY, str(path))
        assert result["total_pop"].sum() == pytest.approx(_TOTAL_IN_QUERY)

    def test_a_geojson_string_is_parsed(self, demography, census):
        """The `readGeoJSONString` fallback branch."""
        result = demography.analysis.calculatePopulationInPolygon(_QUERY, census.to_json())
        assert result["total_pop"].sum() == pytest.approx(_TOTAL_IN_QUERY)


# ==========================================================================
# analysis.createNewArea
# ==========================================================================

@pytest.mark.unit
class TestCreateNewAreaShapelyPolygon:
    def test_it_returns_a_single_row_carrying_the_query_polygon(self, demography, census):
        area = demography.analysis.createNewArea(_QUERY, census).getData()
        assert len(area) == 1
        assert area.geometry.iloc[0].equals(_QUERY)

    def test_the_population_is_the_unscaled_sum_of_the_intersecting_cells(self, demography, census):
        """createNewArea sums whole cells -- no area weighting, unlike
        calculatePopulationInPolygon."""
        area = demography.analysis.createNewArea(_QUERY, census).getData()
        assert area["total_pop"].iloc[0] == 1000 + 2000 + 4000 + 8000

    def test_the_distant_cell_is_excluded(self, demography, census):
        area = demography.analysis.createNewArea(_QUERY, census).getData()
        assert area["total_pop"].iloc[0] != sum(cell["total_pop"] for cell in _CELLS)

    def test_every_population_column_in_the_data_is_summed(self, demography, census):
        area = demography.analysis.createNewArea(_QUERY, census).getData()
        assert area["age_0_14"].iloc[0] == 200 + 400 + 800 + 1600

    def test_columns_the_data_does_not_carry_are_skipped(self, demography, census):
        area = demography.analysis.createNewArea(_QUERY, census).getData()
        assert "age_65_up" not in area.columns

    def test_the_result_inherits_the_crs_of_the_population_data(self, demography, census):
        area = demography.analysis.createNewArea(_QUERY, census).getData()
        assert area.crs.to_epsg() == ITM

    def test_an_explicit_population_type_list_narrows_the_output(self, demography, census):
        area = demography.analysis.createNewArea(
            _QUERY, census, populationTypes=["age_0_14"]
        ).getData()
        assert "age_0_14" in area.columns
        assert "total_pop" not in area.columns

    def test_the_default_save_mode_wraps_the_frame_rather_than_storing_it(self, demography, census):
        from hera.datalayer import nonDBMetadataFrame

        result = demography.analysis.createNewArea(_QUERY, census)
        assert isinstance(result, nonDBMetadataFrame)


@pytest.mark.unit
class TestCreateNewAreaGeoDataFrameShape:
    def test_a_convex_hull_of_the_shape_is_used_when_convex_is_true(self, demography, census):
        """ConvexPolygons buffers by 100 m, so the area grows well past the box."""
        shape = geopandas.GeoDataFrame({"geometry": [_QUERY]}, crs=ITM)
        area = demography.analysis.createNewArea(shape, census, convex=True).getData()
        assert area.geometry.iloc[0].area > _QUERY.area

    def test_the_shape_is_used_verbatim_when_convex_is_false(self, demography, census):
        shape = geopandas.GeoDataFrame({"geometry": [_QUERY]}, crs=ITM)
        area = demography.analysis.createNewArea(shape, census, convex=False).getData()
        assert area.geometry.iloc[0].area == pytest.approx(_QUERY.area)

    def test_a_multi_row_shape_is_unioned_when_convex_is_false(self, demography, census):
        shape = geopandas.GeoDataFrame(
            {"geometry": [box(0, 0, 10, 10), box(20, 20, 30, 30)]}, crs=ITM
        )
        area = demography.analysis.createNewArea(shape, census, convex=False).getData()
        assert area.geometry.iloc[0].area == pytest.approx(200.0)


@pytest.mark.unit
class TestCreateNewAreaSaveModeValidation:
    @pytest.mark.parametrize(
        "saveMode",
        [
            TOOLKIT_SAVEMODE_ONLYFILE,
            TOOLKIT_SAVEMODE_ONLYFILE_REPLACE,
            TOOLKIT_SAVEMODE_FILEANDDB,
            TOOLKIT_SAVEMODE_FILEANDDB_REPLACE,
        ],
    )
    def test_a_region_name_is_required_for_every_saving_mode(self, demography, census, saveMode):
        with pytest.raises(ValueError, match="Must specify regionName"):
            demography.analysis.createNewArea(_QUERY, census, saveMode=saveMode)

    def test_no_region_name_is_needed_when_nothing_is_saved(self, demography, census):
        result = demography.analysis.createNewArea(
            _QUERY, census, saveMode=TOOLKIT_SAVEMODE_NOSAVE, regionName=None
        )
        assert result.getData() is not None


@pytest.mark.unit
class TestCreateNewAreaFileOutput:
    def test_a_region_name_with_an_extension_is_used_verbatim(self, demography, census):
        demography.analysis.createNewArea(
            _QUERY, census, saveMode=TOOLKIT_SAVEMODE_ONLYFILE, regionName="myarea.geojson"
        )
        assert os.path.exists(os.path.join(demography.filesDirectory, "myarea.geojson"))

    def test_a_region_name_without_one_gets_a_shp_suffix(self, demography, census):
        demography.analysis.createNewArea(
            _QUERY, census, saveMode=TOOLKIT_SAVEMODE_ONLYFILE, regionName="MYAREA"
        )
        assert os.path.exists(os.path.join(demography.filesDirectory, "MYAREA.shp"))

    def test_the_written_file_round_trips_the_population(self, demography, census):
        demography.analysis.createNewArea(
            _QUERY, census, saveMode=TOOLKIT_SAVEMODE_ONLYFILE, regionName="myarea.geojson"
        )
        written = geopandas.read_file(
            os.path.join(demography.filesDirectory, "myarea.geojson")
        )
        assert written["total_pop"].iloc[0] == 15000

    @pytest.mark.xfail(
        strict=True,
        reason="B235: createNewArea documents TOOLKIT_SAVEMODE_ONLYFILE as "
               "'raise exception if file exists' and ..._REPLACE as 'replace "
               "the file if it exists', but the two modes share one code path "
               "that never checks os.path.exists -- so ONLYFILE silently "
               "overwrites. See the consolidated findings issue.",
    )
    def test_only_file_refuses_to_clobber_an_existing_file(self, demography, census):
        demography.analysis.createNewArea(
            _QUERY, census, saveMode=TOOLKIT_SAVEMODE_ONLYFILE, regionName="myarea.geojson"
        )
        with pytest.raises((FileExistsError, ValueError)):
            demography.analysis.createNewArea(
                _QUERY, census, saveMode=TOOLKIT_SAVEMODE_ONLYFILE, regionName="myarea.geojson"
            )

    def test_only_file_overwrites_silently_on_the_second_call(self, demography, census):
        """Characterisation of B235."""
        smaller = geopandas.GeoDataFrame([_CELLS[0]], crs=ITM)
        demography.analysis.createNewArea(
            _QUERY, census, saveMode=TOOLKIT_SAVEMODE_ONLYFILE, regionName="myarea.geojson"
        )
        demography.analysis.createNewArea(
            _QUERY, smaller, saveMode=TOOLKIT_SAVEMODE_ONLYFILE, regionName="myarea.geojson"
        )
        written = geopandas.read_file(
            os.path.join(demography.filesDirectory, "myarea.geojson")
        )
        assert written["total_pop"].iloc[0] == 1000

    def test_only_file_and_only_file_replace_behave_identically(self, demography, census):
        """Characterisation of B235: the two modes are indistinguishable."""
        for mode, name in (
            (TOOLKIT_SAVEMODE_ONLYFILE, "plain.geojson"),
            (TOOLKIT_SAVEMODE_ONLYFILE_REPLACE, "replace.geojson"),
        ):
            demography.analysis.createNewArea(
                _QUERY, census, saveMode=mode, regionName=name
            )
            demography.analysis.createNewArea(
                _QUERY, census, saveMode=mode, regionName=name
            )
            assert os.path.exists(os.path.join(demography.filesDirectory, name))


@pytest.mark.unit
class TestCreateNewAreaDatabaseOutput:
    @pytest.mark.xfail(
        strict=True,
        reason="B233: the TOOLKIT_SAVEMODE_FILEANDDB branch reads "
               "toolkit.TOOLKIT_DATASOURCE_NAME, but the module-level name "
               "'toolkit' is rebound on line 10 from hera.toolkit to "
               "hera.measurements.GIS.vector.toolkit, which defines no such "
               "constant -- so saving to the DB raises AttributeError after "
               "the file has already been written. See the consolidated "
               "findings issue.",
    )
    def test_saving_to_the_database_returns_the_stored_document(self, demography, census):
        doc = demography.analysis.createNewArea(
            _QUERY, census, saveMode=TOOLKIT_SAVEMODE_FILEANDDB, regionName="dbarea.geojson"
        )
        assert doc.resource.endswith("dbarea.geojson")

    def test_saving_to_the_database_raises_on_the_rebound_module(self, demography, census):
        """Characterisation of B233."""
        with pytest.raises(AttributeError) as excinfo:
            demography.analysis.createNewArea(
                _QUERY, census, saveMode=TOOLKIT_SAVEMODE_FILEANDDB, regionName="dbarea.geojson"
            )
        assert "TOOLKIT_DATASOURCE_NAME" in str(excinfo.value)

    def test_the_file_is_written_before_the_database_step_fails(self, demography, census):
        """Characterisation of B233: the failure leaves an orphan file behind."""
        with pytest.raises(AttributeError):
            demography.analysis.createNewArea(
                _QUERY, census, saveMode=TOOLKIT_SAVEMODE_FILEANDDB, regionName="dbarea.geojson"
            )
        assert os.path.exists(os.path.join(demography.filesDirectory, "dbarea.geojson"))

    def test_the_rebound_module_carries_none_of_the_datasource_constants(self):
        """Characterisation of B233: the root cause, stated directly."""
        from hera.measurements.GIS.vector import demography as demography_module

        assert demography_module.toolkit.__name__ == "hera.measurements.GIS.vector.toolkit"
        for name in (
            "TOOLKIT_DATASOURCE_NAME",
            "TOOLKIT_TOOLKITNAME_FIELD",
            "TOOLKIT_DATASOURCE_TYPE",
        ):
            assert not hasattr(demography_module.toolkit, name)

    @pytest.mark.xfail(
        strict=True,
        reason="B234: the database branch is guarded by '== "
               "TOOLKIT_SAVEMODE_FILEANDDB' alone, so "
               "TOOLKIT_SAVEMODE_FILEANDDB_REPLACE -- documented as 'save to "
               "a file and store to the DB as a source, replacing the entry "
               "if it exists' -- writes the file and never reaches the "
               "database at all. See the consolidated findings issue.",
    )
    def test_the_replacing_database_mode_also_stores_a_document(self, demography, census):
        demography.analysis.createNewArea(
            _QUERY, census,
            saveMode=TOOLKIT_SAVEMODE_FILEANDDB_REPLACE, regionName="rep.geojson",
        )
        assert len(demography.getCacheDocuments(type="ToolkitDataSource")) == 1

    def test_the_replacing_database_mode_writes_only_a_file(self, demography, census):
        """Characterisation of B234."""
        result = demography.analysis.createNewArea(
            _QUERY, census,
            saveMode=TOOLKIT_SAVEMODE_FILEANDDB_REPLACE, regionName="rep.geojson",
        )
        from hera.datalayer import nonDBMetadataFrame

        assert isinstance(result, nonDBMetadataFrame)
        assert os.path.exists(os.path.join(demography.filesDirectory, "rep.geojson"))
        assert len(demography.getCacheDocuments(type="ToolkitDataSource")) == 0


# ==========================================================================
# presentation
# ==========================================================================

@pytest.mark.unit
class TestPresentationWiring:
    def test_the_datalayer_property_points_back_at_the_toolkit(self, demography):
        assert demography.presentation.datalayer is demography

    def test_the_analysis_layer_points_back_at_the_toolkit_too(self, demography):
        assert demography.analysis.datalayer is demography


@pytest.mark.unit
class TestPlotPopulationDensity:
    def test_it_returns_axes_carrying_one_collection_per_call(self, demography, census):
        ax = demography.presentation.plotPopulationDensity(census)
        assert len(ax.collections) == 1

    def test_it_creates_its_own_figure_when_none_is_given(self, demography, census):
        ax = demography.presentation.plotPopulationDensity(census)
        assert ax.figure is not None

    def test_it_draws_on_the_axes_it_is_handed(self, demography, census):
        fig, ax = plt.subplots()
        assert demography.presentation.plotPopulationDensity(census, ax=ax) is ax

    def test_the_input_data_is_not_mutated(self, demography, census):
        before = list(census.columns)
        demography.presentation.plotPopulationDensity(census)
        assert list(census.columns) == before

    def test_the_colour_bar_names_the_density_unit(self, demography, census):
        ax = demography.presentation.plotPopulationDensity(census)
        label = ax.figure.axes[-1].get_ylabel()
        assert label.startswith("Population density [people/")

    def test_an_explicit_colour_bar_label_wins(self, demography, census):
        ax = demography.presentation.plotPopulationDensity(
            census, colorbar_label="people per square kilometre"
        )
        assert ax.figure.axes[-1].get_ylabel() == "people per square kilometre"

    def test_no_colour_bar_axes_are_added_when_it_is_switched_off(self, demography, census):
        ax = demography.presentation.plotPopulationDensity(census, colorbar=False)
        assert len(ax.figure.axes) == 1

    def test_the_colour_bar_range_matches_the_computed_density(self, demography, census):
        """A 100 m square is 0.01 km2, so 1000 people there is 1e5 per km2."""
        ax = demography.presentation.plotPopulationDensity(census)
        cbar_axes = ax.figure.axes[-1]
        low, high = cbar_axes.get_ylim()
        assert low == pytest.approx(9900.0)     # the FAR cell: 99 / 0.01
        assert high == pytest.approx(800000.0)  # cell D: 8000 / 0.01

    def test_an_explicit_colour_scale_overrides_the_data_range(self, demography, census):
        ax = demography.presentation.plotPopulationDensity(census, vmin=0, vmax=10)
        assert ax.figure.axes[-1].get_ylim() == pytest.approx((0.0, 10.0))

    def test_a_smaller_density_unit_lowers_the_numbers(self, demography, census):
        from hera.utils.unitHandler import ureg

        per_km2 = demography.presentation.plotPopulationDensity(census)
        per_m2 = demography.presentation.plotPopulationDensity(
            census, density_units=ureg.m ** 2
        )
        assert per_m2.figure.axes[-1].get_ylim()[1] < per_km2.figure.axes[-1].get_ylim()[1]

    def test_the_title_and_limits_are_forwarded(self, demography, census):
        ax = demography.presentation.plotPopulationDensity(
            census, title="Density", xlim=(0, 50), ylim=(0, 60)
        )
        assert ax.get_title() == "Density"
        assert ax.get_xlim() == pytest.approx((0.0, 50.0))
        assert ax.get_ylim() == pytest.approx((0.0, 60.0))

    def test_an_input_crs_can_be_stamped_on_unprojected_data(self, demography):
        """set_crs(allow_override=True) then to_crs: the data arrives in ITM."""
        raw = geopandas.GeoDataFrame(_CELLS)  # no CRS at all
        ax = demography.presentation.plotPopulationDensity(
            raw, inputCRS=ITM, outputCRS=ITM
        )
        assert len(ax.collections) == 1

    def test_reprojection_to_wgs84_shrinks_the_drawn_extent(self, demography, census):
        itm = demography.presentation.plotPopulationDensity(census, outputCRS=ITM)
        wgs = demography.presentation.plotPopulationDensity(census, outputCRS=WSG84)
        itm_width = itm.get_xlim()[1] - itm.get_xlim()[0]
        wgs_width = wgs.get_xlim()[1] - wgs.get_xlim()[0]
        assert wgs_width < itm_width


@pytest.mark.unit
class TestPlotPopulation:
    def test_it_returns_axes_with_the_choropleth_drawn(self, demography, census):
        ax = demography.presentation.plotPopulation(census)
        assert len(ax.collections) == 1

    def test_the_colour_bar_is_labelled_population_count(self, demography, census):
        ax = demography.presentation.plotPopulation(census)
        assert ax.figure.axes[-1].get_ylabel() == "Population count"

    def test_the_colour_bar_range_is_the_raw_population_range(self, demography, census):
        ax = demography.presentation.plotPopulation(census)
        assert ax.figure.axes[-1].get_ylim() == pytest.approx((99.0, 8000.0))

    def test_a_different_population_column_can_be_plotted(self, demography, census):
        ax = demography.presentation.plotPopulation(census, populationType="age_0_14")
        assert ax.figure.axes[-1].get_ylim() == pytest.approx((9.0, 1600.0))

    def test_a_missing_column_is_reported_rather_than_silently_ignored(self, demography, census):
        with pytest.raises(Exception):
            demography.presentation.plotPopulation(census, populationType="no_such_column")

    def test_the_title_is_forwarded(self, demography, census):
        ax = demography.presentation.plotPopulation(census, title="Counts")
        assert ax.get_title() == "Counts"

    def test_it_draws_on_a_supplied_axes(self, demography, census):
        fig, ax = plt.subplots()
        assert demography.presentation.plotPopulation(census, ax=ax) is ax


@pytest.mark.unit
class TestPlotPopulationByType:
    def test_it_returns_a_figure_not_axes(self, demography, census):
        fig = demography.presentation.plotPopulationByType(census)
        assert isinstance(fig, plt.Figure)

    def test_the_grid_has_the_requested_number_of_columns(self, demography, census):
        fig = demography.presentation.plotPopulationByType(census, ncols=3)
        # 2 population columns in the data, 3 columns of subplots -> one row.
        assert len(fig.axes) >= 3

    def test_only_the_population_columns_present_in_the_data_get_a_panel(self, demography, census):
        fig = demography.presentation.plotPopulationByType(census)
        titles = [ax.get_title() for ax in fig.axes if ax.get_title()]
        assert set(titles) == {"All", "Children"}

    def test_unused_panels_are_hidden(self, demography, census):
        fig = demography.presentation.plotPopulationByType(census, ncols=3)
        visible = [ax for ax in fig.axes if ax.get_visible() and ax.get_title()]
        assert len(visible) == 2

    def test_an_explicit_column_list_is_honoured(self, demography, census):
        fig = demography.presentation.plotPopulationByType(
            census, populationTypes=["total_pop"], ncols=1
        )
        titles = [ax.get_title() for ax in fig.axes if ax.get_title()]
        assert titles == ["All"]

    def test_a_column_with_no_group_label_keeps_its_own_name_as_the_title(self, demography, census):
        data = census.assign(other_pop=1)
        fig = demography.presentation.plotPopulationByType(
            data, populationTypes=["other_pop"], ncols=1
        )
        titles = [ax.get_title() for ax in fig.axes if ax.get_title()]
        assert titles == ["other_pop"]

    def test_a_single_column_grid_still_works(self, demography, census):
        """np.atleast_1d(axes).flatten() has to cope with a non-array axes."""
        fig = demography.presentation.plotPopulationByType(
            census, populationTypes=["total_pop"], ncols=1
        )
        assert len(fig.axes) >= 1


@pytest.mark.unit
class TestPlotPopulationInPolygon:
    @pytest.fixture()
    def intersection(self, demography, census):
        return demography.analysis.calculatePopulationInPolygon(_QUERY, census)

    def test_it_plots_the_intersection_result(self, demography, intersection):
        ax = demography.presentation.plotPopulationInPolygon(intersection)
        assert len(ax.collections) >= 1

    def test_the_colour_bar_names_the_population_column(self, demography, intersection):
        ax = demography.presentation.plotPopulationInPolygon(intersection)
        assert ax.figure.axes[-1].get_ylabel() == "Population (total_pop)"

    def test_the_colour_bar_range_comes_from_the_scaled_population(self, demography, intersection):
        ax = demography.presentation.plotPopulationInPolygon(intersection)
        assert ax.figure.axes[-1].get_ylim() == pytest.approx((250.0, 2000.0))

    def test_context_data_adds_a_second_layer_underneath(self, demography, intersection, census):
        without = demography.presentation.plotPopulationInPolygon(intersection)
        with_ctx = demography.presentation.plotPopulationInPolygon(
            intersection, contextData=census
        )
        assert len(with_ctx.collections) == len(without.collections) + 1

    def test_the_query_polygon_is_outlined_when_supplied(self, demography, intersection):
        without = demography.presentation.plotPopulationInPolygon(intersection)
        with_poly = demography.presentation.plotPopulationInPolygon(
            intersection, queryPolygon=_QUERY, outputCRS=ITM
        )
        assert len(with_poly.collections) == len(without.collections) + 1

    def test_a_result_without_the_population_column_is_drawn_in_a_flat_colour(self, demography, intersection):
        """The `else` branch: no column to colour by, so no colour bar."""
        geometry_only = intersection[["geometry"]].copy()
        geometry_only.crs = intersection.crs
        ax = demography.presentation.plotPopulationInPolygon(geometry_only)
        assert len(ax.figure.axes) == 1

    def test_the_title_and_limits_are_forwarded(self, demography, intersection):
        ax = demography.presentation.plotPopulationInPolygon(
            intersection, title="Affected", xlim=(0, 10), ylim=(0, 20)
        )
        assert ax.get_title() == "Affected"
        assert ax.get_xlim() == pytest.approx((0.0, 10.0))
        assert ax.get_ylim() == pytest.approx((0.0, 20.0))


@pytest.mark.unit
class TestPlotArea:
    @pytest.fixture()
    def area(self, demography, census):
        return demography.analysis.createNewArea(_QUERY, census).getData()

    def test_it_draws_the_area_polygon(self, demography, area):
        ax = demography.presentation.plotArea(area)
        assert len(ax.collections) == 1

    def test_the_population_is_annotated_at_the_centroid(self, demography, area):
        ax = demography.presentation.plotArea(area, annotate=True)
        texts = [t.get_text() for t in ax.texts]
        assert texts == ["15,000"]

    def test_the_annotation_sits_on_the_polygon_centroid(self, demography, area):
        ax = demography.presentation.plotArea(area, annotate=True)
        centroid = area.geometry.iloc[0].centroid
        assert ax.texts[0].xy == pytest.approx((centroid.x, centroid.y))

    def test_annotation_can_be_switched_off(self, demography, area):
        ax = demography.presentation.plotArea(area, annotate=False)
        assert len(ax.texts) == 0

    def test_nothing_is_annotated_when_the_column_is_absent(self, demography, area):
        ax = demography.presentation.plotArea(area, populationType="age_65_up")
        assert len(ax.texts) == 0

    def test_context_data_adds_a_layer(self, demography, area, census):
        ax = demography.presentation.plotArea(area, contextData=census)
        assert len(ax.collections) == 2

    def test_the_title_is_forwarded(self, demography, area):
        assert demography.presentation.plotArea(area, title="Area").get_title() == "Area"


@pytest.mark.unit
class TestPlotPopulationOnMap:
    def test_the_tile_extent_is_requested_in_wgs84(self, demography, census):
        tiles = _FakeTiles()
        demography.presentation.plotPopulationOnMap(census, tiles)
        call = tiles.corner_calls[0]
        assert call["inputCRS"] == WSG84
        expected = census.to_crs(epsg=WSG84).total_bounds
        assert (call["minx"], call["miny"], call["maxx"], call["maxy"]) == pytest.approx(
            tuple(expected)
        )

    def test_the_zoom_level_is_forwarded(self, demography, census):
        tiles = _FakeTiles()
        demography.presentation.plotPopulationOnMap(census, tiles, zoomlevel=11)
        assert tiles.corner_calls[0]["zoomlevel"] == 11

    def test_the_tile_server_is_forwarded(self, demography, census):
        tiles = _FakeTiles()
        demography.presentation.plotPopulationOnMap(census, tiles, tileServer="OSM")
        assert tiles.corner_calls[0]["tileServer"] == "OSM"

    def test_the_output_crs_is_forwarded_to_the_tile_request(self, demography, census):
        tiles = _FakeTiles()
        demography.presentation.plotPopulationOnMap(census, tiles, outputCRS=ITM)
        assert tiles.corner_calls[0]["outputCRS"] == ITM

    def test_itm_is_used_when_no_output_crs_is_given(self, demography, census):
        tiles = _FakeTiles()
        demography.presentation.plotPopulationOnMap(census, tiles, outputCRS=None)
        assert tiles.corner_calls[0]["outputCRS"] == ITM

    def test_the_fetched_image_is_drawn_before_the_choropleth(self, demography, census):
        tiles = _FakeTiles()
        ax = demography.presentation.plotPopulationOnMap(census, tiles)
        assert tiles.plot_calls[0][0] == "TILE-IMAGE"
        assert tiles.plot_calls[0][1]["ax"] is ax

    def test_the_density_colour_bar_names_the_unit(self, demography, census):
        ax = demography.presentation.plotPopulationOnMap(census, _FakeTiles(), density=True)
        assert ax.figure.axes[-1].get_ylabel().startswith("Population density [people/")

    def test_absolute_counts_are_labelled_population_count(self, demography, census):
        ax = demography.presentation.plotPopulationOnMap(census, _FakeTiles(), density=False)
        assert ax.figure.axes[-1].get_ylabel() == "Population count"

    def test_the_count_colour_bar_uses_the_raw_population_range(self, demography, census):
        ax = demography.presentation.plotPopulationOnMap(census, _FakeTiles(), density=False)
        assert ax.figure.axes[-1].get_ylim() == pytest.approx((99.0, 8000.0))

    def test_an_explicit_colour_bar_label_wins_over_both_defaults(self, demography, census):
        ax = demography.presentation.plotPopulationOnMap(
            census, _FakeTiles(), colorbar_label="mine"
        )
        assert ax.figure.axes[-1].get_ylabel() == "mine"

    def test_the_colour_bar_can_be_suppressed(self, demography, census):
        ax = demography.presentation.plotPopulationOnMap(
            census, _FakeTiles(), colorbar=False
        )
        assert len(ax.figure.axes) == 1

    def test_the_title_is_forwarded(self, demography, census):
        ax = demography.presentation.plotPopulationOnMap(
            census, _FakeTiles(), title="On map"
        )
        assert ax.get_title() == "On map"

    def test_the_input_data_is_not_mutated(self, demography, census):
        before = list(census.columns)
        demography.presentation.plotPopulationOnMap(census, _FakeTiles())
        assert list(census.columns) == before
