"""hera/measurements/GIS/vector/buildings/toolkit.py -- the three
``BuildingsToolkit`` methods that ``test_gis_buildings_toolkit_static.py``
leaves for "an integration test": ``getBuildingsFromRectangle``,
``getBuildingHeightFromRasterTopographyToolkit`` and
``buildingsGeopandasToSTLRasterTopography``.

They turn out not to need an integration environment after all.

How the inputs were produced
----------------------------
*Elevation.*  ``_write_geotiff`` synthesises the SRTM tile with ``rasterio``:
an explicit ``from_origin`` transform (origin lon 35.0 / lat 33.0, pixel
1/16 deg -- an exact binary fraction) and an explicit EPSG, written into
``tmp_path`` under the name the production code computes for itself,
``N32E035.hgt``.  The band is ``value(row, col) = 10*row + 100*col``, which is
linear and therefore reproduced *exactly* by the bilinear interpolation in
``getPointListElevation``, so the expected elevation of a building has the
closed form ``10*rasterx + 100*rastery`` rather than being a recorded number.
Building footprints are placed so their centroids land on exact fractional
pixel offsets.  The folder is registered as a real datasource on a real
(mongomock-backed) raster TopographyToolkit in the same project, and reached
through the *production* path -- the toolkit that
``getBuildingHeightFromRasterTopographyToolkit`` builds for itself via
``toolkitHome`` -- so nothing about the lookup is faked.

*Buildings.*  A four-footprint GeoDataFrame written to GeoJSON in ``tmp_path``
and registered with ``dataFormat="geopandas"`` and ``desc={"crs": ...}`` --
the nesting the repository loader produces and the one
``VectorToolkit.cutRegionFromSource`` reads (see the fixture's own note).

*STL.*  FreeCAD is not installed here.  The conftest stubs the name
``FreeCAD`` with a MagicMock, but the module's import block is
``import FreeCAD; import Part; import Mesh`` inside a single ``try`` -- so
``Part`` and ``Mesh`` are never bound at all (see B238).  The tests therefore
install small *recording* stand-ins for all three over the module globals with
``monkeypatch``, and assert on what the production code asks FreeCAD to build:
how many sketches, at what placement altitude, with what extrusion length, and
what is finally exported.  That is exactly the arithmetic and branching the
method owns; the mesh kernel itself is not under test.

Bugs pinned (each has an xfail(strict) for the intended behaviour plus a
passing characterisation of what actually happens)
--------------------------------------------------------------------------
* **B236** -- ``getBuildingHeightFromRasterTopographyToolkit`` takes a
  ``topographyDataSource`` argument, documents it as "the name of the
  datasource in the topography toolkit", and then never passes it to
  ``getPointListElevation``.  The configured default datasource is used no
  matter what the caller asks for, silently.
* **B237** -- ``buildingsGeopandasToSTLRasterTopography`` lowers each building
  by ``nonFlatTopographyShift`` and compensates by extruding
  ``height + nonFlatTopographyShift``.  In the ``flatTerrain=True`` branch the
  building is *not* lowered (the placement is ``referenceTopography``
  exactly), but the extrusion still adds the shift -- so every flat-terrain
  building comes out ``nonFlatTopographyShift`` metres too tall.
* **B238** -- the module's guarded import is
  ``try: import FreeCAD; import Part; import Mesh except ImportError``.  When
  FreeCAD is importable but ``Part`` is not (the case here, and the case for
  any partial install), ``FreeCAD`` is left bound while ``Part`` and ``Mesh``
  are not -- so the ``try: FreeCAD.newDocument(...) except: raise ValueError``
  guard that exists precisely to give a clear "install FreeCAD" message is
  bypassed, and the method dies much later with a bare ``NameError``.

Deliberately not tested, and why
--------------------------------
* ``cutRegionFromSource``'s CRS-mismatch reprojection: that is B80, already
  pinned in ``test_gis_vector_toolkit.py``, and it is not re-pinned here.  It
  turns out to be *latent* for this method:
  ``geopandas.read_file(..., mask=<gdf with a CRS>)`` reprojects the mask
  itself, so the discarded reprojection costs nothing on the ``mask=`` path.
  One passing test records that, with the caveat that it is a property of the
  reader rather than of the toolkit.
* The flat vs nested ``desc`` question: ``addDataSource(..., crs=4326)``
  stores the key flat, whereas ``cutRegionFromSource`` reads
  ``doc.desc['desc']['crs']``.  That nesting is not a defect -- the
  repository loader forwards a JSON item's own ``desc`` key into
  ``addDataSource``, so real datasources (see the ``BNTL`` entry in
  ``hera/tests/repository/testCases/REPOSITORY_TEST_01.json``) do carry it --
  so the fixture below registers the CRS the same way, rather than pinning a
  bug that is not there.
* The geometry of the STL FreeCAD would produce: no mesh kernel is available,
  and the facets are FreeCAD's responsibility, not the toolkit's.
* ``get_buildings_height`` and ``filter_buildings_in_area``: already covered
  by ``test_gis_buildings_toolkit_static.py``.
"""
import geopandas
import numpy
import pytest
from shapely.geometry import LineString, Point, box

from hera import toolkitHome
from hera.measurements.GIS.utils import ITM, WGS84
from hera.measurements.GIS.vector.buildings import toolkit as buildings_module

# --------------------------------------------------------------------------
# raster synthesis (same construction as test_gis_raster_topography_more.py)
# --------------------------------------------------------------------------

_WEST = 35.0
_NORTH = 33.0
_RES = 0.0625     # 1/16 degree: exact in binary
_SIZE = 8
_TILE = "N32E035.hgt"


def _write_geotiff(path, array, west=_WEST, north=_NORTH, resolution=_RES,
                   crs="EPSG:4326"):
    import rasterio
    from rasterio.transform import from_origin

    with rasterio.open(
        str(path), "w",
        driver="GTiff",
        height=array.shape[0], width=array.shape[1], count=1,
        dtype=array.dtype,
        transform=from_origin(west, north, resolution, resolution),
        crs=crs,
    ) as destination:
        destination.write(array, 1)
    return str(path)


def _ramp_band(size=_SIZE):
    """value(row, col) = 10*row + 100*col -- linear, hence bilinear-exact."""
    rows, cols = numpy.mgrid[0:size, 0:size]
    return (10.0 * rows + 100.0 * cols).astype("float32")


def _flat_band(value, size=_SIZE):
    return numpy.full((size, size), float(value), dtype="float32")


def _expected_elevation(lat, lon):
    rasterx = (lat - _NORTH) / -_RES
    rastery = (lon - _WEST) / _RES
    return 10.0 * rasterx + 100.0 * rastery


# Centroids on exact fractional pixel offsets: (lat, lon).
_CENTROIDS = [
    (32.96875, 35.03125),   # rasterx 0.5, rastery 0.5 ->  55.0
    (32.90625, 35.15625),   # rasterx 1.5, rastery 2.5 -> 265.0
    (32.75000, 35.06250),   # rasterx 4.0, rastery 1.0 -> 140.0
    (32.81250, 35.31250),   # rasterx 3.0, rastery 5.0 -> 530.0
]
_HALF = 0.001   # footprint half-width in degrees; keeps the centroid exact


def _buildings(crs=WGS84):
    """Four square footprints whose centroids are the points above."""
    return geopandas.GeoDataFrame(
        {
            "BLDG_HT": [10.0, 20.0, 30.0, 40.0],
            "geometry": [
                box(lon - _HALF, lat - _HALF, lon + _HALF, lat + _HALF)
                for lat, lon in _CENTROIDS
            ],
        },
        crs=crs,
    )


@pytest.fixture()
def buildings_toolkit(unit_toolkit_factory):
    return unit_toolkit_factory(toolkitHome.GIS_BUILDINGS)


@pytest.fixture()
def raster_topography(unit_toolkit_factory):
    """A raster TopographyToolkit on the same project, so config is shared."""
    return unit_toolkit_factory(toolkitHome.GIS_RASTER_TOPOGRAPHY)


@pytest.fixture()
def ramp_srtm(raster_topography, tmp_path):
    """Register the ramp tile as datasource RAMP and make it the default."""
    folder = tmp_path / "srtm_ramp"
    folder.mkdir()
    _write_geotiff(folder / _TILE, _ramp_band())
    raster_topography.addDataSource("RAMP", str(folder), "string", version=(0, 0, 1))
    raster_topography.setDataSourceDefaultVersion("RAMP", (0, 0, 1))
    raster_topography.setConfig(defaultSRTM="RAMP")
    return "RAMP"


@pytest.fixture()
def flat_srtm(raster_topography, tmp_path):
    """A second, constant-7 tile registered as FLAT but never made default."""
    folder = tmp_path / "srtm_flat"
    folder.mkdir()
    _write_geotiff(folder / _TILE, _flat_band(7.0))
    raster_topography.addDataSource("FLAT", str(folder), "string", version=(0, 0, 1))
    raster_topography.setDataSourceDefaultVersion("FLAT", (0, 0, 1))
    return "FLAT"


# ==========================================================================
# getBuildingHeightFromRasterTopographyToolkit
# ==========================================================================

@pytest.mark.unit
class TestGetBuildingHeightFromRasterTopographyToolkit:
    def test_an_elevation_column_is_joined_onto_the_buildings(self, buildings_toolkit, ramp_srtm):
        result = buildings_toolkit.getBuildingHeightFromRasterTopographyToolkit(_buildings())
        assert "elevation" in result.columns

    def test_the_original_columns_and_geometry_survive_the_join(self, buildings_toolkit, ramp_srtm):
        original = _buildings()
        result = buildings_toolkit.getBuildingHeightFromRasterTopographyToolkit(original)
        assert list(result["BLDG_HT"]) == list(original["BLDG_HT"])
        assert list(result.geometry) == list(original.geometry)

    def test_one_row_per_building_comes_back(self, buildings_toolkit, ramp_srtm):
        result = buildings_toolkit.getBuildingHeightFromRasterTopographyToolkit(_buildings())
        assert len(result) == len(_CENTROIDS)

    def test_each_elevation_is_the_raster_ramp_at_the_footprint_centroid(self, buildings_toolkit, ramp_srtm):
        """Closed form: 10*rasterx + 100*rastery on the linear band."""
        result = buildings_toolkit.getBuildingHeightFromRasterTopographyToolkit(_buildings())
        expected = [_expected_elevation(lat, lon) for lat, lon in _CENTROIDS]
        assert list(result["elevation"]) == pytest.approx(expected)

    def test_the_lat_and_lon_are_not_swapped_on_the_way_in(self, buildings_toolkit, ramp_srtm):
        """The band is asymmetric (10/row vs 100/col), so a swap would show."""
        result = buildings_toolkit.getBuildingHeightFromRasterTopographyToolkit(_buildings())
        lat, lon = _CENTROIDS[1]
        assert result["elevation"].iloc[1] == pytest.approx(_expected_elevation(lat, lon))
        assert result["elevation"].iloc[1] != pytest.approx(_expected_elevation(lon, lat))

    def test_the_centroids_are_reported_in_wgs84(self, buildings_toolkit, ramp_srtm):
        """getPointListElevation stamps lon/lat onto the frame it returns."""
        result = buildings_toolkit.getBuildingHeightFromRasterTopographyToolkit(_buildings())
        assert list(result["lat"]) == pytest.approx([lat for lat, _ in _CENTROIDS])
        assert list(result["lon"]) == pytest.approx([lon for _, lon in _CENTROIDS])

    def test_buildings_given_in_itm_are_reprojected_before_lookup(self, buildings_toolkit, ramp_srtm):
        """`.centroid.to_crs(WGS84)` must handle a projected input too."""
        in_itm = _buildings().to_crs(ITM)
        result = buildings_toolkit.getBuildingHeightFromRasterTopographyToolkit(in_itm)
        expected = [_expected_elevation(lat, lon) for lat, lon in _CENTROIDS]
        assert list(result["elevation"]) == pytest.approx(expected, rel=1e-3)

    @pytest.mark.xfail(
        strict=True,
        reason="B236: getBuildingHeightFromRasterTopographyToolkit accepts "
               "topographyDataSource, documents it as the datasource to read "
               "the topography from, and never forwards it to "
               "getPointListElevation -- the configured default is always "
               "used. See the consolidated findings issue.",
    )
    def test_the_requested_topography_datasource_is_used(self, buildings_toolkit, ramp_srtm, flat_srtm):
        result = buildings_toolkit.getBuildingHeightFromRasterTopographyToolkit(
            _buildings(), topographyDataSource=flat_srtm
        )
        assert list(result["elevation"]) == pytest.approx([7.0] * len(_CENTROIDS))

    def test_the_topography_datasource_argument_is_ignored(self, buildings_toolkit, ramp_srtm, flat_srtm):
        """Characterisation of B236: the default wins even when another is named."""
        result = buildings_toolkit.getBuildingHeightFromRasterTopographyToolkit(
            _buildings(), topographyDataSource=flat_srtm
        )
        expected = [_expected_elevation(lat, lon) for lat, lon in _CENTROIDS]
        assert list(result["elevation"]) == pytest.approx(expected)

    def test_even_a_datasource_that_does_not_exist_is_ignored(self, buildings_toolkit, ramp_srtm):
        """Characterisation of B236: no validation, no error, no effect."""
        result = buildings_toolkit.getBuildingHeightFromRasterTopographyToolkit(
            _buildings(), topographyDataSource="NO_SUCH_DATASOURCE"
        )
        assert "elevation" in result.columns

    def test_a_missing_default_datasource_is_reported(self, buildings_toolkit, raster_topography):
        """No defaultSRTM configured -> the topography toolkit complains."""
        with pytest.raises((ValueError, KeyError)):
            buildings_toolkit.getBuildingHeightFromRasterTopographyToolkit(_buildings())


# ==========================================================================
# getBuildingsFromRectangle
# ==========================================================================

@pytest.fixture()
def buildings_source(buildings_toolkit, tmp_path):
    """Register the four footprints as a geopandas datasource in WGS84.

    The ``desc=dict(crs=...)`` nesting is not decoration: it is the shape the
    repository loader produces (``_handle_DataSource`` forwards the JSON item's
    own ``desc`` key straight into ``addDataSource``, so it lands one level
    down), and ``VectorToolkit.cutRegionFromSource`` reads exactly
    ``doc.desc['desc']['crs']``.  Registering the CRS flat instead -- as
    ``addDataSource(..., crs=...)`` -- puts it where nothing looks for it.
    Compare hera/tests/repository/testCases/REPOSITORY_TEST_01.json, whose
    ``BNTL`` entry carries ``"desc": {"crs": 2039, ...}``.
    """
    path = tmp_path / "buildings.geojson"
    _buildings().to_file(str(path), driver="GeoJSON")
    buildings_toolkit.addDataSource(
        "BNTL", str(path), "geopandas", version=(0, 0, 1), desc=dict(crs=WGS84)
    )
    buildings_toolkit.setDataSourceDefaultVersion("BNTL", (0, 0, 1))
    return "BNTL"


def _rectangle_around(index, pad=0.01):
    lat, lon = _CENTROIDS[index]
    return dict(minx=lon - pad, miny=lat - pad, maxx=lon + pad, maxy=lat + pad)


@pytest.mark.unit
class TestGetBuildingsFromRectangle:
    def test_it_returns_a_geodataframe(self, buildings_toolkit, buildings_source):
        result = buildings_toolkit.getBuildingsFromRectangle(
            dataSourceName=buildings_source, **_rectangle_around(0)
        )
        assert isinstance(result, geopandas.GeoDataFrame)

    def test_a_tight_rectangle_returns_only_the_footprint_inside_it(self, buildings_toolkit, buildings_source):
        result = buildings_toolkit.getBuildingsFromRectangle(
            dataSourceName=buildings_source, **_rectangle_around(0)
        )
        assert len(result) == 1
        assert result["BLDG_HT"].iloc[0] == 10.0

    def test_each_footprint_can_be_selected_in_turn(self, buildings_toolkit, buildings_source):
        heights = []
        for index in range(len(_CENTROIDS)):
            result = buildings_toolkit.getBuildingsFromRectangle(
                dataSourceName=buildings_source, **_rectangle_around(index)
            )
            heights.append(result["BLDG_HT"].iloc[0])
        assert heights == [10.0, 20.0, 30.0, 40.0]

    def test_a_rectangle_covering_everything_returns_every_footprint(self, buildings_toolkit, buildings_source):
        result = buildings_toolkit.getBuildingsFromRectangle(
            minx=34.9, miny=32.5, maxx=35.5, maxy=33.0, dataSourceName=buildings_source
        )
        assert len(result) == len(_CENTROIDS)

    def test_a_rectangle_that_misses_everything_returns_nothing(self, buildings_toolkit, buildings_source):
        result = buildings_toolkit.getBuildingsFromRectangle(
            minx=20.0, miny=20.0, maxx=21.0, maxy=21.0, dataSourceName=buildings_source
        )
        assert len(result) == 0

    def test_the_default_datasource_comes_from_the_config(self, buildings_toolkit, buildings_source):
        buildings_toolkit.setConfig(defaultBuildingDataSource=buildings_source)
        result = buildings_toolkit.getBuildingsFromRectangle(**_rectangle_around(2))
        assert result["BLDG_HT"].iloc[0] == 30.0

    def test_a_missing_config_key_is_reported(self, buildings_toolkit, buildings_source):
        with pytest.raises(KeyError, match="defaultBuildingDataSource"):
            buildings_toolkit.getBuildingsFromRectangle(**_rectangle_around(0))

    def test_elevation_is_not_looked_up_unless_asked_for(self, buildings_toolkit, buildings_source):
        result = buildings_toolkit.getBuildingsFromRectangle(
            dataSourceName=buildings_source, **_rectangle_around(0)
        )
        assert "elevation" not in result.columns

    def test_a_rectangle_in_another_crs_still_selects_the_right_footprint(self, buildings_toolkit, buildings_source):
        """B80 is survivable here, and it is worth recording why.

        B80 (pinned in ``test_gis_vector_toolkit.py``) is that
        ``VectorToolkit.cutRegionFromSource`` builds its ``mask=`` dict before
        the CRS check and then rebinds the name, so the reprojected frame
        never reaches ``getData``.  A caller passing ``inputCRS=ITM`` against
        a WGS84 layer therefore hands an *ITM* mask to
        ``geopandas.read_file(..., mask=...)`` -- which reprojects a
        CRS-carrying mask to the dataset's CRS itself.  The right footprint
        comes back regardless, so B80 is latent for this method rather than
        user-visible; that is a property of the reader, not of the toolkit,
        and it would not save the ``bbox=`` path if the reader ever stopped.
        """
        lat, lon = _CENTROIDS[1]
        centre = geopandas.GeoSeries([Point(lon, lat)], crs=WGS84).to_crs(ITM).iloc[0]
        result = buildings_toolkit.getBuildingsFromRectangle(
            minx=centre.x - 500, miny=centre.y - 500,
            maxx=centre.x + 500, maxy=centre.y + 500,
            dataSourceName=buildings_source, inputCRS=ITM,
        )
        assert list(result["BLDG_HT"]) == [20.0]

    def test_with_elevation_joins_the_topography_height(self, buildings_toolkit, buildings_source, ramp_srtm):
        result = buildings_toolkit.getBuildingsFromRectangle(
            dataSourceName=buildings_source, withElevation=True, **_rectangle_around(1)
        )
        lat, lon = _CENTROIDS[1]
        assert result["elevation"].iloc[0] == pytest.approx(_expected_elevation(lat, lon))


# ==========================================================================
# buildingsGeopandasToSTLRasterTopography
# ==========================================================================

class _Vector:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z

    def __repr__(self):  # pragma: no cover - debugging aid only
        return f"_Vector({self.x}, {self.y}, {self.z})"


class _Rotation:
    def __init__(self, *components):
        self.components = components


class _Placement:
    def __init__(self, vector, rotation):
        self.vector = vector
        self.rotation = rotation


class _Sketch:
    kind = "Sketcher::SketchObject"

    def __init__(self, name):
        self.name = name
        self.Placement = None
        self.segments = []

    def addGeometry(self, segment):
        self.segments.append(segment)


class _Pad:
    kind = "Part::Extrusion"

    def __init__(self, name):
        self.name = name
        self.Base = None
        self.LengthFwd = None
        self.Solid = None
        self.Symmetric = None


class _Document:
    def __init__(self):
        self.Objects = []
        self.recomputes = 0

    def addObject(self, kind, name):
        obj = _Sketch(name) if kind == _Sketch.kind else _Pad(name)
        self.Objects.append(obj)
        return obj

    def recompute(self):
        self.recomputes += 1

    @property
    def sketches(self):
        return [o for o in self.Objects if isinstance(o, _Sketch)]

    @property
    def pads(self):
        return [o for o in self.Objects if isinstance(o, _Pad)]


class _FreeCADStub:
    Vector = _Vector
    Rotation = _Rotation
    Placement = _Placement

    def __init__(self):
        self.documents = []

    def newDocument(self, name):
        document = _Document()
        document.name = name
        self.documents.append(document)
        return document


class _PartStub:
    @staticmethod
    def LineSegment(start, end):
        return (start, end)


class _MeshStub:
    def __init__(self):
        self.exports = []

    def export(self, objects, filename):
        self.exports.append((list(objects), filename))


@pytest.fixture()
def freecad(monkeypatch):
    """Recording stand-ins for FreeCAD, Part and Mesh on the module globals."""
    stub = _FreeCADStub()
    mesh = _MeshStub()
    monkeypatch.setattr(buildings_module, "FreeCAD", stub, raising=False)
    monkeypatch.setattr(buildings_module, "Part", _PartStub, raising=False)
    monkeypatch.setattr(buildings_module, "Mesh", mesh, raising=False)
    stub.mesh = mesh
    return stub


def _stl_input():
    """Buildings with a height column and an elevation column, in ITM."""
    return geopandas.GeoDataFrame(
        {
            "BLDG_HT": [10.0, 25.0],
            "HT_LAND": [100.0, 250.0],
            "geometry": [box(0, 0, 10, 10), box(20, 20, 35, 35)],
        },
        crs=ITM,
    )


def _call(toolkit, data=None, **kwargs):
    data = _stl_input() if data is None else data
    kwargs.setdefault("outputFileName", "/tmp/unused-buildings.stl")
    return toolkit.buildingsGeopandasToSTLRasterTopography(
        buildingData=data,
        buildingHeightColumn="BLDG_HT",
        buildingElevationColumn="HT_LAND",
        **kwargs,
    )


@pytest.mark.unit
class TestBuildingsGeopandasToSTLModuleImports:
    def test_freecad_is_bound_but_part_and_mesh_are_not(self):
        """Characterisation of B238: the guarded import binds only the first name."""
        assert hasattr(buildings_module, "FreeCAD")
        assert not hasattr(buildings_module, "Part")
        assert not hasattr(buildings_module, "Mesh")

    @pytest.mark.xfail(
        strict=True,
        reason="B238: the module's 'try: import FreeCAD; import Part; import "
               "Mesh' leaves FreeCAD bound when only Part/Mesh are missing, "
               "so the method's own 'FreeCAD not found. Install before using "
               "this function' guard never fires and the failure surfaces "
               "later as a bare NameError on Part. See the consolidated "
               "findings issue.",
    )
    def test_a_partial_freecad_install_is_reported_as_such(self, buildings_toolkit):
        with pytest.raises(ValueError, match="FreeCAD not found"):
            _call(buildings_toolkit)

    def test_a_partial_install_dies_on_an_unbound_name_instead(self, buildings_toolkit):
        """Characterisation of B238."""
        with pytest.raises(NameError) as excinfo:
            _call(buildings_toolkit)
        assert "Part" in str(excinfo.value)


@pytest.mark.unit
class TestBuildingsGeopandasToSTLStructure:
    def test_one_sketch_and_one_pad_are_created_per_building(self, buildings_toolkit, freecad):
        _call(buildings_toolkit)
        document = freecad.documents[0]
        assert len(document.sketches) == 2
        assert len(document.pads) == 2

    def test_the_objects_are_named_after_the_dataframe_index(self, buildings_toolkit, freecad):
        _call(buildings_toolkit)
        document = freecad.documents[0]
        assert [s.name for s in document.sketches] == ["Sketch0", "Sketch1"]
        assert [p.name for p in document.pads] == ["building0", "building1"]

    def test_each_pad_is_based_on_its_own_sketch(self, buildings_toolkit, freecad):
        _call(buildings_toolkit)
        document = freecad.documents[0]
        for sketch, pad in zip(document.sketches, document.pads):
            assert pad.Base is sketch

    def test_every_pad_is_a_non_symmetric_solid(self, buildings_toolkit, freecad):
        _call(buildings_toolkit)
        for pad in freecad.documents[0].pads:
            assert pad.Solid is True
            assert pad.Symmetric is False

    def test_the_document_is_recomputed_once_per_building(self, buildings_toolkit, freecad):
        _call(buildings_toolkit)
        assert freecad.documents[0].recomputes == 2

    def test_one_wall_segment_is_added_per_polygon_edge(self, buildings_toolkit, freecad):
        """A closed 10x10 box has five exterior coordinates, hence four edges."""
        _call(buildings_toolkit)
        assert len(freecad.documents[0].sketches[0].segments) == 4

    def test_the_wall_segments_trace_the_footprint_outline(self, buildings_toolkit, freecad):
        _call(buildings_toolkit)
        segments = freecad.documents[0].sketches[0].segments
        drawn = [(seg[0].x, seg[0].y) for seg in segments]
        expected = list(box(0, 0, 10, 10).exterior.coords)[:-1]
        assert set(drawn) == set(expected)

    def test_the_export_receives_every_object_and_the_requested_path(self, buildings_toolkit, freecad, tmp_path):
        target = str(tmp_path / "buildings.stl")
        _call(buildings_toolkit, outputFileName=target)
        objects, filename = freecad.mesh.exports[0]
        assert filename == target
        assert objects == freecad.documents[0].Objects

    def test_rows_whose_geometry_has_no_exterior_are_skipped(self, buildings_toolkit, freecad):
        data = geopandas.GeoDataFrame(
            {
                "BLDG_HT": [10.0, 10.0, 10.0],
                "HT_LAND": [100.0, 100.0, 100.0],
                "geometry": [box(0, 0, 5, 5), Point(1, 1), LineString([(0, 0), (1, 1)])],
            },
            crs=ITM,
        )
        _call(buildings_toolkit, data=data)
        assert len(freecad.documents[0].sketches) == 1

    def test_an_empty_frame_still_exports_an_empty_document(self, buildings_toolkit, freecad):
        empty = _stl_input().iloc[0:0]
        _call(buildings_toolkit, data=empty)
        assert freecad.mesh.exports[0][0] == []

    def test_it_returns_nothing_and_communicates_only_through_the_file(self, buildings_toolkit, freecad):
        assert _call(buildings_toolkit) is None


@pytest.mark.unit
class TestBuildingsGeopandasToSTLAltitudes:
    def test_on_real_terrain_the_base_sits_below_the_ground_by_the_shift(self, buildings_toolkit, freecad):
        _call(buildings_toolkit, flatTerrain=False, nonFlatTopographyShift=10)
        placements = [s.Placement.vector.z for s in freecad.documents[0].sketches]
        assert placements == pytest.approx([100.0 - 10.0, 250.0 - 10.0])

    def test_the_shift_is_configurable(self, buildings_toolkit, freecad):
        _call(buildings_toolkit, flatTerrain=False, nonFlatTopographyShift=3)
        placements = [s.Placement.vector.z for s in freecad.documents[0].sketches]
        assert placements == pytest.approx([97.0, 247.0])

    def test_on_real_terrain_the_building_top_lands_on_ground_plus_height(self, buildings_toolkit, freecad):
        """base + LengthFwd must equal elevation + building height."""
        _call(buildings_toolkit, flatTerrain=False, nonFlatTopographyShift=10)
        document = freecad.documents[0]
        tops = [
            sketch.Placement.vector.z + pad.LengthFwd
            for sketch, pad in zip(document.sketches, document.pads)
        ]
        assert tops == pytest.approx([100.0 + 10.0, 250.0 + 25.0])

    def test_on_flat_terrain_the_base_is_the_reference_height(self, buildings_toolkit, freecad):
        _call(buildings_toolkit, flatTerrain=True, referenceTopography=5.0)
        placements = [s.Placement.vector.z for s in freecad.documents[0].sketches]
        assert placements == pytest.approx([5.0, 5.0])

    def test_on_flat_terrain_the_elevation_column_is_not_read(self, buildings_toolkit, freecad):
        """A frame without the elevation column must still convert."""
        data = _stl_input().drop(columns=["HT_LAND"])
        _call(buildings_toolkit, data=data, flatTerrain=True, referenceTopography=0)
        assert len(freecad.documents[0].sketches) == 2

    def test_the_reference_topography_defaults_to_zero(self, buildings_toolkit, freecad):
        _call(buildings_toolkit, flatTerrain=True)
        placements = [s.Placement.vector.z for s in freecad.documents[0].sketches]
        assert placements == pytest.approx([0.0, 0.0])

    def test_the_placement_is_only_vertical(self, buildings_toolkit, freecad):
        """The footprint carries the x/y position; the placement only lifts it."""
        _call(buildings_toolkit, flatTerrain=True, referenceTopography=5.0)
        for sketch in freecad.documents[0].sketches:
            assert (sketch.Placement.vector.x, sketch.Placement.vector.y) == (0.0, 0.0)

    @pytest.mark.xfail(
        strict=True,
        reason="B237: buildingsGeopandasToSTLRasterTopography extrudes "
               "'wallsheight + nonFlatTopographyShift' unconditionally, but "
               "the flatTerrain branch places the sketch at "
               "referenceTopography without lowering it by the shift -- so "
               "every flat-terrain building is nonFlatTopographyShift metres "
               "too tall. See the consolidated findings issue.",
    )
    def test_on_flat_terrain_the_building_top_is_reference_plus_height(self, buildings_toolkit, freecad):
        _call(buildings_toolkit, flatTerrain=True, referenceTopography=0.0,
              nonFlatTopographyShift=10)
        document = freecad.documents[0]
        tops = [
            sketch.Placement.vector.z + pad.LengthFwd
            for sketch, pad in zip(document.sketches, document.pads)
        ]
        assert tops == pytest.approx([10.0, 25.0])

    def test_flat_terrain_buildings_are_too_tall_by_the_shift(self, buildings_toolkit, freecad):
        """Characterisation of B237."""
        _call(buildings_toolkit, flatTerrain=True, referenceTopography=0.0,
              nonFlatTopographyShift=10)
        document = freecad.documents[0]
        tops = [
            sketch.Placement.vector.z + pad.LengthFwd
            for sketch, pad in zip(document.sketches, document.pads)
        ]
        assert tops == pytest.approx([10.0 + 10.0, 25.0 + 10.0])

    def test_the_extrusion_length_ignores_the_terrain_flag_entirely(self, buildings_toolkit, freecad):
        """Characterisation of B237: identical LengthFwd in both branches."""
        _call(buildings_toolkit, flatTerrain=True, nonFlatTopographyShift=10)
        flat = [pad.LengthFwd for pad in freecad.documents[0].pads]
        _call(buildings_toolkit, flatTerrain=False, nonFlatTopographyShift=10)
        real = [pad.LengthFwd for pad in freecad.documents[1].pads]
        assert flat == real == pytest.approx([20.0, 35.0])

    def test_a_zero_shift_makes_the_two_branches_agree(self, buildings_toolkit, freecad):
        """Characterisation of B237: nonFlatTopographyShift=0 is the safe case."""
        _call(buildings_toolkit, flatTerrain=True, referenceTopography=0.0,
              nonFlatTopographyShift=0)
        document = freecad.documents[0]
        tops = [
            sketch.Placement.vector.z + pad.LengthFwd
            for sketch, pad in zip(document.sketches, document.pads)
        ]
        assert tops == pytest.approx([10.0, 25.0])
