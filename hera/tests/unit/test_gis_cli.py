"""measurements/GIS/CLI.py: the lazy-import bootstrap plus the six argparse
command handlers.

Covered: ``_setup``, ``_lazy_setup`` (and its ``wrapper``),
``topography_vector_list``, ``topography_raster_list``,
``topography_raster_toSTL``, ``buildings_parser_list``,
``buildings_raster_toSTL`` and ``get_landocver``.

Every handler is exercised against a fake toolkitHome injected over the
module global that ``_setup`` would otherwise fill in (with ``_setup``
itself stubbed out, which is also what makes the handler tests survive
B126 below), so the tests assert on the CLI plumbing -- which toolkit was
asked for, which arguments were forwarded, what was printed, what landed
on disk -- and never touch the GIS algorithms, GDAL or a database.  One
test checks the real ``toolkitHome`` still exposes the four constants the
fake stands in for, so the fake cannot silently drift.

Deliberately not tested: the numeric content of the STL string and of the
landcover xarray (both belong to the toolkits, which have their own unit
tests), and the argparse wiring in ``hera/bin/hera-GIS``, which lives under
``if __name__ == "__main__"`` and cannot be imported.

Three bugs are pinned here:

* B126 -- ``_setup`` imports ``WSG84``/``ITM`` from ``hera.utils``, which
  does not define them, so every hera-GIS command dies with ImportError.
* B127 -- ``get_landocver`` misspells ``windDirection`` as
  ``windDirectionis``.
* B128 -- ``buildings_raster_toSTL`` computes ``dxdy``/``inputCRS``/
  ``outputCRS``/``dataSourceName`` and then forwards none of them.
"""
from argparse import Namespace

import numpy
import pandas
import pytest

from hera.measurements.GIS import CLI


# ---------------------------------------------------------------------------
# fakes
# ---------------------------------------------------------------------------

class _FakeVariable:
    """Stands in for an xarray DataArray: only ``.values`` is ever read."""

    def __init__(self, values):
        self.values = numpy.array(values)


class _FakeXarray:
    def __init__(self, fields, coords=()):
        self._fields = {name: _FakeVariable(values) for name, values in fields.items()}
        self.coords = list(coords)

    def __getitem__(self, name):
        return self._fields[name]


STL_TEXT = "solid Topography\nendsolid Topography\n"


def _plainLandcover():
    return _FakeXarray(
        fields=dict(lat=[[0.0, 1.0]], lon=[[2.0, 3.0]], landcover=[[10, 20]]),
    )


class _FakeToolkit:
    def __init__(self, projectName, landcover=None):
        self.projectName = projectName if projectName is not None else "defaultProject"
        self.calls = []
        self.landcover = landcover if landcover is not None else _plainLandcover()

    # -- listing ---------------------------------------------------------
    def getDataSourceTable(self):
        self.calls.append(("getDataSourceTable", {}))
        return "|datasource|version|"

    # -- topography ------------------------------------------------------
    def createElevationSTL(self, **kwargs):
        self.calls.append(("createElevationSTL", kwargs))
        return STL_TEXT

    # -- buildings -------------------------------------------------------
    def getBuildingsFromRectangle(self, **kwargs):
        self.calls.append(("getBuildingsFromRectangle", kwargs))
        return "<buildings geopandas>"

    def buildingsGeopandasToSTLRasterTopography(self, *args):
        self.calls.append(("buildingsGeopandasToSTLRasterTopography", args))

    # -- landcover -------------------------------------------------------
    def getRoughness(self, *args, **kwargs):
        self.calls.append(("getRoughness", args, kwargs))
        return self.landcover

    def getLandCover(self, *args, **kwargs):
        self.calls.append(("getLandCover", args, kwargs))
        return self.landcover

    def call(self, name):
        """The single recorded call named ``name``."""
        matches = [entry for entry in self.calls if entry[0] == name]
        assert len(matches) == 1, f"expected exactly one {name} call, got {self.calls}"
        return matches[0]


class _FakeToolkitHome:
    """Records every getToolkit request and hands back one fake toolkit."""

    GIS_VECTOR_TOPOGRAPHY = "GIS_Vector_Topography"
    GIS_RASTER_TOPOGRAPHY = "GIS_Raster_Topography"
    GIS_BUILDINGS = "GIS_Buildings"
    GIS_LANDCOVER = "GIS_Raster_Landcover"

    def __init__(self):
        self.requests = []
        self.toolkit = None
        self.landcover = None  # the xarray every toolkit it builds will return

    def getToolkit(self, toolkitName, projectName=None, **kwargs):
        self.requests.append(dict(toolkitName=toolkitName, projectName=projectName, **kwargs))
        self.toolkit = _FakeToolkit(projectName, landcover=self.landcover)
        return self.toolkit


@pytest.fixture()
def home(monkeypatch):
    """Neutralise the deferred imports and inject the fake toolkitHome.

    ``_lazy_setup``'s wrapper looks ``_setup`` up in the module globals at
    call time, so replacing it here keeps the real pandas/hera imports from
    overwriting the fakes below.
    """
    fake = _FakeToolkitHome()
    monkeypatch.setattr(CLI, "_setup", lambda: None)
    monkeypatch.setattr(CLI, "toolkitHome", fake)
    monkeypatch.setattr(CLI, "pd", pandas)
    monkeypatch.setattr(CLI, "WSG84", 4326)
    monkeypatch.setattr(CLI, "ITM", 2039)
    return fake


def _landcoverArguments(tmp_path, **overrides):
    arguments = Namespace(
        minx="1", miny="2", maxx="3", maxy="4",
        dxdy="50",
        isBuilding=False,
        windDirection=None,
        resolution=None,
        roughness=True,
        inputCRS=None,
        outputCRS=None,
        projectName="LANDCOVER_PROJECT",
        dataSourceName=None,
        filePath=str(tmp_path / "landcover.csv"),
    )
    for key, value in overrides.items():
        setattr(arguments, key, value)
    return arguments


# ---------------------------------------------------------------------------
# the bootstrap
# ---------------------------------------------------------------------------

@pytest.fixture()
def uninitialized(monkeypatch):
    """Rewind the module to its pre-``_setup`` state, and restore it after."""
    monkeypatch.setattr(CLI, "_initialized", False)
    monkeypatch.setattr(CLI, "pd", None)
    monkeypatch.setattr(CLI, "toolkitHome", None)
    monkeypatch.setattr(CLI, "WSG84", None)
    monkeypatch.setattr(CLI, "ITM", None)


@pytest.mark.unit
class TestSetup:
    def test_it_binds_pandas_and_the_toolkit_home(self, uninitialized):
        with pytest.raises(ImportError):  # B126, pinned below
            CLI._setup()
        assert CLI.pd is pandas
        assert hasattr(CLI.toolkitHome, "getToolkit")
        assert CLI._initialized is True

    def test_the_real_toolkit_home_exposes_every_constant_the_handlers_use(self):
        """Guards the fake above: a renamed constant must break a test."""
        from hera import toolkitHome as realToolkitHome

        for name in ("GIS_VECTOR_TOPOGRAPHY", "GIS_RASTER_TOPOGRAPHY",
                     "GIS_BUILDINGS", "GIS_LANDCOVER"):
            assert hasattr(realToolkitHome, name)
            assert hasattr(_FakeToolkitHome, name)

    def test_a_second_call_does_not_rebind_anything(self, monkeypatch):
        monkeypatch.setattr(CLI, "_initialized", True)
        monkeypatch.setattr(CLI, "pd", "sentinel")
        CLI._setup()
        assert CLI.pd == "sentinel"

    @pytest.mark.xfail(
        strict=True,
        reason="B126: _setup does `from hera.utils import WSG84, ITM`, but "
               "those constants live in hera.measurements.GIS.utils and are "
               "not re-exported by hera.utils (whose lazy __getattr__ only "
               "searches unitHandler/jsonutils/query/matplotlibCountour/"
               "angle/zipUtils). _setup therefore raises ImportError, and "
               "since every command is wrapped in @_lazy_setup the whole "
               "hera-GIS CLI is dead on the first call. See the consolidated "
               "findings issue.",
    )
    def test_it_should_bind_the_crs_constants(self, uninitialized):
        CLI._setup()
        assert CLI.WSG84 == 4326
        assert CLI.ITM == 2039

    def test_it_currently_raises_importerror_for_the_crs_constants(self, uninitialized):
        """Characterisation of B126."""
        with pytest.raises(ImportError, match="WSG84"):
            CLI._setup()
        assert CLI.WSG84 is None
        assert CLI.ITM is None

    def test_the_crs_constants_really_live_in_the_gis_utils_module(self):
        """Characterisation of B126: where the import should have pointed."""
        from hera.measurements.GIS.utils import ITM, WSG84

        assert (WSG84, ITM) == (4326, 2039)
        with pytest.raises(ImportError):
            from hera.utils import WSG84 as _unused  # noqa: F401


@pytest.mark.unit
class TestEveryCommandIsDeadOnArrival:
    @pytest.mark.xfail(
        strict=True,
        reason="B126: the @_lazy_setup bootstrap raises ImportError on "
               "hera.utils' missing WSG84, so no hera-GIS command can run. "
               "See the consolidated findings issue.",
    )
    def test_listing_the_datasources_should_work(self, uninitialized):
        CLI.topography_raster_list(Namespace(projectName=None))

    def test_the_first_call_raises_importerror(self, uninitialized):
        """Characterisation of B126."""
        with pytest.raises(ImportError, match="WSG84"):
            CLI.topography_raster_list(Namespace(projectName=None))

    def test_the_bootstrap_marks_itself_done_even_though_it_failed(self, uninitialized):
        """Characterisation of B126: _setup sets _initialized = True before
        it performs the imports, so the failure is never retried -- the
        second call returns silently with WSG84/ITM still unbound, which
        turns a loud ImportError into a silent wrong-value hazard."""
        with pytest.raises(ImportError):
            CLI.topography_raster_list(Namespace(projectName=None))

        CLI._setup()  # no exception this time

        assert CLI._initialized is True
        assert CLI.WSG84 is None and CLI.ITM is None


@pytest.mark.unit
class TestLazySetup:
    def test_the_wrapper_runs_setup_before_the_handler(self, monkeypatch):
        events = []
        monkeypatch.setattr(CLI, "_setup", lambda: events.append("setup"))

        @CLI._lazy_setup
        def probe(arguments):
            events.append(("probe", arguments))
            return "returned"

        assert probe("ARGS") == "returned"
        assert events == ["setup", ("probe", "ARGS")]

    def test_it_forwards_positional_and_keyword_arguments(self, monkeypatch):
        monkeypatch.setattr(CLI, "_setup", lambda: None)

        @CLI._lazy_setup
        def probe(*args, **kwargs):
            return args, kwargs

        assert probe(1, 2, key="value") == ((1, 2), {"key": "value"})

    def test_the_decorated_handlers_keep_their_own_identity(self):
        assert CLI.topography_raster_toSTL.__name__ == "topography_raster_toSTL"
        assert "STL" in CLI.topography_raster_toSTL.__doc__
        assert CLI.get_landocver.__wrapped__.__name__ == "get_landocver"


# ---------------------------------------------------------------------------
# the listing handlers
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestListingHandlers:
    @pytest.mark.parametrize("handler,toolkitName,header", [
        ("topography_vector_list", "GIS_VECTOR_TOPOGRAPHY", "topography (vector)"),
        ("topography_raster_list", "GIS_RASTER_TOPOGRAPHY", "topography (raster)"),
        ("buildings_parser_list", "GIS_BUILDINGS", "buildings"),
    ])
    def test_it_prints_a_header_and_the_datasource_table(self, home, capsys,
                                                         handler, toolkitName, header):
        getattr(CLI, handler)(Namespace(projectName="MY_PROJECT"))

        assert home.requests == [dict(toolkitName=getattr(_FakeToolkitHome, toolkitName),
                                      projectName="MY_PROJECT")]
        printed = capsys.readouterr().out.splitlines()
        assert printed[0] == f"Loaded {header} datasources for the project MY_PROJECT"
        assert printed[1] == "|datasource|version|"
        assert home.toolkit.call("getDataSourceTable")

    @pytest.mark.parametrize("handler", [
        "topography_vector_list", "topography_raster_list", "buildings_parser_list",
    ])
    def test_without_a_projectname_argument_the_default_project_is_used(self, home, capsys, handler):
        getattr(CLI, handler)(Namespace())

        assert home.requests[0]["projectName"] is None
        assert "defaultProject" in capsys.readouterr().out

    def test_an_explicit_none_projectname_is_passed_through_as_none(self, home):
        CLI.topography_raster_list(Namespace(projectName=None))
        assert home.requests[0]["projectName"] is None


# ---------------------------------------------------------------------------
# topography -> STL
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestTopographyRasterToSTL:
    def _arguments(self, tmp_path, **overrides):
        arguments = Namespace(
            minx=10, miny=20, maxx=30, maxy=40,
            dxdy=15,
            inputCRS=None,
            dataSourceName=None,
            solidName=None,
            fileName=str(tmp_path / "topography"),
            projectName="STL_PROJECT",
        )
        for key, value in overrides.items():
            setattr(arguments, key, value)
        return arguments

    def test_it_writes_the_toolkit_stl_string_to_the_requested_file(self, home, tmp_path):
        CLI.topography_raster_toSTL(self._arguments(tmp_path))
        assert (tmp_path / "topography.stl").read_text() == STL_TEXT

    def test_the_stl_suffix_is_appended_only_when_it_is_missing(self, home, tmp_path):
        CLI.topography_raster_toSTL(self._arguments(tmp_path, fileName=str(tmp_path / "already.stl")))
        assert (tmp_path / "already.stl").read_text() == STL_TEXT
        assert not (tmp_path / "already.stl.stl").exists()

    def test_a_filename_that_merely_contains_dot_stl_is_left_untouched(self, home, tmp_path):
        """The check is a substring test, not a suffix test."""
        CLI.topography_raster_toSTL(self._arguments(tmp_path, fileName=str(tmp_path / "keep.stl.bak")))
        assert (tmp_path / "keep.stl.bak").read_text() == STL_TEXT

    def test_it_asks_the_raster_topography_toolkit_for_the_project(self, home, tmp_path):
        CLI.topography_raster_toSTL(self._arguments(tmp_path))
        assert home.requests == [dict(toolkitName=_FakeToolkitHome.GIS_RASTER_TOPOGRAPHY,
                                      projectName="STL_PROJECT")]

    def test_the_bounding_box_and_resolution_reach_the_toolkit(self, home, tmp_path):
        CLI.topography_raster_toSTL(self._arguments(tmp_path))
        _, kwargs = home.toolkit.call("createElevationSTL")
        assert kwargs["minx"] == 10 and kwargs["miny"] == 20
        assert kwargs["maxx"] == 30 and kwargs["maxy"] == 40
        assert kwargs["dxdy"] == 15

    def test_a_missing_inputcrs_falls_back_to_wsg84(self, home, tmp_path):
        CLI.topography_raster_toSTL(self._arguments(tmp_path))
        assert home.toolkit.call("createElevationSTL")[1]["inputCRS"] == 4326

    def test_an_explicit_inputcrs_is_forwarded_unchanged(self, home, tmp_path):
        CLI.topography_raster_toSTL(self._arguments(tmp_path, inputCRS=2039))
        assert home.toolkit.call("createElevationSTL")[1]["inputCRS"] == 2039

    def test_a_missing_solidname_defaults_to_topography(self, home, tmp_path):
        CLI.topography_raster_toSTL(self._arguments(tmp_path))
        assert home.toolkit.call("createElevationSTL")[1]["solidName"] == "Topography"

    def test_an_explicit_solidname_and_datasource_are_forwarded(self, home, tmp_path):
        CLI.topography_raster_toSTL(
            self._arguments(tmp_path, solidName="Hill", dataSourceName="SRTM")
        )
        kwargs = home.toolkit.call("createElevationSTL")[1]
        assert kwargs["solidName"] == "Hill"
        assert kwargs["dataSourceName"] == "SRTM"

    def test_without_a_dxdy_argument_the_default_resolution_is_used(self, home, tmp_path):
        arguments = self._arguments(tmp_path)
        del arguments.dxdy
        CLI.topography_raster_toSTL(arguments)
        assert home.toolkit.call("createElevationSTL")[1]["dxdy"] == 30


# ---------------------------------------------------------------------------
# buildings -> STL
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestBuildingsRasterToSTL:
    def _arguments(self, tmp_path, **overrides):
        arguments = Namespace(
            minx=1, miny=2, maxx=3, maxy=4,
            dxdy=5,
            inputCRS=None,
            outputCRS=None,
            dataSourceName=None,
            solidName="Topography",
            fileName=str(tmp_path / "buildings"),
        )
        for key, value in overrides.items():
            setattr(arguments, key, value)
        return arguments

    def test_it_asks_the_buildings_toolkit_for_the_rectangle_with_elevation(self, home, tmp_path):
        CLI.buildings_raster_toSTL(self._arguments(tmp_path))

        assert home.requests[0]["toolkitName"] == _FakeToolkitHome.GIS_BUILDINGS
        _, kwargs = home.toolkit.call("getBuildingsFromRectangle")
        assert kwargs["minx"] == 1 and kwargs["maxy"] == 4
        assert kwargs["withElevation"] is True

    def test_it_hands_the_buildings_the_height_and_elevation_columns(self, home, tmp_path):
        CLI.buildings_raster_toSTL(self._arguments(tmp_path))
        _, args = home.toolkit.call("buildingsGeopandasToSTLRasterTopography")
        assert args == ("<buildings geopandas>", "BLDG_HT", "elevation",
                        str(tmp_path / "buildings.stl"))

    def test_the_stl_suffix_is_appended_before_the_toolkit_writes(self, home, tmp_path):
        CLI.buildings_raster_toSTL(self._arguments(tmp_path, fileName=str(tmp_path / "town.stl")))
        _, args = home.toolkit.call("buildingsGeopandasToSTLRasterTopography")
        assert args[-1] == str(tmp_path / "town.stl")

    def test_this_command_has_no_projectname_so_the_default_project_is_used(self, home, tmp_path):
        """hera-GIS never defines --projectName for `buildings toSTL`; the
        handler's `"projectName" in arguments` guard turns that into None."""
        CLI.buildings_raster_toSTL(self._arguments(tmp_path))
        assert home.requests[0]["projectName"] is None

    @pytest.mark.xfail(
        strict=True,
        reason="B128: buildings_raster_toSTL computes dxdy, inputCRS, "
               "outputCRS and dataSourceName from the arguments and then "
               "calls getBuildingsFromRectangle with the bounding box and "
               "withElevation only, so --dataSourceName and --inputCRS are "
               "silently ignored and the default datasource in WGS84 is "
               "always used. See the consolidated findings issue.",
    )
    def test_an_explicit_datasource_and_inputcrs_should_reach_the_toolkit(self, home, tmp_path):
        CLI.buildings_raster_toSTL(
            self._arguments(tmp_path, dataSourceName="BLDG_2020", inputCRS=2039)
        )
        _, kwargs = home.toolkit.call("getBuildingsFromRectangle")
        assert kwargs["dataSourceName"] == "BLDG_2020"
        assert kwargs["inputCRS"] == 2039

    def test_it_currently_forwards_only_the_bounding_box(self, home, tmp_path):
        """Characterisation of B128."""
        CLI.buildings_raster_toSTL(
            self._arguments(tmp_path, dataSourceName="BLDG_2020", inputCRS=2039, dxdy=7)
        )
        _, kwargs = home.toolkit.call("getBuildingsFromRectangle")
        assert set(kwargs) == {"minx", "miny", "maxx", "maxy", "withElevation"}


# ---------------------------------------------------------------------------
# landcover
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetLandcover:
    def test_the_roughness_flag_selects_getroughness(self, home, tmp_path):
        CLI.get_landocver(_landcoverArguments(tmp_path, roughness=True))
        assert [entry[0] for entry in home.toolkit.calls] == ["getRoughness"]

    def test_without_the_roughness_flag_getlandcover_is_used(self, home, tmp_path):
        CLI.get_landocver(_landcoverArguments(tmp_path, roughness=False))
        assert [entry[0] for entry in home.toolkit.calls] == ["getLandCover"]

    def test_the_bounding_box_is_passed_positionally(self, home, tmp_path):
        CLI.get_landocver(_landcoverArguments(tmp_path))
        _, args, _ = home.toolkit.call("getRoughness")
        assert args == ("1", "2", "3", "4")

    def test_the_resolution_and_crs_strings_are_converted_to_numbers(self, home, tmp_path):
        CLI.get_landocver(_landcoverArguments(tmp_path, dxdy="60", inputCRS="2039",
                                              resolution="12.5"))
        _, _, kwargs = home.toolkit.call("getRoughness")
        assert kwargs["dxdy"] == 60 and isinstance(kwargs["dxdy"], int)
        assert kwargs["inputCRS"] == 2039 and isinstance(kwargs["inputCRS"], int)
        assert kwargs["resolution"] == pytest.approx(12.5)

    def test_a_missing_inputcrs_falls_back_to_wsg84(self, home, tmp_path):
        CLI.get_landocver(_landcoverArguments(tmp_path, inputCRS=None))
        assert home.toolkit.call("getRoughness")[2]["inputCRS"] == 4326

    def test_the_building_flag_and_datasource_reach_the_toolkit(self, home, tmp_path):
        CLI.get_landocver(_landcoverArguments(tmp_path, isBuilding=True,
                                              dataSourceName="LANDCOVER_2020"))
        kwargs = home.toolkit.call("getRoughness")[2]
        assert kwargs["isBuilding"] is True
        assert kwargs["dataSourceName"] == "LANDCOVER_2020"

    def test_getlandcover_is_called_without_the_roughness_only_arguments(self, home, tmp_path):
        CLI.get_landocver(_landcoverArguments(tmp_path, roughness=False))
        kwargs = home.toolkit.call("getLandCover")[2]
        assert set(kwargs) == {"dxdy", "inputCRS", "dataSourceName"}

    def test_it_writes_a_csv_of_the_flattened_landcover_grid(self, home, tmp_path):
        target = tmp_path / "out.csv"
        CLI.get_landocver(_landcoverArguments(tmp_path, filePath=str(target)))

        frame = pandas.read_csv(target)
        assert list(frame.columns) == ["lat", "lon", "landcover"]
        assert list(frame["lat"]) == [0.0, 1.0]
        assert list(frame["landcover"]) == [10, 20]

    def test_a_z0_coordinate_adds_a_z0_column(self, home, tmp_path):
        target = tmp_path / "withZ0.csv"
        home.landcover = _FakeXarray(
            fields=dict(lat=[[0.0, 1.0]], lon=[[2.0, 3.0]],
                        landcover=[[10, 20]], z0=[[0.1, 0.4]]),
            coords=("z0",),
        )
        CLI.get_landocver(_landcoverArguments(tmp_path, filePath=str(target)))

        frame = pandas.read_csv(target)
        assert list(frame.columns) == ["lat", "lon", "landcover", "z0"]
        assert list(frame["z0"]) == [pytest.approx(0.1), pytest.approx(0.4)]

    def test_a_missing_projectname_is_defaulted_onto_the_arguments(self, home, tmp_path):
        arguments = _landcoverArguments(tmp_path)
        del arguments.projectName
        CLI.get_landocver(arguments)
        assert arguments.projectName is None
        assert home.requests[0]["projectName"] is None

    def test_a_missing_datasourcename_is_defaulted_onto_the_arguments(self, home, tmp_path):
        arguments = _landcoverArguments(tmp_path)
        del arguments.dataSourceName
        CLI.get_landocver(arguments)
        assert arguments.dataSourceName is None

    def test_it_asks_for_the_landcover_toolkit(self, home, tmp_path):
        CLI.get_landocver(_landcoverArguments(tmp_path))
        assert home.requests == [dict(toolkitName=_FakeToolkitHome.GIS_LANDCOVER,
                                      projectName="LANDCOVER_PROJECT")]

    @pytest.mark.xfail(
        strict=True,
        reason="B127: get_landocver reads the wind direction from the "
               "misspelled attribute `arguments.windDirectionis` instead of "
               "`arguments.windDirection`, so passing --windDirection -- "
               "mandatory for the isBuilding=True roughness case -- raises "
               "AttributeError before the toolkit is ever called. See the "
               "consolidated findings issue.",
    )
    def test_a_wind_direction_should_reach_the_roughness_calculation(self, home, tmp_path):
        CLI.get_landocver(_landcoverArguments(tmp_path, isBuilding=True,
                                              windDirection="270", resolution="10"))
        assert home.toolkit.call("getRoughness")[2]["windMeteorologicalDirection"] == \
            pytest.approx(270.0)

    def test_supplying_a_wind_direction_currently_raises_attributeerror(self, home, tmp_path):
        """Characterisation of B127."""
        with pytest.raises(AttributeError, match="windDirectionis"):
            CLI.get_landocver(_landcoverArguments(tmp_path, windDirection="270"))
        assert home.toolkit.calls == []

    def test_omitting_the_wind_direction_keeps_the_command_working(self, home, tmp_path):
        CLI.get_landocver(_landcoverArguments(tmp_path, windDirection=None))
        assert home.toolkit.call("getRoughness")[2]["windMeteorologicalDirection"] is None
