"""simulations/windProfile/toolkit.py: WindProfileToolkit and its
presentation layer.

Three defects surfaced while probing:

* B110: ``getWindProfile`` recomputes the friction velocity from the very
  height it is about to evaluate --
  ``U_star = ws * KARMAN / log(z / z0)`` then
  ``U_z = (U_star / KARMAN) * log(z / z0)`` -- so the two logarithms cancel
  algebraically and ``U_z == ws`` for every height. The returned "wind
  profile" is perfectly flat: it has no shear, and the roughness length
  ``z0`` has no effect on it whatsoever (a 50x change in z0 yields an
  identical profile). ``U_star`` has to come from a fixed reference height
  (the station's measurement height) for the log law to extrapolate
  anything.
* B111: the height loop is ``numpy.arange(0, height + dz, dz)``, so it
  always evaluates the log law at z=0, where it is singular. The
  non-urban branch emits ``NaN`` for u/v/U_z in that row; the urban
  branch emits ``-0.0`` for the same input. Every returned profile
  therefore starts with a junk row, and the two branches disagree about
  what that row contains.
* B112: ``_getStationsInRegion`` opens ``'wind_stations.json'`` as a bare
  relative path, so it resolves against the process working directory and
  raises ``FileNotFoundError`` anywhere but the module's own directory.
  Its sibling ``getSpatialWind`` reads the same file correctly, via
  ``os.path.dirname(os.path.abspath(__file__))``.

``getSpatialWind`` itself is not exercised: it needs a populated
LandCover toolkit plus a live IMS API token, which is integration-shaped.
"""
import numpy
import pytest
import xarray

from hera.measurements.GIS.utils import KARMAN


@pytest.fixture()
def tk(unit_toolkit_factory):
    from hera import toolkitHome

    return unit_toolkit_factory(toolkitHome.WINDPROFILE)


def _grid(ws=5.0, wd=270.0, z0=0.1, hc=None, ll=5.0, n=1):
    """A minimal landcover-shaped DataArray: 2D lat/lon coords plus the
    wind and roughness fields getWindProfile reads."""
    coords = {
        "lat": (("i", "j"), numpy.full((n, n), 32.0)),
        "lon": (("i", "j"), numpy.full((n, n), 34.0)),
        "ws": (("i", "j"), numpy.full((n, n), ws)),
        "wd": (("i", "j"), numpy.full((n, n), wd)),
        "z0": (("i", "j"), numpy.full((n, n), z0)),
    }
    if hc is not None:
        coords["hc"] = (("i", "j"), numpy.full((n, n), hc))
        coords["ll"] = (("i", "j"), numpy.full((n, n), ll))
    return xarray.DataArray(numpy.zeros((n, n)), dims=("i", "j"), coords=coords)


@pytest.mark.unit
class TestConstruction:
    def test_it_builds_with_a_presentation_layer(self, tk):
        assert tk.presentation is not None

    def test_the_presentation_layer_points_back_at_the_toolkit(self, tk):
        assert tk.presentation.datalayer is tk


@pytest.mark.unit
class TestFindLatLonIndex:
    def test_it_returns_the_index_of_the_nearest_grid_cell(self, tk):
        grid = xarray.DataArray(
            numpy.zeros((2, 2)), dims=("i", "j"),
            coords={
                "lat": (("i", "j"), numpy.array([[32.0, 32.0], [32.1, 32.1]])),
                "lon": (("i", "j"), numpy.array([[34.0, 34.1], [34.0, 34.1]])),
            },
        )
        assert tk._find_lat_lon_index_in_xarray(32.1, 34.1, grid) == (1, 1)
        assert tk._find_lat_lon_index_in_xarray(32.0, 34.0, grid) == (0, 0)

    def test_an_off_grid_point_snaps_to_the_closest_cell(self, tk):
        grid = xarray.DataArray(
            numpy.zeros((2, 2)), dims=("i", "j"),
            coords={
                "lat": (("i", "j"), numpy.array([[32.0, 32.0], [32.1, 32.1]])),
                "lon": (("i", "j"), numpy.array([[34.0, 34.1], [34.0, 34.1]])),
            },
        )
        assert tk._find_lat_lon_index_in_xarray(32.099, 34.098, grid) == (1, 1)


@pytest.mark.unit
class TestGetWindProfileShape:
    def test_it_returns_one_row_per_cell_and_height(self, tk):
        df = tk.getWindProfile(_grid(n=2), height=20.0, dz=10.0)
        # 4 cells x heights [0, 10, 20]
        assert len(df) == 4 * 3

    def test_it_reports_the_expected_columns(self, tk):
        df = tk.getWindProfile(_grid(), height=10.0, dz=10.0)
        assert set(df.columns) == {"lat", "lon", "height", "u", "v", "U_z", "direction"}

    def test_the_wind_direction_is_carried_through_unchanged(self, tk):
        df = tk.getWindProfile(_grid(wd=123.0), height=10.0, dz=10.0)
        assert (df["direction"] == 123.0).all()

    def test_the_components_follow_the_meteorological_convention(self, tk):
        """wd=0 means wind FROM the north, so it blows toward the south:
        v must be negative and u must vanish."""
        df = tk.getWindProfile(_grid(wd=0.0), height=10.0, dz=10.0)
        row = df[df["height"] == 10.0].iloc[0]
        assert row["u"] == pytest.approx(0.0, abs=1e-9)
        assert row["v"] == pytest.approx(-row["U_z"])

    def test_the_components_are_consistent_with_the_speed(self, tk):
        df = tk.getWindProfile(_grid(wd=225.0), height=10.0, dz=10.0)
        row = df[df["height"] == 10.0].iloc[0]
        assert numpy.hypot(row["u"], row["v"]) == pytest.approx(abs(row["U_z"]))


@pytest.mark.unit
class TestGetWindProfileHasNoShear:
    """B110: see the module docstring."""

    @pytest.mark.xfail(
        strict=True,
        reason="B110: U_star is recomputed from the loop's own z, so the "
               "two logarithms cancel and U_z == ws at every height -- the "
               "profile has no shear at all. See the consolidated findings "
               "issue.",
    )
    def test_the_wind_should_increase_with_height(self, tk):
        df = tk.getWindProfile(_grid(ws=5.0, z0=0.1), height=100.0, dz=25.0)
        above_ground = df[df["height"] > 0].sort_values("height")
        assert above_ground["U_z"].iloc[-1] > above_ground["U_z"].iloc[0]

    @pytest.mark.xfail(
        strict=True,
        reason="B110: z0 cancels out of the profile entirely, so a rougher "
               "surface produces exactly the same wind speeds. See the "
               "consolidated findings issue.",
    )
    def test_a_rougher_surface_should_slow_the_wind_near_the_ground(self, tk):
        smooth = tk.getWindProfile(_grid(z0=0.03), height=50.0, dz=25.0)
        rough = tk.getWindProfile(_grid(z0=1.5), height=50.0, dz=25.0)
        at_25 = lambda df: df[df["height"] == 25.0]["U_z"].iloc[0]
        assert at_25(rough) < at_25(smooth)

    def test_the_profile_is_currently_flat_at_the_input_speed(self, tk):
        """Characterisation of B110: U_z is exactly the station wind speed
        at every height above ground."""
        df = tk.getWindProfile(_grid(ws=5.0), height=100.0, dz=25.0)
        above_ground = df[df["height"] > 0]
        assert above_ground["U_z"].values == pytest.approx([5.0] * len(above_ground))

    def test_roughness_currently_makes_no_difference_at_all(self, tk):
        """Characterisation of B110: a 50x change in z0 leaves the profile
        numerically unchanged."""
        smooth = tk.getWindProfile(_grid(z0=0.03), height=50.0, dz=25.0)
        rough = tk.getWindProfile(_grid(z0=1.5), height=50.0, dz=25.0)
        smooth_above = smooth[smooth["height"] > 0]["U_z"].values
        rough_above = rough[rough["height"] > 0]["U_z"].values
        assert rough_above == pytest.approx(smooth_above)


@pytest.mark.unit
class TestGetWindProfileGroundRow:
    """B111: see the module docstring."""

    @pytest.mark.xfail(
        strict=True,
        reason="B111: the height loop starts at z=0, where the log law is "
               "singular, so the first row of every non-urban profile is "
               "NaN instead of a usable value. See the consolidated "
               "findings issue.",
    )
    def test_the_ground_row_should_be_a_real_number(self, tk):
        df = tk.getWindProfile(_grid(), height=10.0, dz=10.0)
        ground = df[df["height"] == 0.0]
        assert not ground["U_z"].isna().any()

    def test_the_ground_row_is_currently_nan_without_a_canopy(self, tk):
        """Characterisation of B111."""
        df = tk.getWindProfile(_grid(), height=10.0, dz=10.0)
        ground = df[df["height"] == 0.0]
        assert ground["U_z"].isna().all()
        assert ground["u"].isna().all()
        assert ground["v"].isna().all()

    def test_the_ground_row_is_currently_zero_with_a_canopy(self, tk):
        """Characterisation of B111: the urban branch produces -0.0 where
        the non-urban branch produces NaN, for the same singular height."""
        df = tk.getWindProfile(_grid(hc=10.0), height=20.0, dz=10.0)
        ground = df[df["height"] == 0.0]
        assert not ground["U_z"].isna().any()
        assert ground["U_z"].values == pytest.approx([0.0] * len(ground))


@pytest.mark.unit
class TestGetWindProfileCanopyBranch:
    def test_a_canopy_below_two_metres_is_skipped_entirely(self, tk):
        """The urban branch `continue`s when hc <= 2, so a short canopy
        produces no rows at all rather than falling back to the
        non-urban log law."""
        df = tk.getWindProfile(_grid(hc=1.0), height=20.0, dz=10.0)
        assert len(df) == 0

    def test_a_tall_canopy_produces_a_row_for_every_height(self, tk):
        df = tk.getWindProfile(_grid(hc=10.0), height=20.0, dz=10.0)
        assert sorted(df["height"].tolist()) == [0.0, 10.0, 20.0]


@pytest.mark.unit
class TestInterpolatedWindFields:
    @pytest.fixture()
    def grid(self):
        return xarray.Dataset(
            {"elevation": (("i", "j"), numpy.zeros((2, 2)))},
            coords={
                "lon": (("i", "j"), numpy.array([[34.0, 34.1], [34.0, 34.1]])),
                "lat": (("i", "j"), numpy.array([[32.0, 32.0], [32.1, 32.1]])),
            },
        )

    STATIONS = [[32.0, 34.0, 10.0, [5.0, 270.0]], [32.1, 34.1, 10.0, [7.0, 180.0]]]

    def test_it_adds_ws_and_wd_fields_to_the_grid(self, tk, grid):
        result = tk.add_interpolated_ws_wd(grid, self.STATIONS)
        assert "ws" in result and "wd" in result

    def test_the_station_values_are_reproduced_at_the_station_locations(self, tk, grid):
        result = tk.add_interpolated_ws_wd(grid, self.STATIONS)
        ws = numpy.asarray(result["ws"])
        assert ws[0, 0] == pytest.approx(5.0)
        assert ws[1, 1] == pytest.approx(7.0)

    def test_intermediate_cells_are_interpolated_between_the_stations(self, tk, grid):
        result = tk.add_interpolated_ws_wd(grid, self.STATIONS)
        ws = numpy.asarray(result["ws"])
        assert 5.0 < ws[0, 1] < 7.0

    def test_the_pointwise_helper_returns_a_speed_direction_pair(self, tk):
        ws, wd = tk._interpolate_wd_ws(34.0, 32.0, 0.0, self.STATIONS)
        assert ws == pytest.approx(5.0)
        assert wd == pytest.approx(270.0)


@pytest.mark.unit
class TestGetStationsInRegionCannotFindItsDataFile:
    """B112: see the module docstring."""

    @pytest.mark.xfail(
        strict=True,
        reason="B112: _getStationsInRegion opens 'wind_stations.json' as a "
               "bare relative path, so it resolves against the process cwd "
               "instead of the module directory. See the consolidated "
               "findings issue.",
    )
    def test_it_should_read_the_station_list_shipped_beside_the_module(self, tk, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        assert tk._getStationsInRegion(34.0, 32.0, 35.0, 33.0, 4326) is not None

    def test_it_currently_raises_filenotfounderror_from_any_other_directory(self, tk, tmp_path, monkeypatch):
        """Characterisation of B112."""
        monkeypatch.chdir(tmp_path)
        with pytest.raises(FileNotFoundError, match="wind_stations.json"):
            tk._getStationsInRegion(34.0, 32.0, 35.0, 33.0, 4326)

    def test_the_data_file_it_wants_really_does_ship_with_the_module(self):
        """So B112 is purely a path bug, not a missing resource."""
        import os

        import hera.simulations.windProfile.toolkit as module

        shipped = os.path.join(os.path.dirname(os.path.abspath(module.__file__)), "wind_stations.json")
        assert os.path.isfile(shipped)


class _FakePlot:
    """Stands in for the matplotlib image handle plotWindDirections wants."""

    def get_array(self):
        return numpy.zeros((4, 4))

    def get_extent(self):
        return (0.0, 1.0, 0.0, 1.0)


@pytest.mark.unit
class TestPlotWindDirections:
    @pytest.fixture()
    def landcover(self):
        return xarray.Dataset(coords={
            "lat": (("i", "j"), numpy.array([[32.0, 32.0], [32.1, 32.1]])),
            "lon": (("i", "j"), numpy.array([[34.0, 34.1], [34.0, 34.1]])),
            "wd": (("i", "j"), numpy.full((2, 2), 270.0)),
        })

    def test_it_returns_a_figure_when_not_showing(self, tk, landcover):
        fig = tk.presentation.plotWindDirections(_FakePlot(), landcover, show_photo=False)
        assert fig is not None

    def test_it_draws_one_arrow_per_grid_cell(self, tk, landcover):
        fig = tk.presentation.plotWindDirections(_FakePlot(), landcover, show_photo=False)
        ax = fig.axes[0]
        assert len(ax.patches) == 4

    def test_the_axes_are_clamped_to_the_underlying_image_extent(self, tk, landcover):
        fig = tk.presentation.plotWindDirections(_FakePlot(), landcover, show_photo=False)
        ax = fig.axes[0]
        assert ax.get_xlim() == pytest.approx((0.0, 1.0))
        assert ax.get_ylim() == pytest.approx((0.0, 1.0))
