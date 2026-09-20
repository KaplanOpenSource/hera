"""simulations/windProfile/toolkit.py: the two network-facing members.

``test_windprofile_toolkit.py`` covers ``getWindProfile`` (B110/B111),
``_find_lat_lon_index_in_xarray``, the interpolated wind fields,
``_getStationsInRegion`` (B112) and the presentation layer, and records
``getSpatialWind`` as "integration-shaped, not exercised".  It is in fact
drivable hermetically -- the landcover call and the IMS call are both
single method calls that can be replaced on their classes -- so this file
covers the pair:

* ``getSpatialWind``          -- landcover + roughness + station list
* ``_getWindSpeedDirection``  -- the IMS Envista polling loop

``requests.request`` is monkeypatched for every test that touches the
second one; the unit layer's autouse socket guard would refuse a real
call anyway, and the method retries 20 times per station before giving
up, so an unpatched call would also be slow.

Three defects surfaced:

* B250: ``getSpatialWind`` calls
  ``landcover_tk.getRoughnessFromLandcover(xarray, dxdy)``, but that
  method's signature is ``(self, landcover,
  windMeteorologicalDirection=None, resolution=None, ...)``.  The grid
  resolution is therefore passed as a *wind direction*: with the default
  ``dxdy=30`` the roughness is computed for a wind from 30 degrees, and
  the ``resolution`` argument the caller meant to set is left at None.
* B251: ``_getWindSpeedDirection`` ends with
  ``max(datetime_objects)`` over a list that is only appended to for
  stations that returned usable data.  When none does -- every station
  offline, a bad token, or every reading out of range -- the call raises
  ``ValueError: max() iterable argument is empty`` instead of returning
  an empty list.
* B252: ``ws`` and ``wd`` are bound only inside
  ``for channel in data['data'][0]['channels']``, and are then read
  unconditionally.  A station that reports data but no WS/WD channel --
  a rain-only station, which the shipped ``wind_stations.json`` list does
  contain -- raises ``UnboundLocalError``.  Worse, when only one of the
  two is missing the loop silently reuses the *previous station's* value.
* B258: the retry loop assigns ``data`` before the subscript that
  validates it, so an error payload -- from a rejected token, say --
  survives all 20 attempts, passes the ``if data:`` guard that exists to
  skip unreachable stations, and is then dereferenced outside the try,
  raising ``KeyError: 'data'`` instead of skipping the station.
"""
import json

import pytest


@pytest.fixture()
def tk(unit_toolkit_factory):
    from hera import toolkitHome

    return unit_toolkit_factory(toolkitHome.WINDPROFILE)


TOKEN = {"Authorization": "ApiToken 00000000-0000-0000-0000-000000000000"}
FRESH = "2024-01-01T12:00:00+02:00"
STALE = "2024-01-01T10:00:00+02:00"


def _station(stationId=2, lat=32.0, lon=35.0, height=375):
    """The shape ``wind_stations.json`` really has: the station id is
    attributes[0] and the location attributes[2]."""
    return {
        "name": f"STATION {stationId}",
        "description": "",
        "height_above_sea": height,
        "attributes": [
            {"name": "stationId", "value": stationId},
            {"name": "stationsTag", "value": "(null)"},
            {"name": "location", "value": {"latitude": lat, "longitude": lon}},
        ],
    }


class _Response:
    """Only ``.text`` is read, and only through json.loads."""

    def __init__(self, payload):
        self.text = json.dumps(payload)


def _reading(ws=5.0, wd=270.0, datetime_str=FRESH, channels=None):
    if channels is None:
        channels = [{"name": "WS", "value": ws}, {"name": "WD", "value": wd}]
    return {"data": [{"datetime": datetime_str, "channels": channels}]}


def _fixedResponse(monkeypatch, payload):
    import requests

    monkeypatch.setattr(requests, "request", lambda *args, **kwargs: _Response(payload))


def _responseSequence(monkeypatch, payloads):
    """One payload per station, in order."""
    import requests

    remaining = list(payloads)

    def _request(*args, **kwargs):
        return _Response(remaining.pop(0))

    monkeypatch.setattr(requests, "request", _request)


@pytest.mark.unit
class TestGetWindSpeedDirection:
    def test_a_station_with_data_is_returned_as_lat_lon_height_reading(
        self, tk, monkeypatch
    ):
        _fixedResponse(monkeypatch, _reading(ws=5.0, wd=270.0))
        assert tk._getWindSpeedDirection([_station()], TOKEN) == [
            [32.0, 35.0, 375, [5.0, 270.0]]
        ]

    def test_the_reading_order_is_speed_then_direction(self, tk, monkeypatch):
        """The channels arrive in either order, so the pairing must come from
        the channel names rather than from their position."""
        _fixedResponse(
            monkeypatch,
            _reading(channels=[{"name": "WD", "value": 123.0},
                               {"name": "WS", "value": 4.0}]),
        )
        assert tk._getWindSpeedDirection([_station()], TOKEN)[0][3] == [4.0, 123.0]

    def test_unrelated_channels_are_ignored(self, tk, monkeypatch):
        """An IMS station reports rain, temperature and gusts alongside the
        mean wind."""
        _fixedResponse(
            monkeypatch,
            _reading(channels=[{"name": "Rain", "value": 0.2},
                               {"name": "WSmax", "value": 19.0},
                               {"name": "WS", "value": 4.0},
                               {"name": "TD", "value": 21.0},
                               {"name": "WD", "value": 123.0}]),
        )
        assert tk._getWindSpeedDirection([_station()], TOKEN)[0][3] == [4.0, 123.0]

    def test_the_authorization_header_is_forwarded(self, tk, monkeypatch):
        import requests

        seen = {}

        def _request(method, url, headers=None, **kwargs):
            seen["method"] = method
            seen["url"] = url
            seen["headers"] = headers
            return _Response(_reading())

        monkeypatch.setattr(requests, "request", _request)
        tk._getWindSpeedDirection([_station(stationId=17)], TOKEN)
        assert seen["method"] == "GET"
        assert seen["headers"] == {"Authorization": TOKEN["Authorization"]}
        assert seen["url"].endswith("/stations/17/data/latest")

    def test_every_station_is_polled(self, tk, monkeypatch):
        _responseSequence(monkeypatch, [_reading(ws=1.0), _reading(ws=2.0),
                                        _reading(ws=3.0)])
        result = tk._getWindSpeedDirection(
            [_station(1), _station(2), _station(3)], TOKEN
        )
        assert [row[3][0] for row in result] == [1.0, 2.0, 3.0]

    def test_an_implausible_wind_speed_is_dropped(self, tk, monkeypatch):
        """The filter is |ws| < 20 m/s and |wd| <= 360; 25 m/s at a surface
        station is a sensor fault, not weather."""
        _responseSequence(monkeypatch, [_reading(ws=25.0), _reading(ws=5.0)])
        result = tk._getWindSpeedDirection([_station(1), _station(2)], TOKEN)
        assert [row[3][0] for row in result] == [5.0]

    def test_an_out_of_range_direction_is_dropped(self, tk, monkeypatch):
        _responseSequence(monkeypatch, [_reading(wd=999.0), _reading(wd=90.0)])
        result = tk._getWindSpeedDirection([_station(1), _station(2)], TOKEN)
        assert [row[3][1] for row in result] == [90.0]

    def test_a_station_that_never_answers_is_skipped(self, tk, monkeypatch):
        """The retry loop swallows every exception; after 20 attempts data is
        still None and the station is left out."""
        import requests

        calls = {"n": 0}

        def _request(*args, **kwargs):
            calls["n"] += 1
            if calls["n"] <= 20:
                raise RuntimeError("connection reset")
            return _Response(_reading(ws=6.0))

        monkeypatch.setattr(requests, "request", _request)
        result = tk._getWindSpeedDirection([_station(1), _station(2)], TOKEN)
        assert calls["n"] == 21
        assert [row[3][0] for row in result] == [6.0]

    def test_readings_older_than_fifteen_minutes_are_dropped(self, tk, monkeypatch):
        """A wind field must be synoptic: only stations within 15 minutes of
        the freshest reading are kept."""
        _responseSequence(
            monkeypatch,
            [_reading(ws=5.0, datetime_str=FRESH),
             _reading(ws=6.0, datetime_str=STALE)],
        )
        result = tk._getWindSpeedDirection([_station(1), _station(2)], TOKEN)
        assert [row[3][0] for row in result] == [5.0]

    def test_a_reading_just_inside_the_window_is_kept(self, tk, monkeypatch):
        _responseSequence(
            monkeypatch,
            [_reading(ws=5.0, datetime_str="2024-01-01T12:00:00+02:00"),
             _reading(ws=6.0, datetime_str="2024-01-01T11:46:00+02:00")],
        )
        result = tk._getWindSpeedDirection([_station(1), _station(2)], TOKEN)
        assert sorted(row[3][0] for row in result) == [5.0, 6.0]

    def test_the_freshest_reading_sets_the_reference_time(self, tk, monkeypatch):
        """It is max(datetimes), not the first station's, so order does not
        matter."""
        _responseSequence(
            monkeypatch,
            [_reading(ws=6.0, datetime_str=STALE),
             _reading(ws=5.0, datetime_str=FRESH)],
        )
        result = tk._getWindSpeedDirection([_station(1), _station(2)], TOKEN)
        assert [row[3][0] for row in result] == [5.0]


@pytest.mark.unit
class TestGetWindSpeedDirectionWithNoUsableStations:
    """B251: see the module docstring."""

    @pytest.mark.xfail(
        strict=True,
        reason="B251: _getWindSpeedDirection ends with max(datetime_objects) "
               "over a list appended to only for stations that returned usable "
               "data. If none does -- everything offline, a rejected token, or "
               "every reading out of range -- it raises ValueError('max() "
               "iterable argument is empty') rather than returning the empty "
               "list its filtering implies. See the consolidated findings issue.",
    )
    def test_an_empty_station_list_returns_no_stations(self, tk):
        assert tk._getWindSpeedDirection([], TOKEN) == []

    def test_an_empty_station_list_currently_raises(self, tk):
        """Characterisation of B251 in its simplest form: no HTTP needed."""
        with pytest.raises(ValueError, match="max"):
            tk._getWindSpeedDirection([], TOKEN)

    @pytest.mark.xfail(
        strict=True,
        reason="B251: the same max() over an empty list when every station "
               "answers but every reading is out of range. "
               "See the consolidated findings issue.",
    )
    def test_all_readings_rejected_returns_no_stations(self, tk, monkeypatch):
        _fixedResponse(monkeypatch, _reading(ws=25.0))
        assert tk._getWindSpeedDirection([_station()], TOKEN) == []

    def test_all_readings_rejected_currently_raises(self, tk, monkeypatch):
        """Characterisation of B251: two stations answer, both readings are
        out of range, and the empty list reaches max()."""
        _fixedResponse(monkeypatch, _reading(ws=25.0))
        with pytest.raises(ValueError, match="max\\(\\) iterable argument is empty"):
            tk._getWindSpeedDirection([_station(1), _station(2)], TOKEN)


@pytest.mark.unit
class TestGetWindSpeedDirectionWithAnErrorPayload:
    """B258: an error response from the API is not skipped.

    The retry loop assigns ``data = json.loads(...)`` and only then reads
    ``data['data'][0]['datetime']``.  For an error payload -- ``{"errorMessage":
    "unauthorized"}`` from a rejected token -- the parse succeeds and the
    subscript fails, so after 20 rejected attempts ``data`` is still the error
    dict.  It is truthy, so the ``if data:`` guard that exists to skip
    unreachable stations lets it through, and the very next line dereferences
    ``data['data']`` outside the try, raising ``KeyError: 'data'``.
    """

    UNAUTHORIZED = {"errorMessage": "unauthorized"}

    @pytest.mark.xfail(
        strict=True,
        reason="B258: the retry loop assigns `data` before the subscript that "
               "validates it, so a response it rejected 20 times still passes the "
               "`if data:` guard; the error payload is then dereferenced outside "
               "the try/except and raises KeyError: 'data'. A rejected IMS token "
               "crashes the call instead of yielding no stations. "
               "See the consolidated findings issue.",
    )
    def test_an_error_payload_is_treated_as_no_data(self, tk, monkeypatch):
        _fixedResponse(monkeypatch, self.UNAUTHORIZED)
        assert tk._getWindSpeedDirection([_station()], TOKEN) == []

    def test_an_error_payload_currently_raises_keyerror(self, tk, monkeypatch):
        """Characterisation of B258."""
        _fixedResponse(monkeypatch, self.UNAUTHORIZED)
        with pytest.raises(KeyError, match="data"):
            tk._getWindSpeedDirection([_station()], TOKEN)

    def test_an_unreachable_station_is_skipped_by_contrast(self, tk, monkeypatch):
        """Characterisation of B258's boundary: when the request itself
        raises, `data` stays None and the guard does work."""
        import requests

        def _request(*args, **kwargs):
            raise RuntimeError("connection reset")

        monkeypatch.setattr(requests, "request", _request)
        with pytest.raises(ValueError, match="max"):
            tk._getWindSpeedDirection([_station()], TOKEN)


@pytest.mark.unit
class TestGetWindSpeedDirectionWithMissingChannels:
    """B252: see the module docstring."""

    RAIN_ONLY = _reading(channels=[{"name": "Rain", "value": 0.0}])

    @pytest.mark.xfail(
        strict=True,
        reason="B252: ws and wd are bound only inside the channel loop and "
               "then read unconditionally, so a station that answers with no "
               "WS/WD channel raises UnboundLocalError. The shipped "
               "wind_stations.json does list stations whose monitor set is not "
               "guaranteed to include both. See the consolidated findings issue.",
    )
    def test_a_station_without_wind_channels_is_skipped(self, tk, monkeypatch):
        _fixedResponse(monkeypatch, self.RAIN_ONLY)
        assert tk._getWindSpeedDirection([_station()], TOKEN) == []

    def test_a_station_without_wind_channels_currently_raises(self, tk, monkeypatch):
        """Characterisation of B252."""
        _fixedResponse(monkeypatch, self.RAIN_ONLY)
        with pytest.raises(UnboundLocalError, match="ws"):
            tk._getWindSpeedDirection([_station()], TOKEN)

    @pytest.mark.xfail(
        strict=True,
        reason="B252, the silent half: because ws/wd survive between loop "
               "iterations, a station missing only one of the two channels is "
               "recorded with the PREVIOUS station's value for it instead of "
               "being skipped. See the consolidated findings issue.",
    )
    def test_a_station_missing_one_channel_does_not_inherit_the_previous_one(
        self, tk, monkeypatch
    ):
        _responseSequence(
            monkeypatch,
            [_reading(ws=5.0, wd=270.0),
             _reading(channels=[{"name": "WS", "value": 7.0}])],
        )
        result = tk._getWindSpeedDirection([_station(1), _station(2)], TOKEN)
        assert [row[3] for row in result] == [[5.0, 270.0]]

    def test_a_station_missing_one_channel_currently_inherits_it(
        self, tk, monkeypatch
    ):
        """Characterisation of B252's silent half: station 2 reports no
        direction and is recorded with station 1's 270 degrees."""
        _responseSequence(
            monkeypatch,
            [_reading(ws=5.0, wd=270.0),
             _reading(channels=[{"name": "WS", "value": 7.0}])],
        )
        result = tk._getWindSpeedDirection([_station(1), _station(2)], TOKEN)
        assert [row[3] for row in result] == [[5.0, 270.0], [7.0, 270.0]]


@pytest.mark.unit
class TestGetSpatialWind:
    """The composition step: landcover -> roughness -> station readings."""

    @pytest.fixture()
    def wired(self, monkeypatch):
        """Replace the two outward-facing calls on their CLASSES, and record
        what they were handed."""
        from hera.measurements.GIS.raster.landcover import LandCoverToolkit
        from hera.simulations.windProfile.toolkit import WindProfileToolkit

        seen = {}

        def _getLandCover(self, *args, **kwargs):
            seen["landCover"] = (args, kwargs)
            return "landcover-grid"

        def _getRoughness(self, *args, **kwargs):
            seen["roughness"] = (args, kwargs)
            return "roughness-grid"

        def _getWind(self, stations, token):
            seen["stations"] = stations
            seen["token"] = token
            return [[32.0, 35.0, 375, [5.0, 270.0]]]

        monkeypatch.setattr(LandCoverToolkit, "getLandCover", _getLandCover)
        monkeypatch.setattr(
            LandCoverToolkit, "getRoughnessFromLandcover", _getRoughness
        )
        monkeypatch.setattr(
            WindProfileToolkit, "_getWindSpeedDirection", _getWind
        )
        return seen

    def test_it_returns_the_roughness_grid_and_the_station_readings(self, tk, wired):
        grid, stations = tk.getSpatialWind(32.0, 35.0, 33.0, 36.0, TOKEN)
        assert grid == "roughness-grid"
        assert stations == [[32.0, 35.0, 375, [5.0, 270.0]]]

    def test_the_bounding_box_and_resolution_reach_the_landcover_toolkit(
        self, tk, wired
    ):
        tk.getSpatialWind(32.0, 35.0, 33.0, 36.0, TOKEN, dxdy=50)
        args, kwargs = wired["landCover"]
        assert args == (32.0, 35.0, 33.0, 36.0, 50)
        assert kwargs["inputCRS"] == 4326

    def test_the_landcover_datasource_name_is_forwarded(self, tk, wired):
        tk.getSpatialWind(
            32.0, 35.0, 33.0, 36.0, TOKEN, landcover_DataSource="LC_2020"
        )
        assert wired["landCover"][1]["dataSourceName"] == "LC_2020"

    def test_the_roughness_is_derived_from_the_landcover_grid(self, tk, wired):
        tk.getSpatialWind(32.0, 35.0, 33.0, 36.0, TOKEN)
        assert wired["roughness"][0][0] == "landcover-grid"

    def test_the_token_is_handed_to_the_station_poller(self, tk, wired):
        tk.getSpatialWind(32.0, 35.0, 33.0, 36.0, TOKEN)
        assert wired["token"] is TOKEN

    def test_the_shipped_station_list_is_read_from_the_module_directory(
        self, tk, wired, tmp_path, monkeypatch
    ):
        """Unlike _getStationsInRegion (B112), this one builds the path from
        __file__, so it works from any working directory."""
        monkeypatch.chdir(tmp_path)
        tk.getSpatialWind(32.0, 35.0, 33.0, 36.0, TOKEN)
        assert len(wired["stations"]) == 83

    def test_every_shipped_station_has_the_fields_the_poller_reads(self, tk, wired):
        """attributes[0] as the id, attributes[2] as the location, and a
        top-level height_above_sea."""
        tk.getSpatialWind(32.0, 35.0, 33.0, 36.0, TOKEN)
        for station in wired["stations"]:
            assert "height_above_sea" in station
            assert station["attributes"][0]["name"] == "stationId"
            assert set(station["attributes"][2]["value"]) == {"latitude", "longitude"}


@pytest.mark.unit
class TestGetSpatialWindPassesResolutionAsADirection:
    """B250: see the module docstring."""

    @pytest.fixture()
    def wired(self, monkeypatch):
        from hera.measurements.GIS.raster.landcover import LandCoverToolkit
        from hera.simulations.windProfile.toolkit import WindProfileToolkit

        seen = {}
        monkeypatch.setattr(
            LandCoverToolkit, "getLandCover", lambda self, *a, **k: "landcover-grid"
        )

        def _getRoughness(
            self, landcover, windMeteorologicalDirection=None, resolution=None,
            **kwargs
        ):
            seen["windMeteorologicalDirection"] = windMeteorologicalDirection
            seen["resolution"] = resolution
            return "roughness-grid"

        monkeypatch.setattr(
            LandCoverToolkit, "getRoughnessFromLandcover", _getRoughness
        )
        monkeypatch.setattr(
            WindProfileToolkit, "_getWindSpeedDirection", lambda self, s, t: []
        )
        return seen

    @pytest.mark.xfail(
        strict=True,
        reason="B250: getSpatialWind calls getRoughnessFromLandcover(xarray, "
               "dxdy), but that method's second parameter is "
               "windMeteorologicalDirection and its third is resolution. The grid "
               "spacing is therefore consumed as a wind direction -- 30 degrees "
               "for the default dxdy=30 -- and the resolution it was meant to set "
               "stays None. See the consolidated findings issue.",
    )
    def test_the_resolution_is_passed_as_the_resolution(self, tk, wired):
        tk.getSpatialWind(32.0, 35.0, 33.0, 36.0, TOKEN, dxdy=50)
        assert wired["resolution"] == 50
        assert wired["windMeteorologicalDirection"] is None

    def test_the_resolution_currently_arrives_as_a_wind_direction(self, tk, wired):
        """Characterisation of B250."""
        tk.getSpatialWind(32.0, 35.0, 33.0, 36.0, TOKEN, dxdy=50)
        assert wired["windMeteorologicalDirection"] == 50
        assert wired["resolution"] is None

    def test_the_signature_it_calls_into(self):
        """Characterisation of B250's mechanism, so the diagnosis does not
        rest on the stub."""
        import inspect

        from hera.measurements.GIS.raster.landcover import LandCoverToolkit

        parameters = list(
            inspect.signature(LandCoverToolkit.getRoughnessFromLandcover).parameters
        )
        assert parameters[:4] == [
            "self", "landcover", "windMeteorologicalDirection", "resolution",
        ]
