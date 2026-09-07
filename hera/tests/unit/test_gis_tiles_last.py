"""TilesToolkit: the slippy-map tile *fetching* half of the toolkit.

``test_gis_tiles.py`` covers the pure Web Mercator arithmetic
(``deg2tile``/``tile2deg``/``tileScaleAtLatLonZoom``) and the config
plumbing, and leaves everything that talks to a tile server for
"integration tests".  That is the majority of the module, and none of it
actually needs a server: the only outward call is a single
``requests.get`` per tile.  This file patches that one seam -- on the
``requests`` module, not on any toolkit instance -- and hands back real
PNG bytes built with PIL, so what is under test is the toolkit's own
work: which zoom/x/y each URL carries, how many tiles are fetched, where
each tile is pasted in the mosaic, how the extent is projected back into
the caller's CRS, and which branch each argument shape selects.

Nothing here opens a socket (the unit conftest makes that an error
anyway), reads a real raster, or contacts MongoDB beyond the in-memory
database the fixtures provide.

Covered here
------------
* ``TilesToolkit.getImageFromCorners`` -- tile-server resolution (explicit
  argument, project config, registered datasource), the input/output CRS
  round trips and the returned extent;
* ``TilesToolkit._getImageFromTiles`` -- squaring, mosaic assembly, URL
  formatting and the returned WGS84 bounds;
* ``TilesToolkit.listImages`` -- the doctype-scoped Measurements query;
* ``presentation.datalayer`` and ``presentation.plot`` -- every input
  shape (name / (image, extent) tuple / bare array with list, ``minX``
  dict, ``left`` dict or rubbish extents) and both output-CRS branches.

Expected numbers are the published OpenStreetMap slippy-map values, not
hera's output: at zoom 10 the bounding box 34.7E..35.2E / 31.9N..32.3N
has its north-west corner in tile (610, 414) and its south-east corner in
tile (612, 416), and tile x's western edge is at ``x / 2**z * 360 - 180``
degrees.

Deliberately not covered
------------------------
* ``deg2tile``/``tile2deg``/``tileScaleAtLatLonZoom``/``doctype``/
  ``setDefaultTileServer`` -- already in ``test_gis_tiles.py``;
* real tile imagery: the fake server returns flat single-colour PNGs, so
  the assertions are about *which* tile landed *where*, never about
  pixel content beyond that.

Bugs pinned here (each a strict xfail for the intended behaviour plus a
passing characterisation of today's behaviour):

* B314: ``_getImageFromTiles`` sizes the mosaic from the *difference*
  of the corner tile indices instead of the inclusive count
  (``lr - ul``, not ``lr - ul + 1``), so the eastern-most and
  southern-most tile row of the requested box is always dropped -- and a
  box that lies inside a single tile produces a 0x0 image with a
  zero-width extent and no HTTP request at all.
* B315: ``presentation.plot`` looks up a named image with
  ``self.datalayer.getImage(...)``, but no ``getImage`` method exists
  anywhere in hera, so the documented "name of datasource image in DB"
  input always raises AttributeError.
* B322: ``presentation.plot`` builds the matplotlib extent as
  ``[minX, maxX, minY, maxY]`` only when ``outputCRS`` is ITM; for every
  other CRS -- including its own ``WSG84`` default for ``inputCRS`` -- it
  emits ``[minY, maxY, minX, maxX]``, so the image is stretched over a
  transposed box.
* B323: ``_getImageFromTiles`` writes the squared corner back into the
  ``lrTiles`` list it was handed, mutating the caller's argument.
"""
import io

import pytest
from PIL import Image

from hera import toolkitHome
from hera.measurements.GIS import ITM, WSG84

TEMPLATE = "https://tiles.invalid/{z}/{x}/{y}.png"

# The published slippy-map tile numbers for the box used throughout.
ZOOM = 10
BOX = dict(minx=34.7, miny=31.9, maxx=35.2, maxy=32.3)
UL_TILE = (610, 414)   # north-west corner: (34.7E, 32.3N)
SE_TILE = (612, 416)   # south-east corner: (35.2E, 31.9N)

# The same two corners expressed in Israeli TM (EPSG:2039).
BOX_ITM = dict(minx=171945.12, miny=645229.69, maxx=219036.24, maxy=689697.39)


def _tileWestEdge(xtile, zoom=ZOOM):
    """Longitude of a tile's western edge -- the slippy-map inverse."""
    return xtile / 2.0**zoom * 360.0 - 180.0


def _png(colour):
    """A real 256x256 single-colour PNG, as bytes."""
    buffer = io.BytesIO()
    Image.new("RGB", (256, 256), colour).save(buffer, format="PNG")
    return buffer.getvalue()


class _Response:
    def __init__(self, content):
        self.content = content


@pytest.fixture()
def tileServer(monkeypatch):
    """Patch ``requests.get`` and record every URL asked for.

    Each tile comes back a different colour, keyed by its x/y, so the
    mosaic can be read back pixel by pixel.  The patch is on the
    ``requests`` module object -- the name ``tiles.py`` resolves at call
    time -- so monkeypatch restores it cleanly.
    """
    import requests

    calls = []

    def fake_get(url, *args, **kwargs):
        calls.append(url)
        tail = url.rsplit("/", 2)
        x, y = int(tail[-2]), int(tail[-1].split(".")[0])
        return _Response(_png((x % 256, y % 256, 7)))

    monkeypatch.setattr(requests, "get", fake_get)
    return calls


@pytest.fixture()
def tiles(unit_toolkit_factory):
    return unit_toolkit_factory(toolkitHome.GIS_TILES)


def _requestedTiles(calls):
    """The (z, x, y) triples the toolkit asked the server for."""
    out = []
    for url in calls:
        parts = url.split("/")
        out.append((int(parts[-3]), int(parts[-2]), int(parts[-1].split(".")[0])))
    return sorted(out)


# ---------------------------------------------------------------------------
# getImageFromCorners
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetImageFromCornersTileServerResolution:
    def test_with_no_server_anywhere_it_refuses(self, tiles, tileServer):
        with pytest.raises(ValueError, match="tile server"):
            tiles.getImageFromCorners(zoomlevel=ZOOM, **BOX)
        assert tileServer == []

    def test_the_project_config_supplies_the_default_server(self, tiles, tileServer):
        tiles.setDefaultTileServer(TEMPLATE)
        tiles.getImageFromCorners(zoomlevel=ZOOM, **BOX)
        assert tileServer
        assert all(url.startswith("https://tiles.invalid/") for url in tileServer)

    def test_an_explicit_server_overrides_the_configured_default(self, tiles, tileServer):
        tiles.setDefaultTileServer(TEMPLATE)
        tiles.getImageFromCorners(
            zoomlevel=ZOOM, tileServer="https://other.invalid/{z}/{x}/{y}.png", **BOX
        )
        assert all(url.startswith("https://other.invalid/") for url in tileServer)

    def test_a_registered_datasource_name_is_resolved_to_its_url_template(
        self, tiles, tileServer
    ):
        tiles.addDataSource("MYSERVER", "https://fromdb.invalid/{z}/{x}/{y}.png",
                            "string", version=(0, 0, 1))
        tiles.getImageFromCorners(zoomlevel=ZOOM, tileServer="MYSERVER", **BOX)
        assert tileServer
        assert all(url.startswith("https://fromdb.invalid/") for url in tileServer)

    def test_an_unknown_name_is_used_as_the_url_template_itself(self, tiles, tileServer):
        """getDataSourceData returns None for an unknown name, and the string
        is then treated as the template."""
        tiles.getImageFromCorners(zoomlevel=ZOOM, tileServer=TEMPLATE, **BOX)
        assert tileServer


@pytest.mark.unit
class TestGetImageFromCornersRequestedTiles:
    def test_the_zoom_level_is_passed_through_to_every_url(self, tiles, tileServer):
        tiles.getImageFromCorners(zoomlevel=ZOOM, tileServer=TEMPLATE, **BOX)
        assert {z for z, _, _ in _requestedTiles(tileServer)} == {ZOOM}

    def test_the_north_west_tile_of_the_box_is_the_first_one_requested(
        self, tiles, tileServer
    ):
        tiles.getImageFromCorners(zoomlevel=ZOOM, tileServer=TEMPLATE, **BOX)
        first = _requestedTiles(tileServer)[0]
        assert first == (ZOOM, UL_TILE[0], UL_TILE[1])

    def test_every_tile_requested_lies_inside_the_requested_box(self, tiles, tileServer):
        tiles.getImageFromCorners(zoomlevel=ZOOM, tileServer=TEMPLATE, **BOX)
        for _, x, y in _requestedTiles(tileServer):
            assert UL_TILE[0] <= x <= SE_TILE[0]
            assert UL_TILE[1] <= y <= SE_TILE[1]

    def test_each_tile_is_fetched_exactly_once(self, tiles, tileServer):
        tiles.getImageFromCorners(zoomlevel=ZOOM, tileServer=TEMPLATE, **BOX)
        assert len(tileServer) == len(set(tileServer))

    def test_itm_input_coordinates_reach_the_same_tiles_as_their_wgs84_twins(
        self, tiles, tileServer
    ):
        tiles.getImageFromCorners(zoomlevel=ZOOM, tileServer=TEMPLATE,
                                  inputCRS=ITM, **BOX_ITM)
        assert _requestedTiles(tileServer)[0] == (ZOOM, UL_TILE[0], UL_TILE[1])


@pytest.mark.unit
class TestGetImageFromCornersReturnValue:
    def test_it_returns_a_pillow_image_and_a_four_element_extent(self, tiles, tileServer):
        image, extent = tiles.getImageFromCorners(
            zoomlevel=ZOOM, tileServer=TEMPLATE, **BOX
        )
        assert isinstance(image, Image.Image)
        assert len(extent) == 4

    def test_the_mosaic_is_a_whole_number_of_256_pixel_tiles(self, tiles, tileServer):
        image, _ = tiles.getImageFromCorners(zoomlevel=ZOOM, tileServer=TEMPLATE, **BOX)
        assert image.size[0] % 256 == 0
        assert image.size[1] % 256 == 0
        assert image.size[0] * image.size[1] // (256 * 256) == len(tileServer)

    def test_the_wgs84_extent_is_left_right_bottom_top(self, tiles, tileServer):
        _, extent = tiles.getImageFromCorners(zoomlevel=ZOOM, tileServer=TEMPLATE, **BOX)
        assert extent[0] < extent[1]
        assert extent[2] < extent[3]

    def test_the_extent_starts_at_the_western_edge_of_the_north_west_tile(
        self, tiles, tileServer
    ):
        _, extent = tiles.getImageFromCorners(zoomlevel=ZOOM, tileServer=TEMPLATE, **BOX)
        assert extent[0] == pytest.approx(_tileWestEdge(UL_TILE[0]))

    def test_the_extent_is_projected_into_the_requested_output_crs(
        self, tiles, tileServer
    ):
        _, degrees = tiles.getImageFromCorners(
            zoomlevel=ZOOM, tileServer=TEMPLATE, **BOX
        )
        _, metres = tiles.getImageFromCorners(
            zoomlevel=ZOOM, tileServer=TEMPLATE, outputCRS=ITM, **BOX
        )
        # Israeli TM eastings over Israel are ~1.2e5..2.5e5 metres.
        assert 1e5 < metres[0] < 3e5
        assert metres[0] < metres[1]
        assert metres[2] < metres[3]
        assert abs(degrees[0]) < 180


# ---------------------------------------------------------------------------
# _getImageFromTiles
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetImageFromTiles:
    def test_the_url_template_receives_z_x_and_y(self, tiles, tileServer):
        tiles._getImageFromTiles([610, 414], [612, 416], 10, TEMPLATE)
        assert "https://tiles.invalid/10/610/414.png" in tileServer

    def test_a_two_by_two_span_is_a_512_pixel_square(self, tiles, tileServer):
        image, _ = tiles._getImageFromTiles([610, 414], [612, 416], 10, TEMPLATE)
        assert image.size == (512, 512)

    def test_squaring_extends_the_short_side_to_the_long_one(self, tiles, tileServer):
        image, _ = tiles._getImageFromTiles([610, 414], [613, 415], 10, TEMPLATE)
        assert image.size == (256 * 3, 256 * 3)

    def test_without_squaring_each_side_keeps_its_own_span(self, tiles, tileServer):
        image, _ = tiles._getImageFromTiles(
            [610, 414], [613, 415], 10, TEMPLATE, square=False
        )
        assert image.size == (256 * 3, 256 * 1)

    def test_each_tile_is_pasted_at_its_own_offset(self, tiles, tileServer):
        """The fake server colours tile (x, y) as (x%256, y%256, 7)."""
        image, _ = tiles._getImageFromTiles([610, 414], [612, 416], 10, TEMPLATE)
        assert image.getpixel((10, 10)) == (610 % 256, 414 % 256, 7)
        assert image.getpixel((256 + 10, 10)) == (611 % 256, 414 % 256, 7)
        assert image.getpixel((10, 256 + 10)) == (610 % 256, 415 % 256, 7)

    def test_the_returned_bounds_run_west_to_east_and_north_to_south(
        self, tiles, tileServer
    ):
        _, extent = tiles._getImageFromTiles([610, 414], [612, 416], 10, TEMPLATE)
        left, right, top, bottom = extent
        assert left < right
        assert top > bottom

    def test_the_bounds_are_the_tile_edges_of_the_assembled_mosaic(
        self, tiles, tileServer
    ):
        _, extent = tiles._getImageFromTiles([610, 414], [612, 416], 10, TEMPLATE)
        assert extent[0] == pytest.approx(_tileWestEdge(610))
        assert extent[1] == pytest.approx(_tileWestEdge(612))

    def test_a_larger_zoom_level_reaches_a_smaller_span_of_degrees(
        self, tiles, tileServer
    ):
        _, coarse = tiles._getImageFromTiles([610, 414], [612, 416], 10, TEMPLATE)
        _, fine = tiles._getImageFromTiles([1220, 828], [1222, 830], 11, TEMPLATE)
        assert (fine[1] - fine[0]) < (coarse[1] - coarse[0])


# ---------------------------------------------------------------------------
# B314: the tile count is a difference, not an inclusive count
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestTheRequestedBoxIsFullyCovered:
    @pytest.mark.xfail(
        strict=True,
        reason="B314: _getImageFromTiles sizes the mosaic as lr - ul rather "
               "than lr - ul + 1, so the eastern-most and southern-most tile "
               "row of the requested box is never fetched. Requesting "
               "34.7E..35.2E / 31.9N..32.3N at zoom 10 spans tiles 610..612 x "
               "414..416 -- nine tiles -- but only four are fetched. See the "
               "consolidated findings issue.",
    )
    def test_every_tile_touching_the_box_is_fetched(self, tiles, tileServer):
        tiles.getImageFromCorners(zoomlevel=ZOOM, tileServer=TEMPLATE, **BOX)
        expected = {
            (ZOOM, x, y)
            for x in range(UL_TILE[0], SE_TILE[0] + 1)
            for y in range(UL_TILE[1], SE_TILE[1] + 1)
        }
        assert set(_requestedTiles(tileServer)) == expected

    def test_only_the_interior_tiles_are_fetched_today(self, tiles, tileServer):
        """Characterisation of B314."""
        tiles.getImageFromCorners(zoomlevel=ZOOM, tileServer=TEMPLATE, **BOX)
        assert set(_requestedTiles(tileServer)) == {
            (ZOOM, 610, 414), (ZOOM, 610, 415),
            (ZOOM, 611, 414), (ZOOM, 611, 415),
        }

    @pytest.mark.xfail(
        strict=True,
        reason="B314: a box that lies inside a single tile makes the corner "
               "tile indices equal, so lr - ul is zero and the toolkit builds "
               "a 0x0 image without asking the server for anything. A caller "
               "zoomed in on one street gets an empty image and no error. See "
               "the consolidated findings issue.",
    )
    def test_a_box_inside_a_single_tile_still_returns_that_tile(self, tiles, tileServer):
        image, _ = tiles.getImageFromCorners(
            minx=34.7805, miny=32.0800, maxx=34.7810, maxy=32.0805,
            zoomlevel=14, tileServer=TEMPLATE,
        )
        assert image.size == (256, 256)

    def test_a_box_inside_a_single_tile_yields_an_empty_image_today(
        self, tiles, tileServer
    ):
        """Characterisation of B314."""
        image, extent = tiles.getImageFromCorners(
            minx=34.7805, miny=32.0800, maxx=34.7810, maxy=32.0805,
            zoomlevel=14, tileServer=TEMPLATE,
        )
        assert image.size == (0, 0)
        assert tileServer == []
        assert extent[0] == pytest.approx(extent[1])
        assert extent[2] == pytest.approx(extent[3])


# ---------------------------------------------------------------------------
# B323: the lrTiles argument is mutated
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestTheCornerArgumentsAreLeftAlone:
    @pytest.mark.xfail(
        strict=True,
        reason="B323: _getImageFromTiles assigns the squared corner back "
               "into the caller's list (`lrTiles[0] = ulTiles[0] + sqrx`), so "
               "the argument is mutated in place. A caller that reuses its "
               "corner list -- to fetch the same box at two zoom levels, say "
               "-- silently gets a different box the second time. See the "
               "consolidated findings issue.",
    )
    def test_the_lower_right_corner_list_survives_the_call(self, tiles, tileServer):
        lower_right = [613, 415]
        tiles._getImageFromTiles([610, 414], lower_right, 10, TEMPLATE)
        assert lower_right == [613, 415]

    def test_the_lower_right_corner_is_overwritten_with_the_squared_one(
        self, tiles, tileServer
    ):
        """Characterisation of B323."""
        lower_right = [613, 415]
        tiles._getImageFromTiles([610, 414], lower_right, 10, TEMPLATE)
        assert lower_right == [613, 417]


# ---------------------------------------------------------------------------
# listImages
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestListImages:
    def test_an_empty_project_lists_nothing(self, tiles):
        assert list(tiles.listImages()) == []

    def test_a_stored_image_document_is_listed(self, tiles, tmp_path):
        tiles.addMeasurementsDocument(
            resource=str(tmp_path / "a.png"), dataFormat="string",
            type=tiles.doctype, desc=dict(imageName="a"),
        )
        assert len(tiles.listImages()) == 1

    def test_documents_of_another_type_are_not_listed(self, tiles, tmp_path):
        tiles.addMeasurementsDocument(
            resource=str(tmp_path / "b.png"), dataFormat="string",
            type="SomethingElse", desc=dict(imageName="b"),
        )
        assert list(tiles.listImages()) == []

    def test_extra_filters_narrow_the_result(self, tiles, tmp_path):
        for name in ("a", "b"):
            tiles.addMeasurementsDocument(
                resource=str(tmp_path / f"{name}.png"), dataFormat="string",
                type=tiles.doctype, desc=dict(imageName=name),
            )
        assert len(tiles.listImages()) == 2
        assert len(tiles.listImages(imageName="a")) == 1

    def test_a_filter_matching_nothing_gives_an_empty_list(self, tiles, tmp_path):
        tiles.addMeasurementsDocument(
            resource=str(tmp_path / "a.png"), dataFormat="string",
            type=tiles.doctype, desc=dict(imageName="a"),
        )
        assert list(tiles.listImages(imageName="nosuch")) == []


# ---------------------------------------------------------------------------
# presentation
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestPresentationDatalayer:
    def test_the_presentation_layer_points_back_at_its_toolkit(self, tiles):
        assert tiles.presentation.datalayer is tiles

    def test_the_presentation_layer_is_not_a_toolkit_subclass(self, tiles):
        assert not isinstance(tiles.presentation, type(tiles))


def _image():
    import numpy

    return numpy.arange(12, dtype=float).reshape(3, 4)


@pytest.mark.unit
class TestPresentationPlotInputShapes:
    def test_a_tuple_of_image_and_extent_is_accepted(self, tiles):
        result = tiles.presentation.plot(
            (_image(), [34.7, 35.2, 31.9, 32.3]), outputCRS=ITM
        )
        import matplotlib.image

        assert isinstance(result, matplotlib.image.AxesImage)

    def test_a_tuple_together_with_extents_is_refused(self, tiles):
        with pytest.raises(ValueError, match="extents must be None"):
            tiles.presentation.plot(
                (_image(), [34.7, 35.2, 31.9, 32.3]), extents=[0, 1, 0, 1]
            )

    def test_a_bare_image_without_extents_is_refused(self, tiles):
        with pytest.raises(ValueError, match="extents must be supplied"):
            tiles.presentation.plot(_image())

    def test_a_bare_image_with_a_list_extent_is_accepted(self, tiles):
        result = tiles.presentation.plot(
            _image(), extents=[34.7, 35.2, 31.9, 32.3], outputCRS=ITM
        )
        assert result.get_extent()[0] < result.get_extent()[1]

    def test_a_minx_dictionary_extent_is_accepted(self, tiles):
        byList = tiles.presentation.plot(
            _image(), extents=[34.7, 35.2, 31.9, 32.3], outputCRS=ITM
        ).get_extent()
        byDict = tiles.presentation.plot(
            _image(),
            extents=dict(minX=34.7, maxX=35.2, minY=31.9, maxY=32.3),
            outputCRS=ITM,
        ).get_extent()
        assert byDict == pytest.approx(byList)

    def test_a_left_right_bottom_top_dictionary_extent_is_accepted(self, tiles):
        byList = tiles.presentation.plot(
            _image(), extents=[34.7, 35.2, 31.9, 32.3], outputCRS=ITM
        ).get_extent()
        byDict = tiles.presentation.plot(
            _image(),
            extents=dict(left=34.7, right=35.2, bottom=31.9, top=32.3),
            outputCRS=ITM,
        ).get_extent()
        assert byDict == pytest.approx(byList)

    def test_an_extent_that_is_neither_list_nor_dict_is_refused(self, tiles):
        with pytest.raises(ValueError, match="extents is either a list"):
            tiles.presentation.plot(_image(), extents=(34.7, 35.2, 31.9, 32.3))


@pytest.mark.unit
class TestPresentationPlotProjection:
    def test_the_itm_extent_is_in_metres(self, tiles):
        extent = tiles.presentation.plot(
            _image(), extents=[34.7, 35.2, 31.9, 32.3], outputCRS=ITM
        ).get_extent()
        assert 1e5 < extent[0] < 3e5
        assert extent[0] < extent[1]
        assert extent[2] < extent[3]

    def test_an_explicit_axis_is_used_rather_than_a_new_figure(self, tiles):
        import matplotlib.pyplot as plt

        figure, axis = plt.subplots()
        before = len(plt.get_fignums())
        tiles.presentation.plot(
            _image(), extents=[34.7, 35.2, 31.9, 32.3], outputCRS=ITM, ax=axis
        )
        assert len(plt.get_fignums()) == before

    def test_with_display_off_no_figure_is_left_open(self, tiles):
        import matplotlib.pyplot as plt

        plt.close("all")
        result = tiles.presentation.plot(
            _image(), extents=[34.7, 35.2, 31.9, 32.3], outputCRS=ITM, display=False
        )
        assert plt.get_fignums() == []
        assert result.get_extent()[0] < result.get_extent()[1]


# ---------------------------------------------------------------------------
# B315: plot cannot look an image up by name
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestPresentationPlotByName:
    @pytest.mark.xfail(
        strict=True,
        reason="B315: presentation.plot's documented first argument is the "
               "'name of datasource image in DB', and it resolves it with "
               "self.datalayer.getImage(name), but no getImage method exists "
               "on TilesToolkit -- or anywhere in hera -- so the named-image "
               "branch always raises AttributeError. See the consolidated "
               "findings issue.",
    )
    def test_a_stored_image_can_be_plotted_by_name(self, tiles, tmp_path):
        import numpy

        target = tmp_path / "img.png"
        Image.fromarray(numpy.zeros((4, 4, 3), dtype="uint8")).save(target)
        tiles.addMeasurementsDocument(
            resource=str(target), dataFormat="image", type=tiles.doctype,
            desc=dict(imageName="img", minX=34.7, maxX=35.2, minY=31.9, maxY=32.3),
        )
        tiles.presentation.plot("img", outputCRS=ITM)

    def test_the_toolkit_has_no_getimage_method(self, tiles):
        """Characterisation of B315."""
        assert not hasattr(tiles, "getImage")

    def test_plotting_by_name_raises_attribute_error(self, tiles):
        """Characterisation of B315."""
        with pytest.raises(AttributeError, match="getImage"):
            tiles.presentation.plot("anyName", outputCRS=ITM)


# ---------------------------------------------------------------------------
# B322: the extent is transposed for every output CRS but ITM
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestPresentationPlotNonItmOutput:
    EXTENTS = [34.7, 35.2, 31.9, 32.3]

    @pytest.mark.xfail(
        strict=True,
        reason="B322: presentation.plot orders the matplotlib extent as "
               "[x_min, x_max, y_min, y_max] only under `if outputCRS==ITM`; "
               "every other output CRS takes the else branch, which emits "
               "[y_min, y_max, x_min, x_max]. Plotting in WGS84 therefore "
               "stretches the image over a transposed box -- latitudes on the "
               "x axis and longitudes on the y axis. See the consolidated "
               "findings issue.",
    )
    def test_a_wgs84_round_trip_keeps_the_extent_it_was_given(self, tiles):
        extent = tiles.presentation.plot(
            _image(), extents=self.EXTENTS, inputCRS=WSG84, outputCRS=WSG84
        ).get_extent()
        assert extent == pytest.approx(self.EXTENTS)

    def test_a_wgs84_round_trip_transposes_the_extent_today(self, tiles):
        """Characterisation of B322."""
        extent = tiles.presentation.plot(
            _image(), extents=self.EXTENTS, inputCRS=WSG84, outputCRS=WSG84
        ).get_extent()
        assert extent == pytest.approx([31.9, 32.3, 34.7, 35.2])
