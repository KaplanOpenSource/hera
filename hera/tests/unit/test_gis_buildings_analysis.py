"""Blocks and analysis.ConvexPolygons: urban-canopy morphology parameters.

Splits a domain into a grid of blocks and computes, per block, the plan
area fraction (lambda_P), frontal area density (lambda_F) and the
area-weighted mean building height (h_c) -- the standard BNTL-style urban
canopy parameters used to drive roughness in dispersion models.

One defect surfaced in ``Blocks.getHc``:

* B79: it guards on ``if self._LambdaP is not None`` -- but ``_LambdaP``
  (no parentheses) is the *method itself*, a bound method object that is
  never ``None``. The guard can never take the ``else`` branch, so
  ``getHc()`` always returns ``self._hc`` (0 until ``_LambdaP()`` has
  actually run) regardless of whether the calculation happened. It is
  latent rather than currently harmful -- the one real caller,
  ``Lambda()``, always calls ``_LambdaP()`` before ``getHc()`` -- but the
  documented "returns None if lambda_P is not set" contract cannot happen.
"""
import geopandas
import pytest
from shapely.geometry import box

from hera.measurements.GIS.vector.buildings.analysis import Blocks, analysis


def _domain(size=200):
    return geopandas.GeoDataFrame({"geometry": [box(0, 0, size, size)]})


def _buildings(*rows):
    """rows: (xmin, ymin, xmax, ymax, height, ht_land, ftype, hi_pnt_z)."""
    geoms = [box(r[0], r[1], r[2], r[3]) for r in rows]
    return geopandas.GeoDataFrame({
        "geometry": geoms,
        "BLDG_HT": [r[4] for r in rows],
        "HT_LAND": [r[5] for r in rows],
        "FTYPE": [r[6] for r in rows],
        "HI_PNT_Z": [r[7] for r in rows],
    })


@pytest.mark.unit
class TestBlocksConstruction:
    def test_size_sets_a_fixed_grid_spacing_on_both_axes(self):
        blocks = Blocks(level=0, df=_domain(), size=200)
        assert blocks._Division == {"x": 200, "y": 200}
        assert blocks._DivisionType == {"x": "size", "y": "size"}

    def test_npxy_sets_a_fixed_block_count_on_both_axes(self):
        blocks = Blocks(level=0, df=_domain(), npxy=3)
        assert blocks._Division == {"x": 4, "y": 4}
        assert blocks._DivisionType == {"x": "count", "y": "count"}

    def test_size_and_npxy_together_raises(self):
        with pytest.raises(ValueError, match="size with npxy"):
            Blocks(level=0, df=_domain(), size=200, npxy=3)

    def test_width_and_npx_together_raises(self):
        with pytest.raises(ValueError, match="width with npx"):
            Blocks(level=0, df=_domain(), width=200, npx=3, height=200)

    def test_no_division_parameters_raises(self):
        with pytest.raises(ValueError, match="not set"):
            Blocks(level=0, df=_domain())

    def test_width_and_height_can_be_set_independently(self):
        blocks = Blocks(level=0, df=_domain(), width=150, height=50)
        assert blocks._Division == {"x": 150, "y": 50}


@pytest.mark.unit
class TestIterBlocksGrid:
    def test_a_200m_domain_split_by_100_gives_four_blocks(self):
        root = Blocks(level=0, df=_domain(200), size=100)
        level1 = root.iterBlocks(size=100)
        assert len(level1._listOfDicts) == 4

    def test_each_block_records_its_own_bounds(self):
        root = Blocks(level=0, df=_domain(200), size=100)
        level1 = root.iterBlocks(size=100)
        first = level1._listOfDicts[0]
        assert (first["xMin1"][0], first["yMin1"][0]) == (0.0, 0.0)
        assert (first["xMax1"][0], first["yMax1"][0]) == (100.0, 100.0)


@pytest.mark.unit
class TestLambdaPipeline:
    """One building per 100x100 block; verified against the hand-computed
    plan/frontal area fractions and area-weighted height."""

    @pytest.fixture()
    def blocks(self):
        root = Blocks(level=0, df=_domain(200), size=100)
        return root.iterBlocks(size=100)

    @pytest.fixture()
    def buildings(self):
        # 20x20 building (area 400) in the (0,0) block, height 10;
        # 20x20 building in the (1,1) block, height 15.
        return _buildings(
            (10, 10, 30, 30, 10.0, 5.0, 1, 15.0),
            (150, 150, 170, 170, 15.0, 5.0, 1, 20.0),
        )

    def test_lambda_p_is_the_plan_area_fraction_per_block(self, blocks, buildings):
        result = blocks.Lambda(buildings, windDirection=270)
        occupied = result[(result["i0"] == 0) & (result["j0"] == 0)].iloc[0]
        assert occupied["lambdaP"] == pytest.approx(400 / 10000)

    def test_lambda_f_is_the_frontal_area_fraction_per_block(self, blocks, buildings):
        result = blocks.Lambda(buildings, windDirection=270)
        occupied = result[(result["i0"] == 0) & (result["j0"] == 0)].iloc[0]
        # a 20x20 square rotated by a multiple of 90 degrees keeps its bounds
        assert occupied["lambdaF"] == pytest.approx(20 * 10 / 10000)

    def test_hc_is_the_building_height_when_only_one_building_is_present(self, blocks, buildings):
        result = blocks.Lambda(buildings, windDirection=270)
        occupied = result[(result["i0"] == 1) & (result["j0"] == 1)].iloc[0]
        assert occupied["hc"] == pytest.approx(15.0)

    def test_empty_blocks_get_zero_for_every_parameter(self, blocks, buildings):
        result = blocks.Lambda(buildings, windDirection=270)
        empty = result[(result["i0"] == 0) & (result["j0"] == 1)].iloc[0]
        assert (empty["lambdaP"], empty["lambdaF"], empty["hc"]) == (0, 0, 0)

    def test_ftype_16_buildings_are_excluded_from_lambda_p(self, blocks):
        buildings = _buildings((10, 10, 30, 30, 10.0, 5.0, 16, 15.0))
        result = blocks.Lambda(buildings, windDirection=270)
        occupied = result[(result["i0"] == 0) & (result["j0"] == 0)].iloc[0]
        assert occupied["lambdaP"] == 0

    def test_getHc_currently_returns_the_area_weighted_height_never_none(self, blocks, buildings):
        """Characterisation of B79."""
        block = Blocks(level=0, df=_domain(100), size=100)
        assert block.getHc() == 0
        block._hc = 42.0
        assert block.getHc() == 42.0


@pytest.mark.unit
class TestConvexPolygons:
    def test_nearby_buildings_are_grouped_into_one_convex_hull(self):
        data = geopandas.GeoDataFrame({"geometry": [box(0, 0, 10, 10), box(12, 0, 22, 10)]})
        result = analysis.ConvexPolygons(_FakeAnalysisSelf(), data, buffer=5)
        assert len(result) == 1

    def test_far_apart_buildings_stay_in_separate_groups(self):
        data = geopandas.GeoDataFrame({"geometry": [
            box(0, 0, 10, 10), box(12, 0, 22, 10),
            box(500, 500, 510, 510), box(512, 500, 522, 510),
        ]})
        result = analysis.ConvexPolygons(_FakeAnalysisSelf(), data, buffer=1)
        assert len(result) == 2

    def test_the_result_is_sorted_by_descending_area(self):
        data = geopandas.GeoDataFrame({"geometry": [
            box(0, 0, 5, 5),
            box(500, 500, 530, 530),
        ]})
        result = analysis.ConvexPolygons(_FakeAnalysisSelf(), data, buffer=1)
        assert list(result["area"]) == sorted(result["area"], reverse=True)


class _FakeAnalysisSelf:
    """ConvexPolygons only uses `self` for the recursive call; a bare
    instance lets us call it as an unbound function without a real
    BuildingsToolkit/datalayer behind it."""
