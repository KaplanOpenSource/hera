"""Aerodynamic roughness length from land-cover class.

z0 feeds straight into the wind profile and from there into every dispersion
calculation, so the numbers matter physically.  The module carries TWO tables
for the same IGBP codes: one cited to Floors et al. (2021), WES 6, 1379,
table a2, and one labelled in the source as "Example values".

The published values are the reference these tests use.
"""
import pathlib

import pytest

from hera.measurements.GIS.raster.landcover import LandCoverToolkit

# Floors et al. (2021), table a2 -- the table _handleType1 cites.
PUBLISHED = {
    0: 0.0001,   # water
    1: 1.0,      # evergreen needleleaf forest
    2: 1.0,      # evergreen broadleaf forest
    3: 1.0,      # deciduous needleleaf forest
    4: 1.0,      # deciduous broadleaf forest
    5: 1.0,      # mixed forest
    6: 0.05,     # closed shrubland
    7: 0.06,     # open shrubland
    8: 0.05,     # woody savanna
    9: 0.15,     # savanna
    10: 0.12,    # grassland
    11: 0.3,     # permanent wetland
    12: 0.15,    # cropland
    13: 0.8,     # urban and built-up
    14: 0.14,    # cropland / natural vegetation mosaic
    15: 0.001,   # snow and ice
    16: 0.01,    # barren or sparsely vegetated
}


@pytest.fixture()
def toolkit(unit_toolkit_factory):
    from hera import toolkitHome

    return unit_toolkit_factory(toolkitHome.GIS_LANDCOVER)


@pytest.mark.unit
class TestPublishedRoughnessTable:
    """_handleType1, checked entry by entry against the cited paper."""

    @pytest.mark.parametrize("code, expected", sorted(PUBLISHED.items()))
    def test_each_class_has_its_published_roughness(self, toolkit, code, expected):
        assert toolkit._handleType1(code) == pytest.approx(expected)

    def test_water_is_the_smoothest_surface(self, toolkit):
        """Open water is the physical floor for z0, by four orders of magnitude."""
        water = toolkit._handleType1(0)
        assert water == pytest.approx(0.0001)
        assert all(
            water <= toolkit._handleType1(code) for code in PUBLISHED if code != 0
        )

    def test_forest_is_the_roughest(self, toolkit):
        """Closed-canopy forest tops the scale at z0 = 1 m."""
        forest = toolkit._handleType1(1)
        assert forest == pytest.approx(1.0)
        assert all(toolkit._handleType1(code) <= forest for code in PUBLISHED)

    def test_urban_is_rougher_than_grassland(self, toolkit):
        assert toolkit._handleType1(13) > toolkit._handleType1(10)

    def test_snow_and_ice_are_nearly_as_smooth_as_water(self, toolkit):
        assert toolkit._handleType1(15) == pytest.approx(0.001)
        assert toolkit._handleType1(15) < toolkit._handleType1(16)

    def test_an_unknown_code_falls_back(self, toolkit):
        """Documented fallback of 0.05, roughly short vegetation."""
        assert toolkit._handleType1(999) == pytest.approx(0.05)

    def test_every_value_is_a_physical_roughness_length(self, toolkit):
        """z0 must be positive and, for natural surfaces, below a few metres."""
        for code in PUBLISHED:
            value = toolkit._handleType1(code)
            assert 0 < value <= 3.0


@pytest.mark.unit
class TestTheSecondTable:
    """The inline table in getRoughnessAtPoint disagrees with the published one."""

    SOURCE = pathlib.Path("hera/measurements/GIS/raster/landcover.py")

    # The ramp actually written in getRoughnessAtPoint.
    EXAMPLE_RAMP = {
        0: 0.01, 1: 0.1, 2: 0.15, 3: 0.2, 4: 0.25, 5: 0.3, 6: 0.4, 7: 0.5,
        8: 0.6, 9: 0.7, 10: 0.8, 11: 0.9, 12: 1.0, 13: 1.1, 14: 1.2,
        15: 1.3, 16: 1.4,
    }

    def test_the_source_labels_it_as_example_values(self):
        assert "# Example values" in self.SOURCE.read_text(encoding="utf-8")

    def test_it_is_an_arithmetic_ramp_not_a_physical_table(self):
        """Codes 2 to 16 step by a constant 0.05 or 0.1 -- land cover does not."""
        steps = [
            round(self.EXAMPLE_RAMP[code + 1] - self.EXAMPLE_RAMP[code], 4)
            for code in range(2, 16)
        ]
        assert set(steps) == {0.05, 0.1}

    def test_it_contradicts_the_published_table_on_water(self, toolkit):
        """0.01 against a published 0.0001 -- a factor of 100."""
        assert self.EXAMPLE_RAMP[0] == pytest.approx(100 * toolkit._handleType1(0))

    def test_it_contradicts_the_published_table_on_snow(self, toolkit):
        """1.3 against a published 0.001 -- a factor of 1300."""
        assert self.EXAMPLE_RAMP[15] / toolkit._handleType1(15) == pytest.approx(1300.0)

    def test_it_inverts_the_ordering_of_water_and_forest(self, toolkit):
        """The ramp makes forest smoother than snow, which is backwards."""
        assert self.EXAMPLE_RAMP[1] < self.EXAMPLE_RAMP[15]
        assert toolkit._handleType1(1) > toolkit._handleType1(15)

    @pytest.mark.xfail(
        strict=True,
        reason="B43: getRoughnessAtPoint carries its own IGBP table, labelled "
               "'# Example values', which is an arithmetic ramp rather than a "
               "physical lookup. It disagrees with the published table in "
               "_handleType1 by a factor of 100 for water and 1300 for snow and "
               "ice, and reverses the forest/snow ordering. z0 feeds the wind "
               "profile, so a caller reaching this branch gets a materially wrong "
               "simulation. See the consolidated findings issue.",
    )
    @pytest.mark.parametrize("code", sorted(PUBLISHED))
    def test_the_two_tables_agree(self, toolkit, code):
        assert self.EXAMPLE_RAMP[code] == pytest.approx(toolkit._handleType1(code))
