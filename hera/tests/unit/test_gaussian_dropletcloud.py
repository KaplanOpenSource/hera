"""gaussian/DropletCloud.py: FixedPositionDropletsCloud and its subclasses
all construct FallingNonEvaporatingDroplets internally, so B96 (already
documented: the constructor always raises AttributeError, see
FallingNonEvaporatingDroplets) cascades to every class in this file --
none of them can be constructed either.

B100: FixedPointClippedDropletCloud.__init__ calls super().__init__()
(which dispatches to this subclass's own _initDropletPosition override,
using self.clippedDiameter) before setting self.clippedDiameter -- so it
crashes with a numpy TypeError reading a still-None clippedDiameter,
before ever reaching the B96 crash its siblings hit cleanly.
"""
import pytest

from hera.simulations.gaussian.DropletCloud import (
    CirclePositionClippedDropletsCloud,
    FixedPointClippedDropletCloud,
    FixedPositionDropletsCloud,
    LinePositionDropletsCloud,
)
from hera.utils import ureg


@pytest.mark.unit
class TestDropletCloudCascadesB96:
    def test_fixed_position_cannot_be_constructed(self):
        with pytest.raises(AttributeError, match="getMeteorology"):
            FixedPositionDropletsCloud(
                mmd=200 * ureg.um, geometricstd=1.5,
                position=(0 * ureg.m, 0 * ureg.m, 10 * ureg.m), Q=1 * ureg.kg, clouds=3,
            )

    def test_line_position_cannot_be_constructed(self):
        with pytest.raises(AttributeError, match="getMeteorology"):
            LinePositionDropletsCloud(
                mmd=200 * ureg.um, geometricstd=1.5,
                position=(0 * ureg.m, 0 * ureg.m, 10 * ureg.m), Q=1 * ureg.kg,
                linelength=10 * ureg.m, clouds=3, linepositions=2,
            )

    @pytest.mark.xfail(
        strict=True,
        reason="B100: __init__ calls super().__init__() (which dispatches "
               "to this subclass's _initDropletPosition override) before "
               "setting self.clippedDiameter -- so it reads clippedDiameter "
               "while it is still None, raising a numpy TypeError before "
               "ever reaching the B96 crash a sibling class hits. "
               "See the consolidated findings issue.",
    )
    def test_clipped_hits_b96_like_its_siblings(self):
        with pytest.raises(AttributeError, match="getMeteorology"):
            FixedPointClippedDropletCloud(
                mmd=200 * ureg.um, geometricstd=1.5,
                position=(0 * ureg.m, 0 * ureg.m, 10 * ureg.m), Q=1 * ureg.kg,
                clippedDiameter=1 * ureg.mm, clouds=3,
            )

    def test_clipped_currently_raises_earlier_for_a_different_reason(self):
        """Characterisation of B100."""
        with pytest.raises(TypeError, match="log"):
            FixedPointClippedDropletCloud(
                mmd=200 * ureg.um, geometricstd=1.5,
                position=(0 * ureg.m, 0 * ureg.m, 10 * ureg.m), Q=1 * ureg.kg,
                clippedDiameter=1 * ureg.mm, clouds=3,
            )

    def test_circle_position_clipped_cannot_be_constructed(self):
        with pytest.raises(AttributeError, match="getMeteorology"):
            CirclePositionClippedDropletsCloud(
                mmd=200 * ureg.um, geometricstd=1.5,
                position=(0 * ureg.m, 0 * ureg.m, 10 * ureg.m), Q=1 * ureg.kg,
                clippedDiameter=1 * ureg.mm, clouds=3, circlepositions=2,
            )


@pytest.mark.unit
class TestClippedDiameterProperty:
    """The one piece of pure logic in this file that doesn't need
    construction: the clippedDiameter setter, exercised via __new__."""

    def test_a_string_is_parsed_as_a_pint_expression(self):
        obj = FixedPointClippedDropletCloud.__new__(FixedPointClippedDropletCloud)
        obj.clippedDiameter = "5 mm"
        assert obj.clippedDiameter.m_as(ureg.mm) == pytest.approx(5.0)

    def test_a_quantity_is_converted_to_millimeters(self):
        obj = FixedPointClippedDropletCloud.__new__(FixedPointClippedDropletCloud)
        obj.clippedDiameter = 1 * ureg.cm
        assert obj.clippedDiameter.m_as(ureg.mm) == pytest.approx(10.0)
