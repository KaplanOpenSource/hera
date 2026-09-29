"""Angle conversions.

The two conventions, per CLAUDE.md:

* meteorological — 0 deg = North, increasing clockwise.  Used for wind
  direction input.
* mathematical  — 0 deg = East, increasing counter-clockwise.  Used
  internally for geometry.

Every assertion below is derived from those definitions, not from what
hera/utils/angle.py happens to return.
"""
import pytest

from hera.utils import toAzimuthAngle, toMathematicalAngle, toMeteorologicalAngle


@pytest.mark.unit
class TestMeteorologicalConversion:
    """(270 - x) mod 360 is the transform between the two conventions."""

    @pytest.mark.parametrize(
        "mathematical, meteorological",
        [
            (0, 270),    # East      -> wind from the West
            (90, 180),   # North     -> wind from the South
            (180, 90),   # West      -> wind from the East
            (270, 0),    # South     -> wind from the North
        ],
    )
    def test_cardinal_directions(self, mathematical, meteorological):
        assert toMeteorologicalAngle(mathematical) == meteorological

    def test_result_is_always_in_range(self):
        for angle in range(-720, 1080, 7):
            assert 0 <= toMeteorologicalAngle(angle) < 360

    def test_is_an_involution(self):
        """Applying the conversion twice must return the original angle.

        This is why one function can serve both directions -- see
        TestConversionAliasing.  A test is the only thing stopping someone
        from "fixing" one of the two names and silently breaking the other.
        """
        for angle in range(0, 360, 3):
            assert toMeteorologicalAngle(toMeteorologicalAngle(angle)) == angle


@pytest.mark.unit
class TestConversionAliasing:
    """toMathematicalAngle and toMeteorologicalAngle are one object.

    Mathematically sound: (270 - x) mod 360 is self-inverse, so the same
    transform converts in both directions.  It is undocumented though, and it
    means the two names give no protection against converting the wrong way --
    the precise risk CLAUDE.md calls out when it requires explicit conversion.
    """

    def test_the_two_names_are_the_same_object(self):
        assert toMathematicalAngle is toMeteorologicalAngle

    def test_round_trip_between_the_conventions(self):
        """Whatever the aliasing, a round trip must be identity."""
        for angle in (0, 45, 90, 137, 180, 275, 359):
            assert toMathematicalAngle(toMeteorologicalAngle(angle)) == angle


@pytest.mark.unit
class TestAzimuthConversion:
    """Azimuth: 0 deg = North, increasing clockwise, always in [0, 360)."""

    @pytest.mark.parametrize(
        "mathematical, azimuth",
        [
            (90, 0),     # North
            (0, 90),     # East
            (270, 180),  # South
            (180, 270),  # West
        ],
    )
    def test_cardinal_directions(self, mathematical, azimuth):
        assert toAzimuthAngle(mathematical) == azimuth

    def test_in_range_for_a_single_cycle(self):
        for angle in range(0, 360):
            assert 0 <= toAzimuthAngle(angle) < 360

    @pytest.mark.xfail(
        strict=True,
        reason="B1: toAzimuthAngle has no % 360 before its branch, so any input "
               "outside the first cycle yields an invalid azimuth "
               "(az(800) = -350, az(-400) = 490). See the consolidated findings issue.",
    )
    @pytest.mark.parametrize("angle", [800, -400, 1000, -721])
    def test_normalises_input_outside_one_cycle(self, angle):
        """An azimuth must be in [0, 360) for ANY real input, not just [0, 360)."""
        assert 0 <= toAzimuthAngle(angle) < 360
