"""LandCoverToolkit.roughnesslength2sandgrainroughness: Ks = z0 * 30 (Desmond et al. 2017, eq. 5).

The method is defined twice in the class body (once around line 561, once
around line 717) -- both @staticmethod, both computing the same formula.
The second definition simply shadows the first in Python's normal
class-body semantics; the earlier one is dead code, but the behaviour
itself is not affected, so this is noted rather than filed as a defect.
"""
import pytest

from hera.measurements.GIS.raster.landcover import LandCoverToolkit


@pytest.mark.unit
class TestRoughnessLength2SandGrainRoughness:
    def test_it_multiplies_by_thirty(self):
        assert LandCoverToolkit.roughnesslength2sandgrainroughness(0.1) == pytest.approx(3.0)

    def test_zero_roughness_length_gives_zero(self):
        assert LandCoverToolkit.roughnesslength2sandgrainroughness(0.0) == 0.0

    def test_it_is_callable_without_an_instance(self):
        """Declared @staticmethod, so the class itself can call it."""
        assert LandCoverToolkit.roughnesslength2sandgrainroughness(1.0) == pytest.approx(30.0)
