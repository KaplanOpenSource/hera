"""machineLearningDeepLearning/LightingModule.py: LightningModuleHera is a
plain mixin class (no torch base), so its container-holding setter is
testable without any real model/dataset."""
import pytest

from hera.simulations.machineLearningDeepLearning.LightingModule import LightningModuleHera


@pytest.mark.unit
class TestSetModelContainer:
    def test_it_defaults_to_none(self):
        assert LightningModuleHera().modelContainer is None

    def test_it_stores_whatever_is_passed(self):
        module = LightningModuleHera()
        module.setModelContainer("some-container")
        assert module.modelContainer == "some-container"
