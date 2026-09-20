"""simulations/gaussian/MeshUtils.py: GaussianToMesh.dxdy, its only
externally-visible state not covered by grid-construction tests
elsewhere."""
import pytest

from hera.simulations.gaussian.MeshUtils import GaussianToMesh


@pytest.mark.unit
class TestDxDy:
    def test_it_defaults_to_none(self):
        assert GaussianToMesh().dxdy is None

    def test_the_constructor_argument_is_stored(self):
        assert GaussianToMesh(dxdy=5.0).dxdy == 5.0

    def test_it_can_be_reassigned(self):
        mesh = GaussianToMesh()
        mesh.dxdy = 2.5
        assert mesh.dxdy == 2.5
