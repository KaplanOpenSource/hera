"""simulations/deposition/models.py: the canopy branch of the Petroff
deposition model, plus the property accessors themselves.

`test_deposition_models.py` (an earlier batch) already covers
construction, the tangled `ustar`/`diameter`/`heatFlux` properties
(B46-B49) and the dispatch, but only ever exercises the bare-ground path
-- so the nested `PhiM`/`PsiM` stability functions, which live inside the
`obstacleShape is not None` branch, were never reached. This file drives
that branch with a leaf canopy.

Note the tests deliberately pass `diameter=0.3`: because of B46 the
surviving `ustar` property returns `_diameter`, so the diameter is what
the model actually uses as friction velocity. Setting it to a realistic
u* keeps the canopy arithmetic inside the stability range PhiM/PsiM
accept, and is the only way to steer that value from outside.
"""
import pytest

DESERT = {"name": "Desert", "zrough": 0.04}
FOREST = {
    "name": "Forest", "zrough": 0.4, "hcanop": 10.0, "displacementHeight": 6.0,
    "obstacleShape": "leaf", "LAI": 3.0, "plagiophile": 0.5, "xb": 1.0,
    "xim": 1.0, "xin": 1.0, "xit": 1.0, "ObstacleSize": 0.05,
}


@pytest.fixture()
def make(unit_project):
    """depositionModels builds a real Project, so it must be pointed at
    the mongomock-backed test project, never its own default."""
    from hera.simulations.deposition.models import depositionModels

    def _make(**kwargs):
        kwargs.setdefault("surface", DESERT)
        kwargs.setdefault("diameter", 0.3)
        return depositionModels(projectName=unit_project.projectName, **kwargs)

    return _make


@pytest.mark.unit
class TestPropertyAccessors:
    def test_the_deposition_model_name_round_trips(self, make):
        model = make()
        model.depositionModel = "SomethingElse"
        assert model.depositionModel == "SomethingElse"

    def test_reading_ustar_goes_through_the_property(self, make):
        """B46: the accessor resolves to the diameter getter, so this is
        also the value the model computes with."""
        assert make(diameter=0.25).ustar == pytest.approx(0.25)

    def test_writing_ustar_goes_through_the_property(self, make):
        """B46: the write lands on `_diameter`."""
        model = make(diameter=0.25)
        model.ustar = 0.4
        assert model.diameter == pytest.approx(0.4)


@pytest.mark.unit
class TestPetroffCanopyBranch:
    def test_a_leaf_canopy_returns_a_positive_deposition_velocity(self, make):
        assert make(surface=FOREST).depositionRate_Petroff() > 0

    def test_a_canopy_deposits_faster_than_bare_ground(self, make):
        """Foliage adds collection surface, so the canopy velocity must
        exceed the bare-ground one for the same particle."""
        bare = make(surface=DESERT).depositionRate_Petroff()
        canopy = make(surface=FOREST).depositionRate_Petroff()
        assert canopy > bare

    def test_a_denser_canopy_deposits_faster(self, make):
        thin = make(surface=dict(FOREST, LAI=1.0)).depositionRate_Petroff()
        thick = make(surface=dict(FOREST, LAI=6.0)).depositionRate_Petroff()
        assert thick > thin

    def test_the_canopy_path_is_reached_only_via_obstacle_shape(self, make):
        """Dropping obstacleShape from an otherwise identical canopy
        surface falls back to the bare-ground formula, which proves the
        branch above is what produced the difference."""
        with_shape = make(surface=FOREST).depositionRate_Petroff()
        without = dict(FOREST)
        del without["obstacleShape"]
        assert make(surface=without).depositionRate_Petroff() != pytest.approx(with_shape)

    def test_an_obstacle_shape_without_a_full_efficiency_set_raises(self, make):
        """Only 'leaf' appears in all four efficiency maps; 'needle' is
        present in `eb` alone."""
        with pytest.raises(KeyError):
            make(surface=dict(FOREST, obstacleShape="needle")).depositionRate_Petroff()

    def test_an_out_of_range_stability_is_rejected_by_the_nested_functions(self, make):
        """PhiM/PsiM accept only -2 <= xi <= 1. Since
        xi = (hcanop - hdepl) / L and L is inversely proportional to the
        heat flux, a strongly convective surface pushes xi below -2, and
        the nested guards raise rather than returning nonsense."""
        with pytest.raises(KeyError, match="out of stability range"):
            make(surface=FOREST, heatFlux=5000.0).depositionRate_Petroff()
