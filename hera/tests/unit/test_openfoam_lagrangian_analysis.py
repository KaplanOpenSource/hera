"""analysis.calcConcentrationPointWise: binning particle mass into a
regular grid to get a concentration field (mass per cell volume)."""
import pandas
import pytest

from hera.simulations.openFoam.lagrangian.abstractLagrangianSolver import analysis


@pytest.fixture()
def calc():
    return analysis.__new__(analysis)


@pytest.mark.unit
class TestCalcConcentrationPointWise:
    def test_particles_in_the_same_cell_have_their_mass_summed(self, calc):
        df = pandas.DataFrame({
            "x": [0.1, 0.4], "y": [0.1, 0.4], "z": [0.0, 0.0],
            "mass": [1.0, 2.0], "time": [0, 0],
        })
        result = calc.calcConcentrationPointWise(df, dxdydz=0.5)
        assert len(result) == 1
        assert result["C"].iloc[0] == pytest.approx((1.0 + 2.0) / 0.5**3)

    def test_particles_in_different_cells_stay_separate(self, calc):
        df = pandas.DataFrame({
            "x": [0.1, 1.2], "y": [0.1, 0.1], "z": [0.0, 0.0],
            "mass": [1.0, 3.0], "time": [0, 0],
        })
        result = calc.calcConcentrationPointWise(df, dxdydz=0.5)
        assert len(result) == 2

    def test_concentration_scales_inversely_with_cell_volume(self, calc):
        df = pandas.DataFrame({"x": [0.0], "y": [0.0], "z": [0.0], "mass": [1.0], "time": [0]})
        coarse = calc.calcConcentrationPointWise(df, dxdydz=1.0)["C"].iloc[0]
        fine = calc.calcConcentrationPointWise(df, dxdydz=0.5)["C"].iloc[0]
        assert fine > coarse

    def test_different_timesteps_are_kept_as_separate_groups(self, calc):
        df = pandas.DataFrame({
            "x": [0.0, 0.0], "y": [0.0, 0.0], "z": [0.0, 0.0],
            "mass": [1.0, 5.0], "time": [0, 1],
        })
        result = calc.calcConcentrationPointWise(df, dxdydz=1.0)
        assert len(result) == 2
